# simulation_runner/harmonic/harmonic_simulation.py
"""Section 4 frequency-domain harmonic (steady-state) analysis."""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional

import numpy as np
from scipy.sparse import coo_matrix

from processing.common.primary_artifact_manifest import write_primary_artifact_manifest
from processing.harmonic.damping_table import interpolate_zeta_hz, load_zeta_vs_frequency_hz
from processing.harmonic.diagnostics.harmonic_run_diagnostic import log_harmonic_structural_summary
from processing.harmonic.frequency_response import frequency_grid_hz
from processing.boundary_supports import resolve_penalty_fixed_dofs
from processing.harmonic.operations import (
    AssembleHarmonicLoadVector,
    AssembleHarmonicStructuralMatrices,
    BuildHarmonicDampingMatrix,
    HarmonicSweepConfig,
    ModifyHarmonicStructuralMatrices,
    SolveHarmonicFrequencySweep,
)
from processing.static.diagnostics.runtime_monitor_telemetry import RuntimeMonitorTelemetry
from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

logger = logging.getLogger(__name__)

_DEFAULT_F_MIN_HZ = 1.0
_DEFAULT_F_MAX_HZ = 50.0
_DEFAULT_NUM_POINTS = 21
_DEFAULT_ZETA = 0.02


def effective_harmonic_config(simulation_settings: dict) -> Dict[str, Any]:
    """Resolve ``[Harmonic]`` keys with defaults (see docs/conventions/HARMONIC_FREQUENCY_DOMAIN.md)."""
    hm = simulation_settings.get("harmonic") or {}
    f_min = hm.get("frequency_min_hz")
    if f_min is None:
        f_min = _DEFAULT_F_MIN_HZ
    f_max = hm.get("frequency_max_hz")
    if f_max is None:
        f_max = _DEFAULT_F_MAX_HZ
    n_pts = hm.get("num_frequency_points")
    if n_pts is None:
        n_pts = _DEFAULT_NUM_POINTS
    zeta = hm.get("modal_damping_ratio")
    if zeta is None:
        zeta = _DEFAULT_ZETA
    ra = hm.get("rayleigh_alpha")
    rb = hm.get("rayleigh_beta")
    if ra is None:
        ra = 0.0
    if rb is None:
        rb = 0.0
    lp = hm.get("load_phase_rad")
    if lp is None:
        lp = 0.0
    pfs = hm.get("parallel_frequency_sweep")
    if pfs is None:
        pfs = False
    mp_ref = hm.get("mp_damping_reference")
    if mp_ref is None:
        mp_ref = "geometric_mean"
    mp_ref = str(mp_ref).strip().lower().replace("-", "_").replace(" ", "_")
    if mp_ref in ("geometricmean", "sweep_geometric_mean"):
        mp_ref = "geometric_mean"
    if mp_ref in ("currentfrequency", "at_frequency"):
        mp_ref = "current_frequency"
    use_modal = hm.get("use_modal_superposition")
    if use_modal is None:
        use_modal = False
    ns_modal = hm.get("modal_superposition_num_modes")
    if ns_modal is None:
        ns_modal = 10
    prescribed_phase = hm.get("prescribed_motion_phase_rad")
    if prescribed_phase is None:
        prescribed_phase = 0.0
    lin_sol = hm.get("harmonic_linear_solver")
    if lin_sol is None:
        lin_sol = "spsolve"
    lin_sol = str(lin_sol).strip().lower()
    f_min = float(f_min)
    f_max = float(f_max)
    n_pts = int(n_pts)
    zeta = float(zeta)
    ra = float(ra)
    rb = float(rb)
    lp = float(lp)
    pfs = bool(pfs)
    if f_max <= f_min:
        raise ValueError(
            f"harmonic: frequency_max_hz ({f_max}) must be greater than frequency_min_hz ({f_min})"
        )
    if n_pts < 2:
        raise ValueError("harmonic: num_frequency_points must be >= 2")
    if zeta < 0.0:
        raise ValueError("harmonic: modal_damping_ratio must be non-negative")
    if ra < 0.0 or rb < 0.0:
        raise ValueError("harmonic: rayleigh_alpha and rayleigh_beta must be non-negative")
    if mp_ref not in ("geometric_mean", "current_frequency"):
        raise ValueError(
            "harmonic: mp_damping_reference must be 'geometric_mean' or 'current_frequency'"
        )
    if lin_sol not in ("spsolve", "splu"):
        raise ValueError("harmonic: harmonic_linear_solver must be 'spsolve' or 'splu'")

    zeta_file = hm.get("damping_zeta_table_file")
    basis_dir = hm.get("harmonic_modal_basis_dir")
    basis_job = hm.get("harmonic_modal_basis_job_name")
    if zeta_file is not None and not isinstance(zeta_file, str):
        zeta_file = str(zeta_file)
    if basis_dir is not None and not isinstance(basis_dir, str):
        basis_dir = str(basis_dir)
    if basis_job is not None and not isinstance(basis_job, str):
        basis_job = str(basis_job) if basis_job else None
    if basis_job == "":
        basis_job = None

    return {
        "frequency_min_hz": f_min,
        "frequency_max_hz": f_max,
        "num_frequency_points": n_pts,
        "modal_damping_ratio": zeta,
        "rayleigh_alpha": ra,
        "rayleigh_beta": rb,
        "load_phase_rad": lp,
        "parallel_frequency_sweep": pfs,
        "mp_damping_reference": mp_ref,
        "use_modal_superposition": bool(use_modal),
        "modal_superposition_num_modes": int(ns_modal),
        "prescribed_motion_phase_rad": float(prescribed_phase),
        "harmonic_linear_solver": lin_sol,
        "damping_zeta_table_file": zeta_file,
        "harmonic_modal_basis_dir": basis_dir,
        "harmonic_modal_basis_job_name": basis_job,
    }


class HarmonicSimulationRunner:
    """
    Harmonic / frequency-response analysis using direct frequency sweep.

    Staged methods (``collect_harmonic_assembly_inputs``, ``assemble_harmonic_structural_matrices``, …)
    mirror :class:`~simulation_runner.static.linear_static_simulation.LinearStaticSimulationRunner` orchestration.
    """

    def __init__(self, settings: dict, job_name: str):
        self.settings = settings
        self.job_name = job_name
        self.simulation_settings = settings.get("simulation_settings") or {}

        raw_el = self.settings.get("elements", np.array([]))
        self.elements = (
            raw_el if isinstance(raw_el, np.ndarray) else np.asarray(raw_el, dtype=object)
        )
        self.mesh_dictionary = self.settings.get("mesh_dictionary") or {}
        if self.elements.size == 0 or not self.mesh_dictionary:
            raise ValueError("Harmonic simulation requires elements and mesh_dictionary in settings")

        self.job_results_dir = self.settings.get("job_results_dir")
        if not self.job_results_dir:
            raise ValueError("Harmonic simulation requires job_results_dir")

        self.primary_results_dir = os.path.join(self.job_results_dir, "primary_results")
        self.secondary_results_dir = os.path.join(self.job_results_dir, "secondary_results")
        self.tertiary_results_dir = os.path.join(self.job_results_dir, "tertiary_results")
        self.diagnostics_dir = os.path.join(self.job_results_dir, "diagnostics")
        self.logs_dir = os.path.join(self.job_results_dir, "logs")

        self.primary_results: dict = {"global": {}, "element": {"data": []}}
        self.secondary_results: dict = {"global": {}, "element": {"data": []}}

        self.results_root = self.job_results_dir
        self.prescribed_displacement_dict = self.settings.get("prescribed_displacement_dict")
        self.grid_dictionary = self.settings.get("grid_dictionary") or {}
        self.element_dictionary = self.settings.get("element_dictionary") or {}
        self.material_dictionary = self.settings.get("material_dictionary") or {}
        self.section_dictionary = self.settings.get("section_dictionary") or {}
        self.job_dir = self.settings.get("job_dir")

        self.K_global = None
        self.M_global = None
        self.F_global = None
        self.K_mod = None
        self.M_mod = None
        self.C_mod = None
        self.harmonic_frequencies_hz = None
        self.harmonic_displacement = None

    def _dof_per_node(self) -> int:
        ed = self.settings.get("element_dictionary")
        return 7 if ed is not None and mesh_uses_warping_dof(ed) else 6

    def _total_dof(self) -> int:
        num_nodes = len(self.mesh_dictionary["node_ids"])
        return int(num_nodes * self._dof_per_node())

    def _resolved_penalty_fixed_dofs(self, total_dof: int):
        hm = self.simulation_settings.get("harmonic") or {}
        return resolve_penalty_fixed_dofs(
            total_dof=total_dof,
            dof_per_node=self._dof_per_node(),
            prescribed_displacement_dict=self.prescribed_displacement_dict,
            section_settings=hm,
            grid_node_ids=self.mesh_dictionary.get("node_ids"),
        )

    def _ensure_coo_ke(self, matrices) -> np.ndarray:
        out = []
        for matrix in matrices:
            if isinstance(matrix, coo_matrix):
                out.append(matrix)
            else:
                dense = np.asarray(matrix, dtype=np.float64)
                out.append(coo_matrix(dense))
        return np.array(out, dtype=object)

    def _harmonic_prescribed_partition(
        self,
        bc_dofs: np.ndarray,
        n_dof: int,
        *,
        global_prescribed_phase_rad: float = 0.0,
    ) -> tuple[np.ndarray | None, np.ndarray | None]:
        """Return (u_prescribed_template, dof_unknown) for partitioned harmonic solve, or (None, None)."""
        pd = self.prescribed_displacement_dict
        if pd is None or len(pd.get("global_dof", [])) == 0:
            return None, None
        gd = np.asarray(pd["global_dof"], dtype=np.int64)
        val = np.asarray(pd["value"], dtype=np.float64)
        motion_mask = np.abs(val) >= 1e-12
        motion_dofs = gd[motion_mask]
        if motion_dofs.size == 0:
            return None, None
        u_p = np.zeros(n_dof, dtype=np.complex128)
        ph_rows = np.asarray(pd.get("phase_rad", np.zeros(len(val))), dtype=np.float64).ravel()
        if ph_rows.shape[0] != val.shape[0]:
            ph_rows = np.zeros_like(val, dtype=np.float64)
        ph_eff = ph_rows[motion_mask] + float(global_prescribed_phase_rad)
        u_p[motion_dofs] = val[motion_mask] * np.exp(1j * ph_eff)
        bc = np.asarray(bc_dofs, dtype=np.int64)
        unknown = np.setdiff1d(
            np.arange(n_dof, dtype=np.int64),
            np.union1d(bc, motion_dofs),
        )
        if unknown.size == 0:
            raise ValueError(
                "harmonic: no free DOFs left after prescribed displacements; check boundary conditions"
            )
        return u_p, unknown

    def _resolve_job_path(self, p: str) -> str:
        if os.path.isabs(p):
            return p
        base = self.job_dir
        if not base:
            raise ValueError("harmonic: relative path requires job_dir in settings (run_job wiring)")
        return os.path.normpath(os.path.join(base, p))

    def _run_secondary_tertiary_from_formulation_cache(
        self,
        U_global: np.ndarray,
        *,
        results_subdir: str | None = None,
    ) -> None:
        from processing.static.results.containers.formulation_results import (
            FormulationResultSet,
            strict_shape_functions_validation_from_env,
            validate_shape_functions_populated,
        )
        from processing.static.results.postprocess_secondary_tertiary import (
            run_secondary_tertiary_from_formulation_cache,
        )

        eo = self.settings.get("element_objects")
        fo = self.settings.get("force_objects")
        if eo is None or fo is None:
            raise ValueError(
                "post_processing.run_secondary_tertiary_harmonic requires element_objects and force_objects"
            )
        cache = FormulationResultSet(
            element_objects=list(np.asarray(eo, dtype=object).ravel()),
            force_objects=list(np.asarray(fo, dtype=object).ravel()),
        )
        validate_shape_functions_populated(
            cache.element_objects,
            cache.force_objects,
            strict=strict_shape_functions_validation_from_env(),
        )
        U = np.asarray(U_global, dtype=np.float64).ravel()
        el_flat = list(np.asarray(self.elements, dtype=object).ravel())
        sec_dir = self.secondary_results_dir
        ter_dir = self.tertiary_results_dir
        if results_subdir:
            sec_dir = os.path.join(sec_dir, results_subdir)
            ter_dir = os.path.join(ter_dir, results_subdir)
            os.makedirs(sec_dir, exist_ok=True)
            os.makedirs(ter_dir, exist_ok=True)
        run_secondary_tertiary_from_formulation_cache(
            elements=el_flat,
            grid_dictionary=self.grid_dictionary,
            element_dictionary=self.element_dictionary,
            material_dictionary=self.material_dictionary,
            section_dictionary=self.section_dictionary,
            U_global=U,
            formulation_cache=cache,
            results_root=self.results_root,
            secondary_results_dir=sec_dir,
            tertiary_results_dir=ter_dir,
        )

    def _harmonic_post_processing_enabled(self) -> bool:
        pp = self.simulation_settings.get("post_processing") or {}
        return bool(pp.get("run_secondary_tertiary_harmonic"))

    def prepare_harmonic_job_tree(self) -> None:
        """Create primary/secondary/tertiary/diagnostics/logs folders."""
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.secondary_results_dir, exist_ok=True)
        os.makedirs(self.tertiary_results_dir, exist_ok=True)
        os.makedirs(self.diagnostics_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)

    def collect_harmonic_assembly_inputs(self):
        """
        Build element lists and COO matrices for structural and load assembly.

        Returns
        -------
        el_flat, Ke_list, Me_list, Fe_flat, total_dof, modal_assembly_dir, force_assembly_dir
        """
        total_dof = self._total_dof()
        eo = np.asarray(self.settings.get("element_objects"), dtype=object).ravel()
        fo = np.asarray(self.settings.get("force_objects"), dtype=object).ravel()
        if eo.size == 0 or fo.size == 0:
            raise ValueError("Harmonic requires element_objects and force_objects (run_job wiring).")

        Ke_list = self._ensure_coo_ke([obj.K_e for obj in eo])
        Fe_flat = np.empty(fo.shape[0], dtype=object)
        for i, obj in enumerate(fo):
            Fe_flat[i] = np.asarray(obj.F_e).ravel()

        el_flat = list(np.asarray(self.elements, dtype=object).ravel())
        n_el = len(el_flat)

        Me_src = self.settings.get("element_mass_matrices")
        if Me_src is None:
            raise ValueError("Harmonic requires element_mass_matrices in settings")
        Me_arr = np.asarray(Me_src)
        if Me_arr.ndim == 3 and Me_arr.shape[0] == n_el:
            Me_flat = [Me_arr[i] for i in range(n_el)]
        elif Me_arr.dtype == object and Me_arr.size == n_el:
            Me_flat = list(Me_arr.ravel())
        elif Me_arr.dtype == object and Me_arr.ndim == 1 and len(Me_arr) == n_el:
            Me_flat = list(Me_arr)
        else:
            raise ValueError(
                f"Harmonic: cannot interpret element_mass_matrices with shape {Me_arr.shape}, dtype={Me_arr.dtype}"
            )
        if len(Me_flat) != n_el:
            raise ValueError(
                f"Harmonic: mass matrix count {len(Me_flat)} != element count {n_el}"
            )
        Me_list = self._ensure_coo_ke(Me_flat)

        if len(el_flat) != len(Ke_list) or len(el_flat) != len(Me_list):
            raise ValueError(
                "Harmonic assembly count mismatch: "
                f"elements={len(el_flat)}, K_e={len(Ke_list)}, M_e={len(Me_list)}"
            )

        modal_assembly_dir = os.path.join(self.primary_results_dir, "assembly_modal")
        force_assembly_dir = os.path.join(self.primary_results_dir, "assembly_force")
        return el_flat, Ke_list, Me_list, Fe_flat, total_dof, modal_assembly_dir, force_assembly_dir

    def assemble_harmonic_structural_matrices(
        self, el_flat, Ke_list, Me_list, total_dof: int, modal_assembly_dir: str
    ):
        """Global ``K``, ``M`` for harmonic analysis."""
        K_global, M_global, _ = AssembleHarmonicStructuralMatrices(
            elements=el_flat,
            element_stiffness_matrices=list(Ke_list),
            element_mass_matrices=list(Me_list),
            total_dof=total_dof,
            job_results_dir=modal_assembly_dir,
        ).run()
        self.K_global = K_global
        self.M_global = M_global
        return K_global, M_global

    def assemble_harmonic_load_vector(
        self, el_flat, Ke_list, Fe_flat, total_dof: int, force_assembly_dir: str
    ):
        """Global load vector ``F``."""
        F_global = AssembleHarmonicLoadVector(
            elements=el_flat,
            element_stiffness_matrices=list(Ke_list),
            element_force_vectors=list(Fe_flat),
            total_dof=total_dof,
            job_results_dir=force_assembly_dir,
        ).run()
        self.F_global = F_global
        return F_global

    def modify_harmonic_structural_matrices(self, K_global, M_global, F_global, total_dof: int):
        """BCs on structural matrices and load."""
        h_fixed = self._resolved_penalty_fixed_dofs(total_dof)
        K_mod, M_mod, bc_dofs, F_mod = ModifyHarmonicStructuralMatrices(
            prescribed_displacements=self.prescribed_displacement_dict,
            job_results_dir=self.primary_results_dir,
            fixed_dofs=h_fixed,
        ).run(K_global, M_global, F_global)
        self.K_mod = K_mod
        self.M_mod = M_mod
        return K_mod, M_mod, bc_dofs, F_mod

    def build_harmonic_damping_matrix(self, M_mod, K_mod, zeta: float, omega_ref: float, ra: float, rb: float):
        """Damping matrix ``C`` for the sweep."""
        C_mod = BuildHarmonicDampingMatrix(
            job_results_dir=self.primary_results_dir,
        ).run(M_mod, K_mod, zeta, omega_ref, ra, rb)
        self.C_mod = C_mod
        return C_mod

    def log_harmonic_structural_diagnostics(self, K_mod, M_mod, C_mod, F_mod, bc_dofs) -> None:
        log_harmonic_structural_summary(
            K_mod,
            M_mod,
            C_mod,
            F_mod,
            n_bc_dofs=int(np.asarray(bc_dofs).size),
            job_results_dir=self.primary_results_dir,
        )

    def build_harmonic_complex_load_vector(self, F_mod, load_phase_rad: float) -> np.ndarray:
        F_use = F_mod.astype(np.complex128, copy=True)
        if load_phase_rad != 0.0:
            F_use *= np.exp(1j * float(load_phase_rad))
        return F_use

    def solve_harmonic_frequency_sweep(
        self,
        K_mod,
        M_mod,
        C_mod,
        F_use: np.ndarray,
        freqs_hz: np.ndarray,
        u_template,
        dof_unknown,
        sweep_cfg: HarmonicSweepConfig,
    ) -> np.ndarray:
        U = SolveHarmonicFrequencySweep(
            job_results_dir=self.primary_results_dir,
            resolve_job_path=self._resolve_job_path,
        ).run(
            K_mod,
            M_mod,
            C_mod,
            F_use,
            freqs_hz,
            u_template,
            dof_unknown,
            sweep_cfg,
        )
        self.harmonic_displacement = U
        return U

    def save_harmonic_primary_results(self, freqs_hz: np.ndarray, U: np.ndarray) -> str:
        """Write ``harmonic_results/*.txt``; returns output directory path."""
        out_dir = os.path.join(self.primary_results_dir, "harmonic_results")
        os.makedirs(out_dir, exist_ok=True)
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_frequencies_hz.txt"), freqs_hz, fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_real.txt"), np.real(U), fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_imag.txt"), np.imag(U), fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_abs.txt"), np.abs(U), fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_phase_rad.txt"), np.angle(U), fmt="%.8e")
        return out_dir

    def write_harmonic_primary_artifact_manifest(self) -> None:
        hr = "primary_results/harmonic_results"
        jn = self.job_name
        write_primary_artifact_manifest(
            self.results_root,
            family="harmonic",
            job_name=jn,
            artifacts={
                "frequencies_hz": f"{hr}/{jn}_frequencies_hz.txt",
                "displacement_real": f"{hr}/{jn}_displacement_real.txt",
                "displacement_imag": f"{hr}/{jn}_displacement_imag.txt",
                "displacement_abs": f"{hr}/{jn}_displacement_abs.txt",
                "displacement_phase_rad": f"{hr}/{jn}_displacement_phase_rad.txt",
            },
        )

    def run_harmonic_post_processing_if_enabled(self, U: np.ndarray) -> None:
        if not self._harmonic_post_processing_enabled():
            return
        pp = self.simulation_settings.get("post_processing") or {}
        n_col = int(U.shape[1])
        raw_multi = pp.get("harmonic_secondary_tertiary_frequency_indices")
        indices: List[int] = []
        if raw_multi is not None and raw_multi != "":
            if isinstance(raw_multi, str):
                indices = [int(x.strip()) for x in raw_multi.split(",") if x.strip()]
            else:
                indices = [int(x) for x in raw_multi]
        export_all = bool(pp.get("harmonic_secondary_tertiary_all_frequencies", False))
        comp_raw = pp.get("harmonic_secondary_tertiary_displacement_component")
        comp = "real" if comp_raw is None else str(comp_raw).strip().lower()
        if comp not in ("real", "imag", "both"):
            raise ValueError(
                "post_processing.harmonic_secondary_tertiary_displacement_component must be "
                f"'real', 'imag', or 'both', got {comp_raw!r}"
            )

        def run_one_frequency_column(idx: int) -> None:
            Ucol = np.asarray(U[:, idx], dtype=np.complex128)
            if comp in ("real", "both"):
                self._run_secondary_tertiary_from_formulation_cache(
                    np.real(Ucol),
                    results_subdir=os.path.join("harmonic_post", f"freq_{idx:04d}"),
                )
            if comp in ("imag", "both"):
                self._run_secondary_tertiary_from_formulation_cache(
                    np.imag(Ucol),
                    results_subdir=os.path.join("harmonic_post", f"freq_{idx:04d}_imag"),
                )

        if indices:
            for idx in indices:
                if idx < 0 or idx >= n_col:
                    raise ValueError(
                        f"post_processing.harmonic_secondary_tertiary_frequency_indices "
                        f"entry {idx} out of range for {n_col} samples"
                    )
                run_one_frequency_column(idx)
        elif export_all:
            for idx in range(n_col):
                run_one_frequency_column(idx)
        else:
            idx = int(pp.get("harmonic_frequency_index", 0))
            if idx < 0 or idx >= n_col:
                raise ValueError(
                    f"post_processing.harmonic_frequency_index {idx} out of range "
                    f"for {n_col} samples"
                )
            run_one_frequency_column(idx)

    def run(self) -> None:
        cfg = effective_harmonic_config(self.simulation_settings)
        f_min = cfg["frequency_min_hz"]
        f_max = cfg["frequency_max_hz"]
        n_pts = cfg["num_frequency_points"]
        zeta = cfg["modal_damping_ratio"]
        ra = cfg["rayleigh_alpha"]
        rb = cfg["rayleigh_beta"]
        load_phase = cfg["load_phase_rad"]
        use_parallel = cfg["parallel_frequency_sweep"]
        mp_ref = cfg["mp_damping_reference"]
        use_modal = cfg["use_modal_superposition"]
        num_modal = cfg["modal_superposition_num_modes"]
        prescribed_phase = cfg["prescribed_motion_phase_rad"]
        lin_solver = cfg["harmonic_linear_solver"]

        self.prepare_harmonic_job_tree()
        monitor = RuntimeMonitorTelemetry(job_results_dir=self.diagnostics_dir)

        el_flat, Ke_list, Me_list, Fe_flat, total_dof, modal_assembly_dir, force_assembly_dir = (
            self.collect_harmonic_assembly_inputs()
        )

        with monitor.stage("AssembleHarmonicStructuralMatrices"):
            K_global, M_global = self.assemble_harmonic_structural_matrices(
                el_flat, Ke_list, Me_list, total_dof, modal_assembly_dir
            )

        with monitor.stage("AssembleHarmonicLoadVector"):
            F_global = self.assemble_harmonic_load_vector(
                el_flat, Ke_list, Fe_flat, total_dof, force_assembly_dir
            )

        with monitor.stage("ModifyHarmonicStructuralMatrices"):
            K_mod, M_mod, bc_dofs, F_mod = self.modify_harmonic_structural_matrices(
                K_global, M_global, F_global, total_dof
            )

        f_ref_hz = float(np.sqrt(f_min * f_max))
        omega_ref = 2.0 * np.pi * f_ref_hz
        with monitor.stage("BuildHarmonicDampingMatrix"):
            C_mod = self.build_harmonic_damping_matrix(M_mod, K_mod, zeta, omega_ref, ra, rb)

        self.log_harmonic_structural_diagnostics(K_mod, M_mod, C_mod, F_mod, bc_dofs)

        F_use = self.build_harmonic_complex_load_vector(F_mod, load_phase)

        u_template, dof_unknown = self._harmonic_prescribed_partition(
            bc_dofs,
            total_dof,
            global_prescribed_phase_rad=prescribed_phase,
        )

        freqs_hz = frequency_grid_hz(f_min, f_max, n_pts)
        self.harmonic_frequencies_hz = freqs_hz

        ft_tab = zt_tab = None
        zeta_table_path = cfg.get("damping_zeta_table_file")
        if zeta_table_path:
            ft_tab, zt_tab = load_zeta_vs_frequency_hz(self._resolve_job_path(zeta_table_path))
        zpf = (
            interpolate_zeta_hz(freqs_hz, ft_tab, zt_tab) if ft_tab is not None else None
        )

        basis_dir = cfg.get("harmonic_modal_basis_dir")
        sweep_cfg = HarmonicSweepConfig(
            use_modal=use_modal,
            num_modal=num_modal,
            basis_dir=basis_dir,
            basis_job=cfg.get("harmonic_modal_basis_job_name"),
            zeta=zeta,
            use_parallel=use_parallel,
            mp_ref=mp_ref,
            rayleigh_alpha=ra,
            rayleigh_beta=rb,
            lin_solver=lin_solver,
            zpf=zpf,
            ft_tab=ft_tab,
            zt_tab=zt_tab,
        )
        with monitor.stage("SolveHarmonicFrequencySweep"):
            U = self.solve_harmonic_frequency_sweep(
                K_mod,
                M_mod,
                C_mod,
                F_use,
                freqs_hz,
                u_template,
                dof_unknown,
                sweep_cfg,
            )

        out_dir = self.save_harmonic_primary_results(freqs_hz, U)
        self.write_harmonic_primary_artifact_manifest()

        self.primary_results["global"]["harmonic_frequencies_hz"] = freqs_hz
        self.primary_results["global"]["harmonic_displacement"] = U

        self.run_harmonic_post_processing_if_enabled(U)

        logger.info(
            "Harmonic sweep completed: %s frequencies from %.4g to %.4g Hz -> %s",
            n_pts,
            f_min,
            f_max,
            out_dir,
        )
