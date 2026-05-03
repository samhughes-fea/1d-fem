# simulation_runner/harmonic/harmonic_simulation.py
"""§4 Frequency-domain harmonic (steady-state) analysis."""

from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional

import numpy as np
from scipy.sparse import coo_matrix

from processing.eigen.assembly import assemble_global_matrices
from processing.eigen.boundary_conditions import apply_boundary_conditions
from processing.harmonic.damping_table import interpolate_zeta_hz, load_zeta_vs_frequency_hz
from processing.harmonic.eigen_basis_io import load_modal_basis_from_modal_results_dir
from processing.harmonic.frequency_response import (
    frequency_grid_hz,
    geometric_mean_reference_omega_rad,
    harmonic_damping_matrix,
    sweep_displacements,
)
from processing.harmonic.modal_superposition import (
    harmonic_displacement_modal_superposition,
    undamped_natural_modes,
)
from processing.static.operations.assembly import AssembleGlobalSystem

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
    """Harmonic / frequency-response analysis using direct frequency sweep."""

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

    def _total_dof(self) -> int:
        num_nodes = len(self.mesh_dictionary["node_ids"])
        from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

        ed = self.settings.get("element_dictionary")
        dpn = 7 if ed is not None and mesh_uses_warping_dof(ed) else 6
        return int(num_nodes * dpn)

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

        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.secondary_results_dir, exist_ok=True)
        os.makedirs(self.tertiary_results_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)

        total_dof = self._total_dof()
        eo = np.asarray(self.settings.get("element_objects"), dtype=object).ravel()
        fo = np.asarray(self.settings.get("force_objects"), dtype=object).ravel()
        if eo.size == 0 or fo.size == 0:
            raise ValueError("Harmonic requires element_objects and force_objects (run_job wiring).")

        # Do not ravel pre-stacked ``element_stiffness_matrices`` numpy arrays: ``np.array([K0,..])``
        # can become numeric (n_el, n, n) and ravel to n_el*n*n scalars. Use formulation K_e per object.
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
        # Prefer numeric (n_el, ndof, ndof) stacks before raveling generic arrays (avoids 576 scalars).
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

        K_global, M_global, _ = assemble_global_matrices(
            elements=el_flat,
            element_stiffness_matrices=list(Ke_list),
            element_mass_matrices=list(Me_list),
            total_dof=total_dof,
            job_results_dir=os.path.join(self.primary_results_dir, "assembly_modal"),
        )

        assembler = AssembleGlobalSystem(
            elements=el_flat,
            element_stiffness_matrices=list(Ke_list),
            element_force_vectors=list(Fe_flat),
            total_dof=total_dof,
            job_results_dir=os.path.join(self.primary_results_dir, "assembly_force"),
        )
        _, F_global, _local_map, _assembly_map = assembler.assemble()
        F_global = np.asarray(F_global).ravel()

        K_mod, M_mod, bc_dofs = apply_boundary_conditions(
            K_global,
            M_global,
            prescribed_displacements=self.prescribed_displacement_dict,
        )
        F_mod = F_global.copy()
        if bc_dofs.size:
            bci = np.asarray(bc_dofs, dtype=np.int64)
            if np.iscomplexobj(F_mod):
                F_mod[bci] = 0.0 + 0.0j
            else:
                F_mod[bci] = 0.0

        f_ref_hz = float(np.sqrt(f_min * f_max))
        omega_ref = 2.0 * np.pi * f_ref_hz
        C_mod = harmonic_damping_matrix(M_mod, K_mod, zeta, omega_ref, ra, rb)

        F_use = F_mod.astype(np.complex128, copy=True)
        if load_phase != 0.0:
            F_use *= np.exp(1j * float(load_phase))

        u_template, dof_unknown = self._harmonic_prescribed_partition(
            bc_dofs,
            total_dof,
            global_prescribed_phase_rad=prescribed_phase,
        )

        freqs_hz = frequency_grid_hz(f_min, f_max, n_pts)

        ft_tab = zt_tab = None
        zeta_table_path = cfg.get("damping_zeta_table_file")
        if zeta_table_path:
            ft_tab, zt_tab = load_zeta_vs_frequency_hz(self._resolve_job_path(zeta_table_path))
        zpf = (
            interpolate_zeta_hz(freqs_hz, ft_tab, zt_tab) if ft_tab is not None else None
        )

        basis_dir = cfg.get("harmonic_modal_basis_dir")

        if use_modal:
            if dof_unknown is not None:
                raise ValueError(
                    "harmonic: use_modal_superposition is incompatible with non-zero prescribed "
                    "displacements in this version; use the direct frequency sweep instead."
                )
            if ra != 0.0 or rb != 0.0:
                logger.warning(
                    "harmonic: Rayleigh damping (rayleigh_alpha / rayleigh_beta) is not applied on "
                    "the modal superposition path; use the direct frequency sweep for full Rayleigh C."
                )
            if basis_dir:
                fh_hz, Phi_full = load_modal_basis_from_modal_results_dir(
                    self._resolve_job_path(basis_dir),
                    job_name=cfg.get("harmonic_modal_basis_job_name"),
                )
                nm = min(num_modal, int(Phi_full.shape[1]))
                Phi = np.asarray(Phi_full[:, :nm], dtype=np.float64)
                omega_n = (2.0 * np.pi * np.asarray(fh_hz[:nm], dtype=np.float64)).ravel()
            else:
                omega_n, Phi = undamped_natural_modes(K_mod, M_mod, num_modal)
            zeta_pm = None
            if ft_tab is not None:
                zeta_pm = interpolate_zeta_hz(omega_n / (2.0 * np.pi), ft_tab, zt_tab)
            U = harmonic_displacement_modal_superposition(
                omega_n,
                Phi,
                F_use,
                freqs_hz,
                zeta,
                zeta_per_mode=zeta_pm,
            )
        else:
            sweep_kw: Dict[str, Any] = dict(
                parallel=use_parallel,
                mp_damping_reference=mp_ref,
                zeta_mp=zeta,
                rayleigh_alpha_mp=ra,
                rayleigh_beta_mp=rb,
                linear_solver=lin_solver,
            )
            if zpf is not None:
                sweep_kw["zeta_per_frequency"] = zpf
                sweep_kw["geometric_omega_ref_rad"] = geometric_mean_reference_omega_rad(freqs_hz)
            U = sweep_displacements(
                K_mod,
                M_mod,
                C_mod,
                F_use,
                freqs_hz,
                u_prescribed=u_template,
                dof_unknown=dof_unknown,
                **sweep_kw,
            )

        out_dir = os.path.join(self.primary_results_dir, "harmonic_results")
        os.makedirs(out_dir, exist_ok=True)

        np.savetxt(os.path.join(out_dir, f"{self.job_name}_frequencies_hz.txt"), freqs_hz, fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_real.txt"), np.real(U), fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_imag.txt"), np.imag(U), fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_abs.txt"), np.abs(U), fmt="%.8e")
        np.savetxt(os.path.join(out_dir, f"{self.job_name}_displacement_phase_rad.txt"), np.angle(U), fmt="%.8e")

        self.primary_results["global"]["harmonic_frequencies_hz"] = freqs_hz
        self.primary_results["global"]["harmonic_displacement"] = U

        if self._harmonic_post_processing_enabled():
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

        logger.info(
            "Harmonic sweep completed: %s frequencies from %.4g to %.4g Hz -> %s",
            n_pts,
            f_min,
            f_max,
            out_dir,
        )
