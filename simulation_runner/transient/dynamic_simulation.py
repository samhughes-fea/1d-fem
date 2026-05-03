# simulation_runner/transient/dynamic_simulation.py
"""Transient dynamics (§3): time integration M u'' + C u' + K u = F(t)."""

from __future__ import annotations

import logging
import os
import warnings
from typing import Any, Callable, Optional

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from processing.boundary_supports import resolve_penalty_fixed_dofs
from processing.common.primary_artifact_manifest import write_primary_artifact_manifest
from processing.dynamic.diagnostics.transient_run_diagnostic import log_transient_modified_system
from processing.dynamic.operations import (
    AssembleDynamicGlobalSystem,
    IntegrateTransientSystem,
    ModifyDynamicGlobalSystem,
)
from processing.dynamic.transient_forcing import build_transient_force_func
from processing.dynamic.assembly import assemble_global_force_vector
from processing.static.diagnostics.runtime_monitor_telemetry import RuntimeMonitorTelemetry
from pre_processing.element_library.beam_warping import mesh_uses_warping_dof
from pre_processing.parsing.simulation_settings_resolution import effective_transient_config

logger = logging.getLogger(__name__)


class TransientSimulationRunner:
    """
    Transient finite element analysis: assemble K/M/C, apply BCs, Newmark integrate.

    Instance methods ``assemble_dynamic_global_system``, ``modify_dynamic_global_system``,
    ``integrate_transient_system``, etc. mirror the staged style of
    :class:`~simulation_runner.static.linear_static_simulation.LinearStaticSimulationRunner`.
    """

    def __init__(self, settings: dict, job_name: str):
        self.settings = settings
        self.job_name = job_name
        self.elements = self.settings.get("elements", np.array([]))
        self.mesh_dictionary = self.settings.get("mesh_dictionary", {})
        self.grid_dictionary = self.settings.get("grid_dictionary") or self.mesh_dictionary
        if self.elements.size == 0 or not self.mesh_dictionary:
            raise ValueError("Missing elements or mesh data in settings")

        self.element_dictionary = self.settings.get("element_dictionary")
        self.material_dictionary = self.settings.get("material_dictionary")
        self.section_dictionary = self.settings.get("section_dictionary")
        self.point_load_array = self.settings.get("point_load_array", np.empty((0, 9)))
        self.distributed_load_array = self.settings.get("distributed_load_array", np.empty((0, 9)))
        self.element_objects = self.settings.get("element_objects")
        self.force_objects = self.settings.get("force_objects")
        self.prescribed_displacement_dict = self.settings.get("prescribed_displacement_dict")

        self.element_stiffness_matrices = self._ensure_sparse_format(
            self.settings.get("element_stiffness_matrices", None)
        )
        self.element_mass_matrices = self._ensure_sparse_format(
            self.settings.get("element_mass_matrices", None)
        )
        self.element_damping_matrices = self.settings.get("element_damping_matrices")

        job_results_dir = self.settings.get("job_results_dir")
        if job_results_dir:
            self.results_root = job_results_dir
            self.primary_results_dir = os.path.join(job_results_dir, "primary_results")
            self.secondary_results_dir = os.path.join(job_results_dir, "secondary_results")
            self.tertiary_results_dir = os.path.join(job_results_dir, "tertiary_results")
            self.diagnostics_dir = os.path.join(job_results_dir, "diagnostics")
            self.logs_dir = os.path.join(job_results_dir, "logs")
        else:
            self.results_root = os.path.join("post_processing", "results", self.job_name)
            self.primary_results_dir = os.path.join(self.results_root, "primary_results")
            self.secondary_results_dir = os.path.join(self.results_root, "secondary_results")
            self.tertiary_results_dir = os.path.join(self.results_root, "tertiary_results")
            self.diagnostics_dir = os.path.join(self.results_root, "diagnostics")
            self.logs_dir = os.path.join(self.results_root, "logs")

        self.simulation_settings = self.settings.get("simulation_settings", {})

        self.K_global: csr_matrix | None = None
        self.M_global: csr_matrix | None = None
        self.C_global: csr_matrix | None = None
        self.K_mod: csr_matrix | None = None
        self.M_mod: csr_matrix | None = None
        self.C_mod: csr_matrix | None = None
        self.U_time_history: np.ndarray | None = None
        self.V_time_history: np.ndarray | None = None
        self.A_time_history: np.ndarray | None = None
        self.t_grid: np.ndarray | None = None

    def _ensure_sparse_format(self, matrices):
        if matrices is None:
            return None
        return np.array(
            [coo_matrix(m) if not isinstance(m, coo_matrix) else m for m in matrices],
            dtype=object,
        )

    def _dynamic_post_enabled(self) -> bool:
        cfg = self.simulation_settings.get("post_processing") or {}
        return bool(cfg.get("run_secondary_tertiary_dynamic", False))

    def _resolve_time_index(self, n_rows: int) -> int:
        cfg = self.simulation_settings.get("post_processing") or {}
        ti = int(cfg.get("dynamic_time_index", -1))
        if ti < 0:
            ti = n_rows + ti
        if ti < 0 or ti >= n_rows:
            raise IndexError(
                f"dynamic_time_index resolves to {ti} but displacement history has {n_rows} row(s)"
            )
        return ti

    def _dynamic_snapshot_row_indices(self, n_rows: int) -> list[int]:
        """Rows of U(t) for optional formulation-cache post (single or multi-snapshot)."""
        cfg = self.simulation_settings.get("post_processing") or {}
        raw = cfg.get("dynamic_time_indices")
        if raw is not None and not (isinstance(raw, (list, tuple)) and len(raw) == 0):
            if isinstance(raw, str) and not raw.strip():
                return [self._resolve_time_index(n_rows)]
            parts: list[int]
            if isinstance(raw, (list, tuple)):
                parts = [int(x) for x in raw]
            else:
                parts = [int(x.strip()) for x in str(raw).split(",") if str(x).strip()]
            out: list[int] = []
            for ti in parts:
                ri = int(ti)
                if ri < 0:
                    ri = n_rows + ri
                if ri < 0 or ri >= n_rows:
                    raise IndexError(
                        f"dynamic_time_indices entry resolves to {ri} but displacement "
                        f"history has {n_rows} row(s)"
                    )
                out.append(ri)
            return out
        return [self._resolve_time_index(n_rows)]

    def _run_secondary_tertiary_from_cache(
        self, U_global: np.ndarray, results_subdir: str | None = None
    ) -> None:
        from processing.static.results.containers.formulation_results import (
            FormulationResultSet,
            strict_shape_functions_validation_from_env,
            validate_shape_functions_populated,
        )
        from processing.static.results.postprocess_secondary_tertiary import (
            run_secondary_tertiary_from_formulation_cache,
        )

        required = {
            "grid_dictionary": self.grid_dictionary,
            "element_dictionary": self.element_dictionary,
            "material_dictionary": self.material_dictionary,
            "section_dictionary": self.section_dictionary,
        }
        missing = [k for k, v in required.items() if v is None]
        if missing:
            raise ValueError(
                "post_processing.run_secondary_tertiary_dynamic requires job wiring for: "
                + ", ".join(missing)
            )
        if self.element_objects is None or self.force_objects is None:
            raise ValueError(
                "post_processing.run_secondary_tertiary_dynamic requires element_objects and force_objects"
            )

        cache = FormulationResultSet(
            element_objects=list(np.asarray(self.element_objects, dtype=object).ravel()),
            force_objects=list(np.asarray(self.force_objects, dtype=object).ravel()),
        )
        validate_shape_functions_populated(
            cache.element_objects,
            cache.force_objects,
            strict=strict_shape_functions_validation_from_env(),
        )
        U = np.asarray(U_global, dtype=np.float64).ravel()
        sec_dir = self.secondary_results_dir
        ter_dir = self.tertiary_results_dir
        if results_subdir:
            sec_dir = os.path.join(sec_dir, results_subdir)
            ter_dir = os.path.join(ter_dir, results_subdir)
            os.makedirs(sec_dir, exist_ok=True)
            os.makedirs(ter_dir, exist_ok=True)
        run_secondary_tertiary_from_formulation_cache(
            elements=list(self.elements),
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

    def _dof_per_node(self) -> int:
        ed = self.element_dictionary
        return 7 if ed is not None and mesh_uses_warping_dof(ed) else 6

    def _mesh_node_ids(self):
        node_ids = self.mesh_dictionary.get("node_ids")
        if node_ids is None and self.grid_dictionary is not None:
            node_ids = self.grid_dictionary.get("ids", [])
        return node_ids

    def _total_dof(self) -> int:
        node_ids = self._mesh_node_ids()
        if node_ids is None:
            num_nodes = 0
        else:
            num_nodes = len(np.asarray(node_ids, dtype=object).ravel())
        return int(num_nodes * self._dof_per_node())

    def _resolved_penalty_fixed_dofs(self):
        dyn = effective_transient_config(self.simulation_settings)
        return resolve_penalty_fixed_dofs(
            total_dof=self._total_dof(),
            dof_per_node=self._dof_per_node(),
            prescribed_displacement_dict=self.prescribed_displacement_dict,
            section_settings=dyn,
            grid_node_ids=self._mesh_node_ids(),
        )

    def _assemble_reference_force(self, total_dof: int) -> np.ndarray:
        """Constant-in-time equivalent nodal load from element force objects (same as static assembly)."""
        eo = self.element_objects
        fo = self.force_objects
        if eo is None or fo is None:
            return np.zeros(total_dof, dtype=np.float64)
        if self.element_dictionary is None or self.grid_dictionary is None:
            logger.warning(
                "Transient: missing grid_dictionary or element_dictionary; using zero external force"
            )
            return np.zeros(total_dof, dtype=np.float64)
        el_list = list(np.asarray(self.elements, dtype=object).ravel())
        fos = list(np.asarray(fo, dtype=object).ravel())
        Fe_list = [np.asarray(obj.F_e).ravel() for obj in fos]
        return assemble_global_force_vector(
            el_list,
            Fe_list,
            total_dof,
            job_results_dir=os.path.join(self.primary_results_dir, "transient_force_assembly"),
        )

    def setup_simulation(self):
        """Create output directory structure."""
        logger.info("Setting up dynamic simulation for job: %s", self.job_name)
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.secondary_results_dir, exist_ok=True)
        os.makedirs(self.tertiary_results_dir, exist_ok=True)
        os.makedirs(self.diagnostics_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        logger.info("Results will be saved to: %s", self.results_root)

    def assemble_dynamic_global_system(self, total_dof: int) -> tuple[csr_matrix, csr_matrix, Any, Any]:
        """Global ``K``, ``M``, optional ``C``, and assembly metadata. Caches ``self.K_global``, ``self.M_global``, ``self.C_global``."""
        K_global, M_global, C_global, meta = AssembleDynamicGlobalSystem(
            elements=list(self.elements),
            element_stiffness_matrices=self.element_stiffness_matrices,
            element_mass_matrices=self.element_mass_matrices,
            total_dof=total_dof,
            job_results_dir=self.primary_results_dir,
            element_damping_matrices=self.element_damping_matrices,
        ).run()
        self.K_global = K_global
        self.M_global = M_global
        self.C_global = C_global
        return K_global, M_global, C_global, meta

    def apply_rayleigh_damping_if_needed(
        self,
        K_global: csr_matrix,
        M_global: csr_matrix,
        C_global: Any,
        dyn_config: dict,
    ) -> Any:
        """Return effective ``C`` (Rayleigh fill-in when element ``C`` is empty)."""
        ra = dyn_config.get("rayleigh_alpha")
        rb = dyn_config.get("rayleigh_beta")
        ra_f = float(ra) if ra is not None else 0.0
        rb_f = float(rb) if rb is not None else 0.0
        if (C_global is None or getattr(C_global, "nnz", 0) == 0) and (ra_f != 0.0 or rb_f != 0.0):
            C_use = (ra_f * M_global + rb_f * K_global).tocsr()
            logger.info(
                "Transient Rayleigh damping assembled: alpha=%s beta=%s (element C absent)",
                ra_f,
                rb_f,
            )
            self.C_global = C_use
            return C_use
        if (C_global is not None and getattr(C_global, "nnz", 0) > 0) and (ra_f != 0.0 or rb_f != 0.0):
            logger.warning(
                "Transient: Rayleigh alpha/beta are ignored when assembled element C is non-empty "
                "(precedence: element C)."
            )
        return C_global

    def modify_dynamic_global_system(
        self, K_global: csr_matrix, M_global: csr_matrix, C_global: Any
    ) -> tuple[csr_matrix, csr_matrix, Any, np.ndarray]:
        """Apply BCs; cache ``self.K_mod``, ``self.M_mod``, ``self.C_mod``."""
        fdyn = self._resolved_penalty_fixed_dofs()
        K_mod, M_mod, C_mod, bc_dyn = ModifyDynamicGlobalSystem(
            fixed_dofs=fdyn,
            prescribed_displacements=self.prescribed_displacement_dict,
            job_results_dir=self.primary_results_dir,
        ).run(K_global, M_global, C_global)
        self.K_mod = K_mod
        self.M_mod = M_mod
        self.C_mod = C_mod
        log_transient_modified_system(
            K_mod,
            M_mod,
            C_mod,
            n_bc_dofs=int(np.asarray(bc_dyn).size),
            job_results_dir=self.primary_results_dir,
        )
        return K_mod, M_mod, C_mod, bc_dyn

    def integrate_transient_system(
        self,
        K_mod: csr_matrix,
        M_mod: csr_matrix,
        C_mod: Any,
        t_grid: np.ndarray,
        F_func: Callable[[float], np.ndarray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Newmark integration; caches ``self.U_time_history``, ``self.V_time_history``, ``self.A_time_history``, ``self.t_grid``."""
        total_dof = int(K_mod.shape[0])
        u0 = np.zeros(total_dof, dtype=np.float64)
        v0 = np.zeros(total_dof, dtype=np.float64)
        U, V, A = IntegrateTransientSystem(
            t_grid=t_grid,
            force_func=F_func,
            job_results_dir=self.primary_results_dir,
        ).run(K_mod, M_mod, C_mod, u0, v0)
        self.U_time_history = U
        self.V_time_history = V
        self.A_time_history = A
        self.t_grid = t_grid
        return U, V, A

    def save_transient_primary_results(self, t_grid: np.ndarray, U: np.ndarray, V: np.ndarray, A: np.ndarray) -> None:
        """Write ``dynamic_results/*.txt`` time histories."""
        results_dir = os.path.join(self.primary_results_dir, "dynamic_results")
        os.makedirs(results_dir, exist_ok=True)
        np.savetxt(os.path.join(results_dir, f"{self.job_name}_time.txt"), t_grid, fmt="%.6f")
        np.savetxt(os.path.join(results_dir, f"{self.job_name}_displacements.txt"), U, fmt="%.6e")
        np.savetxt(os.path.join(results_dir, f"{self.job_name}_velocities.txt"), V, fmt="%.6e")
        np.savetxt(os.path.join(results_dir, f"{self.job_name}_accelerations.txt"), A, fmt="%.6e")

    def run_dynamic_post_processing_if_enabled(self, U: np.ndarray) -> None:
        """Optional formulation-cache secondary/tertiary from displacement snapshots."""
        if not self._dynamic_post_enabled():
            return
        indices = self._dynamic_snapshot_row_indices(U.shape[0])
        use_subdir = len(indices) > 1
        for ti in indices:
            logger.info("Dynamic post-processing snapshot at time index %s", ti)
            sub = os.path.join("dynamic_post", f"t_{ti:06d}") if use_subdir else None
            self._run_secondary_tertiary_from_cache(U[ti], results_subdir=sub)

    def write_transient_primary_artifact_manifest(self) -> None:
        """``logs/primary_artifacts.json`` index for transient outputs."""
        dr = "primary_results/dynamic_results"
        write_primary_artifact_manifest(
            self.results_root,
            family="transient",
            job_name=self.job_name,
            artifacts={
                "time": f"{dr}/{self.job_name}_time.txt",
                "displacements": f"{dr}/{self.job_name}_displacements.txt",
                "velocities": f"{dr}/{self.job_name}_velocities.txt",
                "accelerations": f"{dr}/{self.job_name}_accelerations.txt",
            },
        )

    def run(self):
        """Execute transient workflow: setup, assemble, BCs, time integrate, save."""
        try:
            self.setup_simulation()
            dyn_config = effective_transient_config(self.simulation_settings)
            dt = float(dyn_config.get("time_step", 0.001))
            end_time = float(dyn_config.get("end_time", 1.0))
            n_steps = max(1, int(round(end_time / dt)))
            t_grid = np.linspace(0.0, end_time, n_steps + 1)

            total_dof = self._total_dof()
            monitor = RuntimeMonitorTelemetry(job_results_dir=self.diagnostics_dir)

            with monitor.stage("AssembleDynamicGlobalSystem"):
                K_global, M_global, C_global, _ = self.assemble_dynamic_global_system(total_dof)

            C_global = self.apply_rayleigh_damping_if_needed(K_global, M_global, C_global, dyn_config)

            with monitor.stage("ModifyDynamicGlobalSystem"):
                K_mod, M_mod, C_mod, _bc_dyn = self.modify_dynamic_global_system(
                    K_global, M_global, C_global
                )

            F_ref = self._assemble_reference_force(total_dof)
            F_func = build_transient_force_func(
                F_ref,
                {**dyn_config, "end_time": end_time},
                total_dof=total_dof,
                job_dir=self.settings.get("job_dir"),
                end_time=end_time,
            )

            with monitor.stage("IntegrateTransientSystem"):
                U, V, A = self.integrate_transient_system(K_mod, M_mod, C_mod, t_grid, F_func)

            self.save_transient_primary_results(t_grid, U, V, A)
            self.run_dynamic_post_processing_if_enabled(U)
            self.write_transient_primary_artifact_manifest()

            logger.info("Dynamic simulation completed successfully -> %s", self.results_root)
        except Exception as exc:
            logger.exception("Dynamic simulation failed")
            raise RuntimeError("Dynamic simulation aborted") from exc


class DynamicSimulationRunner(TransientSimulationRunner):
    """Deprecated alias for :class:`TransientSimulationRunner`."""

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "DynamicSimulationRunner is deprecated; use TransientSimulationRunner.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)
