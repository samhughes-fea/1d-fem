# simulation_runner/transient/dynamic_simulation.py
"""Transient dynamics (§3): time integration M u'' + C u' + K u = F(t)."""

import logging
import os
import numpy as np
from scipy.sparse import coo_matrix

from processing.dynamic.assembly import assemble_global_system
from processing.dynamic.boundary_conditions import apply_boundary_conditions
from processing.dynamic.time_integration import newmark_integrate
from pre_processing.parsing.simulation_settings_resolution import effective_transient_config

logger = logging.getLogger(__name__)


class DynamicSimulationRunner:
    """
    Handles dynamic finite element analysis (time integration).
    """

    def __init__(self, settings: dict, job_name: str):
        self.settings = settings
        self.job_name = job_name
        self.elements = self.settings.get("elements", np.array([]))
        self.mesh_dictionary = self.settings.get("mesh_dictionary", {})
        if self.elements.size == 0 or not self.mesh_dictionary:
            raise ValueError("Missing elements or mesh data in settings")

        self.grid_dictionary = self.settings.get("grid_dictionary")
        self.element_dictionary = self.settings.get("element_dictionary")
        self.material_dictionary = self.settings.get("material_dictionary")
        self.section_dictionary = self.settings.get("section_dictionary")
        self.point_load_array = self.settings.get("point_load_array", np.empty((0, 9)))
        self.distributed_load_array = self.settings.get("distributed_load_array", np.empty((0, 9)))
        self.element_objects = self.settings.get("element_objects")
        self.force_objects = self.settings.get("force_objects")

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
            self.logs_dir = os.path.join(job_results_dir, "logs")
        else:
            self.results_root = os.path.join("post_processing", "results", self.job_name)
            self.primary_results_dir = os.path.join(self.results_root, "primary_results")
            self.secondary_results_dir = os.path.join(self.results_root, "secondary_results")
            self.tertiary_results_dir = os.path.join(self.results_root, "tertiary_results")
            self.logs_dir = os.path.join(self.results_root, "logs")

        self.simulation_settings = self.settings.get("simulation_settings", {})

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

    def _run_secondary_tertiary_from_cache(self, U_global: np.ndarray) -> None:
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
        run_secondary_tertiary_from_formulation_cache(
            elements=list(self.elements),
            grid_dictionary=self.grid_dictionary,
            element_dictionary=self.element_dictionary,
            material_dictionary=self.material_dictionary,
            section_dictionary=self.section_dictionary,
            U_global=U,
            formulation_cache=cache,
            results_root=self.results_root,
            secondary_results_dir=self.secondary_results_dir,
            tertiary_results_dir=self.tertiary_results_dir,
        )

    def setup_simulation(self):
        """Create output directory structure."""
        logger.info("Setting up dynamic simulation for job: %s", self.job_name)
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.secondary_results_dir, exist_ok=True)
        os.makedirs(self.tertiary_results_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        logger.info("Results will be saved to: %s", self.results_root)

    def run(self):
        """Execute transient workflow: setup, assemble, BCs, time integrate, save."""
        try:
            self.setup_simulation()
            dyn_config = effective_transient_config(self.simulation_settings)
            dt = dyn_config.get("time_step", 0.001)
            end_time = dyn_config.get("end_time", 1.0)
            n_steps = max(1, int(round(end_time / dt)))
            t_grid = np.linspace(0.0, end_time, n_steps + 1)

            num_nodes = len(self.mesh_dictionary["node_ids"])
            from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

            ed = self.settings.get("element_dictionary")
            dpn = 7 if ed is not None and mesh_uses_warping_dof(ed) else 6
            total_dof = num_nodes * dpn

            K_global, M_global, C_global, _ = assemble_global_system(
                elements=list(self.elements),
                element_stiffness_matrices=self.element_stiffness_matrices,
                element_mass_matrices=self.element_mass_matrices,
                element_damping_matrices=self.element_damping_matrices,
                total_dof=total_dof,
                job_results_dir=self.primary_results_dir,
            )

            K_mod, M_mod, C_mod, _ = apply_boundary_conditions(
                K_global, M_global, C_global, fixed_dofs=list(range(6))
            )

            u0 = np.zeros(total_dof, dtype=np.float64)
            v0 = np.zeros(total_dof, dtype=np.float64)

            def F_func(t):
                return np.zeros(total_dof, dtype=np.float64)

            U, V, A = newmark_integrate(
                K_mod, M_mod, C_mod, u0, v0, t_grid, F_func
            )

            results_dir = os.path.join(self.primary_results_dir, "dynamic_results")
            os.makedirs(results_dir, exist_ok=True)
            np.savetxt(os.path.join(results_dir, f"{self.job_name}_time.txt"), t_grid, fmt="%.6f")
            np.savetxt(os.path.join(results_dir, f"{self.job_name}_displacements.txt"), U, fmt="%.6e")
            np.savetxt(os.path.join(results_dir, f"{self.job_name}_velocities.txt"), V, fmt="%.6e")
            np.savetxt(os.path.join(results_dir, f"{self.job_name}_accelerations.txt"), A, fmt="%.6e")

            if self._dynamic_post_enabled():
                ti = self._resolve_time_index(U.shape[0])
                logger.info("Dynamic post-processing snapshot at time index %s", ti)
                self._run_secondary_tertiary_from_cache(U[ti])

            logger.info("Dynamic simulation completed successfully -> %s", self.results_root)
        except Exception as exc:
            logger.exception("Dynamic simulation failed")
            raise RuntimeError("Dynamic simulation aborted") from exc
