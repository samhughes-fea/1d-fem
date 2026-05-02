# simulation_runner/dynamic/dynamic_simulation.py
"""Dynamic simulation runner: time integration of M u'' + C u' + K u = F(t)."""

import logging
import os
import numpy as np
from scipy.sparse import coo_matrix

from processing.dynamic.assembly import assemble_global_system
from processing.dynamic.boundary_conditions import apply_boundary_conditions
from processing.dynamic.time_integration import newmark_integrate

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
            self.logs_dir = os.path.join(job_results_dir, "logs")
        else:
            self.results_root = os.path.join("post_processing", "results", self.job_name)
            self.primary_results_dir = os.path.join(self.results_root, "primary_results")
            self.logs_dir = os.path.join(self.results_root, "logs")

        self.simulation_settings = self.settings.get("simulation_settings", {})

    def _ensure_sparse_format(self, matrices):
        if matrices is None:
            return None
        return np.array(
            [coo_matrix(m) if not isinstance(m, coo_matrix) else m for m in matrices],
            dtype=object,
        )

    def setup_simulation(self):
        """Create output directory structure."""
        logger.info("Setting up dynamic simulation for job: %s", self.job_name)
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        logger.info("Results will be saved to: %s", self.results_root)

    def run(self):
        """Execute dynamic workflow: setup, assemble, BCs, time integrate, save."""
        try:
            self.setup_simulation()
            dyn_config = self.simulation_settings.get("dynamic", {})
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
            logger.info("Dynamic simulation completed successfully -> %s", self.results_root)
        except Exception as exc:
            logger.exception("Dynamic simulation failed")
            raise RuntimeError("Dynamic simulation aborted") from exc
