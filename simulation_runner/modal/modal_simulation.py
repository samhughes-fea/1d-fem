import logging
import numpy as np
import os
import datetime
from scipy.sparse import coo_matrix, linalg

from processing.modal.assembly import assemble_global_matrices
from processing.modal.boundary_conditions import apply_boundary_conditions
from simulation_runner.modal.modal_diagnostic import log_modal_diagnostics

logger = logging.getLogger(__name__)


class ModalSimulationRunner:
    """
    Handles modal finite element analysis (natural frequencies and mode shapes).
    """

    def __init__(self, settings, job_name):
        self.settings = settings
        self.job_name = job_name
        self.start_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        self.primary_results = {"global": {}, "element": {"data": []}}
        self.secondary_results = {"global": {}, "element": {"data": []}}

        self.elements = self.settings.get("elements", np.array([]))
        self.mesh_dictionary = self.settings.get("mesh_dictionary", {})

        if self.elements.size == 0 or not self.mesh_dictionary:
            logger.error("Missing elements or mesh data in settings!")
            raise ValueError("Missing elements or mesh data in settings!")

        self.solver_name = self.settings.get("solver_name", "eigen")
        self.element_stiffness_matrices = self._ensure_sparse_format(
            self.settings.get("element_stiffness_matrices", None)
        )
        self.element_mass_matrices = self._ensure_sparse_format(
            self.settings.get("element_mass_matrices", None)
        )

        job_results_dir = self.settings.get("job_results_dir")
        if job_results_dir:
            self.results_root = job_results_dir
            self.primary_results_dir = os.path.join(job_results_dir, "primary_results")
            self.diagnostics_dir = os.path.join(job_results_dir, "diagnostics")
            self.logs_dir = os.path.join(job_results_dir, "logs")
        else:
            self.results_root = os.path.join(
                "post_processing", "results", f"{self.job_name}_{self.start_time}"
            )
            self.primary_results_dir = os.path.join(self.results_root, "primary_results")
            self.diagnostics_dir = os.path.join(self.results_root, "diagnostics")
            self.logs_dir = os.path.join(self.results_root, "logs")

        self.simulation_settings = self.settings.get("simulation_settings", {})

    def _ensure_sparse_format(self, matrices):
        """Converts matrices to sparse COO format if needed."""
        if matrices is None:
            return None
        return np.array([
            coo_matrix(matrix) if not isinstance(matrix, coo_matrix) else matrix
            for matrix in matrices
        ], dtype=object)

    def setup_simulation(self):
        """Create output directory structure under results root."""
        logger.info("Setting up modal simulation for job: %s", self.job_name)
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.diagnostics_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        logger.info("Results will be saved to: %s", self.results_root)

    # -------------------------------------------------------------------------
    # 1) ASSEMBLE GLOBAL MATRICES
    # -------------------------------------------------------------------------

    def _assemble_global_matrices(self, job_results_dir):
        """Assemble global stiffness and mass matrices."""
        logger.info("Assembling global stiffness and mass matrices...")
        num_nodes = len(self.mesh_dictionary["node_ids"])
        total_dof = num_nodes * 6

        K_global, M_global, _ = assemble_global_matrices(
            elements=list(self.elements),
            element_stiffness_matrices=self.element_stiffness_matrices,
            element_mass_matrices=self.element_mass_matrices,
            total_dof=total_dof,
            job_results_dir=job_results_dir,
        )
        if K_global is None or M_global is None:
            raise ValueError("Global matrices could not be assembled")
        log_modal_diagnostics(K_global, M_global, job_results_dir)
        logger.info("Global stiffness and mass matrices assembled.")
        return K_global, M_global

    # -------------------------------------------------------------------------
    # 2) APPLY BOUNDARY CONDITIONS
    # -------------------------------------------------------------------------

    def _modify_global_matrices(self, K_global, M_global, job_results_dir):
        """Apply boundary conditions to the modal system."""
        logger.info("Applying boundary conditions to global matrices...")
        K_mod, M_mod, bc_dofs = apply_boundary_conditions(K_global, M_global)
        log_modal_diagnostics(K_mod, M_mod, job_results_dir)
        logger.info("Boundary conditions applied.")
        return K_mod, M_mod, bc_dofs

    # -------------------------------------------------------------------------
    # 3) SOLVE MODAL SYSTEM
    # -------------------------------------------------------------------------

    def solve_modal(self, K_mod, M_mod, num_modes, job_results_dir):
        """
        Solves the modal system for natural frequencies and mode shapes.

        Parameters:
            K_mod (csr_matrix): Modified stiffness matrix.
            M_mod (csr_matrix): Modified mass matrix.
            num_modes (int): Number of modes to compute.
            job_results_dir (str): Directory for logging results.

        Returns:
            frequencies (np.ndarray): Natural frequencies (Hz).
            mode_shapes (np.ndarray): Mode shape vectors.
        """
        logger.info(f"🔹 Solving for {num_modes} natural frequencies and mode shapes...")

        try:
            eigenvalues, eigenvectors = linalg.eigsh(K_mod, k=num_modes, M=M_mod, which="SM")

            frequencies = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
            mode_shapes = eigenvectors

            logger.info(f"✅ Computed {num_modes} natural frequencies.")

            return frequencies, mode_shapes
        except Exception as e:
            logger.error(f"❌ Modal solver failure: {e}")
            raise

    # -------------------------------------------------------------------------
    # 4) SAVE PRIMARY RESULTS
    # -------------------------------------------------------------------------

    def _save_primary_results(self, frequencies, mode_shapes):
        """Save natural frequencies and mode shapes to files."""
        results_dir = os.path.join(self.primary_results_dir, "modal_results")
        os.makedirs(results_dir, exist_ok=True)
        np.savetxt(os.path.join(results_dir, f"{self.job_name}_frequencies.txt"), frequencies, fmt="%.6f")
        np.savetxt(os.path.join(results_dir, f"{self.job_name}_mode_shapes.txt"), mode_shapes, fmt="%.6f")
        logger.info("Saved modal results.")

    # -------------------------------------------------------------------------
    # 5) COMPUTE SECONDARY RESULTS (PLACEHOLDERS)
    # -------------------------------------------------------------------------

    def _compute_secondary_results(self, frequencies, mode_shapes):
        """Compute secondary modal results (placeholder)."""
        self.secondary_results["global"]["modal_participation"] = np.array([0.0])
        logger.info("Computed secondary modal results.")

    def _save_secondary_results(self):
        """Save secondary modal results."""
        results_dir = os.path.join(self.primary_results_dir, "modal_results")
        os.makedirs(results_dir, exist_ok=True)
        np.savetxt(
            os.path.join(results_dir, f"{self.job_name}_modal_participation.txt"),
            self.secondary_results["global"]["modal_participation"],
            fmt="%.6f",
        )
        logger.info("Saved secondary modal results.")

    # -------------------------------------------------------------------------
    # RUN
    # -------------------------------------------------------------------------

    def run(self):
        """Execute the full modal workflow: setup, assemble, BCs, solve, save."""
        try:
            self.setup_simulation()
            num_modes = self.simulation_settings.get("modal", {}).get("num_modes", 10)
            job_results_dir = self.primary_results_dir

            K_global, M_global = self._assemble_global_matrices(job_results_dir)
            K_mod, M_mod, _ = self._modify_global_matrices(K_global, M_global, job_results_dir)
            frequencies, mode_shapes = self.solve_modal(K_mod, M_mod, num_modes, job_results_dir)

            self._save_primary_results(frequencies, mode_shapes)
            self._compute_secondary_results(frequencies, mode_shapes)
            self._save_secondary_results()

            logger.info("Modal simulation completed successfully -> %s", self.results_root)
        except Exception as exc:
            logger.exception("Modal simulation failed")
            raise RuntimeError("Modal simulation aborted") from exc