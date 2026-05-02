from __future__ import annotations

import copy
import logging
import numpy as np
import os
import datetime
from scipy.sparse import coo_matrix, linalg

from processing.modal.assembly import assemble_global_matrices
from processing.modal.boundary_conditions import apply_boundary_conditions
from processing.modal.buckling import (
    apply_buckling_boundary_conditions,
    assemble_global_geometric_stiffness,
    solve_linear_buckling_eigenpairs,
)
from simulation_runner.modal.modal_diagnostic import log_modal_diagnostics

logger = logging.getLogger(__name__)

_LINEAR_TO_NONLINEAR_PRESTRESS_TWINS = {
    "LinearTimoshenkoBeamElement3D": "NonlinearTimoshenkoBeamElement3D",
    "LinearEulerBernoulliBeamElement3D": "NonlinearEulerBernoulliBeamElement3D",
}


def _nonlinear_twin_element_dictionary(element_dictionary: dict) -> dict:
    """Return a copy of *element_dictionary* with beam types swapped to nonlinear twins."""
    ed = dict(element_dictionary)
    # Must use dtype=object so longer nonlinear type names are not truncated (Unicode arrays).
    types = np.asarray(ed["types"], dtype=object).copy()
    for i, t in enumerate(types):
        ts = str(t)
        if ts not in _LINEAR_TO_NONLINEAR_PRESTRESS_TWINS:
            raise ValueError(
                "modal.buckling_nonlinear_prestress_twins: unsupported element type "
                f"{ts!r}; registered twins for "
                f"{sorted(_LINEAR_TO_NONLINEAR_PRESTRESS_TWINS.keys())}"
            )
        types[i] = _LINEAR_TO_NONLINEAR_PRESTRESS_TWINS[ts]
    ed["types"] = types
    return ed


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

    def _total_dof(self) -> int:
        num_nodes = len(self.mesh_dictionary["node_ids"])
        from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

        ed = self.settings.get("element_dictionary")
        dpn = 7 if ed is not None and mesh_uses_warping_dof(ed) else 6
        return int(num_nodes * dpn)

    def _assemble_global_matrices(self, job_results_dir):
        """Assemble global stiffness and mass matrices."""
        logger.info("Assembling global stiffness and mass matrices...")
        total_dof = self._total_dof()

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

    def solve_modal_vibration(self, K_mod, M_mod, num_modes, job_results_dir):
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

    def _save_buckling_results(self, load_factors: np.ndarray, mode_shapes: np.ndarray) -> None:
        results_dir = os.path.join(self.primary_results_dir, "modal_results")
        os.makedirs(results_dir, exist_ok=True)
        np.savetxt(
            os.path.join(results_dir, f"{self.job_name}_buckling_load_factors.txt"),
            load_factors,
            fmt="%.8e",
        )
        np.savetxt(
            os.path.join(results_dir, f"{self.job_name}_buckling_mode_shapes.txt"),
            mode_shapes,
            fmt="%.8e",
        )
        logger.info("Saved linear buckling results (load factors and modes).")

    def _instantiate_nonlinear_prestress_twins(self, twin_root_dir: str):
        """
        Nonlinear beam elements plus cached ``element_objects`` / ``force_objects`` for prestress.

        Used when ``modal.buckling_nonlinear_prestress_twins`` is true so ``element.txt`` may list
        linear types while the nonlinear static prestress run uses registered nonlinear twins.
        """
        from pre_processing.element_library.element_factory import ElementFactory

        ed = self.settings["element_dictionary"]
        nl_ed = _nonlinear_twin_element_dictionary(ed)
        twin_jrd = os.path.join(twin_root_dir, "twin_elements")
        for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
            os.makedirs(os.path.join(twin_jrd, sub), exist_ok=True)
        factory = ElementFactory(job_results_dir=twin_jrd)
        eids = np.asarray(ed["ids"], dtype=np.int64)
        nl_elements = factory.create_elements_batch(
            element_ids=eids,
            element_dictionary=nl_ed,
            grid_dictionary=self.settings["grid_dictionary"],
            material_dictionary=self.settings["material_dictionary"],
            section_dictionary=self.settings["section_dictionary"],
            point_load_array=self.settings.get("point_load_array", np.empty((0, 9))),
            distributed_load_array=self.settings.get("distributed_load_array", np.empty((0, 9))),
            simulation_settings=self.simulation_settings,
        )
        nl_eos = [e.element_stiffness_matrix() for e in nl_elements]
        nl_fos = [e.element_force_vector() for e in nl_elements]
        return list(nl_elements), nl_eos, nl_fos

    def _run_linear_buckling(
        self,
        K_global,
        job_results_dir: str,
        num_modes: int,
        modal_cfg: dict,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Prestress via linear static solve, assemble :math:`\\mathbf{K}_g`, solve buckling eigenpairs."""
        prestress = str(modal_cfg.get("buckling_prestress", "linear_static")).lower()
        if prestress == "none":
            raise ValueError(
                "modal.buckling_prestress='none' cannot define reference internal forces for linear buckling "
                "(geometric stiffness K_sigma would be zero). Use buckling_prestress=linear_static or nonlinear_static."
            )
        if prestress not in ("linear_static", "nonlinear_static"):
            raise ValueError(
                "modal.buckling_prestress must be 'linear_static' or 'nonlinear_static' "
                "(reference prestress for K_sigma)."
            )
        eo = self.settings.get("element_objects")
        fo = self.settings.get("force_objects")
        if eo is None or fo is None:
            raise ValueError(
                "Buckling requires element_objects and force_objects in modal settings (run_job wiring)."
            )
        element_objects_list = list(np.asarray(eo, dtype=object).ravel())
        force_objects_list = list(np.asarray(fo, dtype=object).ravel())
        load_scale = float(modal_cfg.get("buckling_load_factor", 1.0))

        if prestress == "linear_static":
            from simulation_runner.static.linear_static_simulation import LinearStaticSimulationRunner

            static_runner = LinearStaticSimulationRunner(
                elements=list(self.elements),
                grid_dictionary=self.settings["grid_dictionary"],
                element_dictionary=self.settings["element_dictionary"],
                material_dictionary=self.settings["material_dictionary"],
                section_dictionary=self.settings["section_dictionary"],
                point_load_array=self.settings.get("point_load_array", np.empty((0, 9))),
                distributed_load_array=self.settings.get("distributed_load_array", np.empty((0, 9))),
                element_objects=np.asarray(element_objects_list, dtype=object),
                force_objects=np.asarray(force_objects_list, dtype=object),
                job_name=f"{self.job_name}_prestress",
                job_results_dir=os.path.join(job_results_dir, "linear_prestress"),
                simulation_settings=self.simulation_settings,
            )
            pd = self.settings.get("prescribed_displacement_dict")
            if pd is not None:
                static_runner.prescribed_displacements = pd

            static_runner.solve_linear_system_only(force_scale=load_scale)
            U_global = static_runner.U_global
        else:
            use_twins = bool(modal_cfg.get("buckling_nonlinear_prestress_twins", False))
            nl_elements_run = list(self.elements)
            nl_eo = element_objects_list
            nl_fo = force_objects_list
            nl_element_dictionary = self.settings["element_dictionary"]
            if use_twins:
                nl_elements_run, nl_eo, nl_fo = self._instantiate_nonlinear_prestress_twins(
                    os.path.join(job_results_dir, "nonlinear_prestress")
                )
                nl_element_dictionary = _nonlinear_twin_element_dictionary(
                    self.settings["element_dictionary"]
                )
            else:
                for i, elem in enumerate(self.elements):
                    if not callable(getattr(elem, "tangent_stiffness_matrix", None)):
                        raise ValueError(
                            "modal.buckling_prestress='nonlinear_static' requires nonlinear beam elements "
                            "(tangent_stiffness_matrix), or set modal.buckling_nonlinear_prestress_twins=true "
                            "to build nonlinear twins from linear types. "
                            f"Element {i} ({type(elem).__name__}) is not supported."
                        )

            from simulation_runner.static.nonlinear_static_simulation import NonlinearStaticSimulationRunner

            nl_settings = copy.deepcopy(self.simulation_settings)
            nl_cfg = nl_settings.setdefault("nonlinear", {})
            nl_cfg["num_increments"] = 1
            nl_cfg["load_factors"] = [load_scale]

            nl_runner = NonlinearStaticSimulationRunner(
                elements=nl_elements_run,
                grid_dictionary=self.settings["grid_dictionary"],
                element_dictionary=nl_element_dictionary,
                material_dictionary=self.settings["material_dictionary"],
                section_dictionary=self.settings["section_dictionary"],
                point_load_array=self.settings.get("point_load_array", np.empty((0, 9))),
                distributed_load_array=self.settings.get("distributed_load_array", np.empty((0, 9))),
                element_objects=np.asarray(nl_eo, dtype=object),
                force_objects=np.asarray(nl_fo, dtype=object),
                job_name=f"{self.job_name}_prestress_nl",
                job_results_dir=os.path.join(job_results_dir, "nonlinear_prestress"),
                simulation_settings=nl_settings,
            )
            pd = self.settings.get("prescribed_displacement_dict")
            if pd is not None:
                nl_runner.prescribed_displacements = pd

            nl_runner.run()
            if not getattr(nl_runner, "newton_converged", False):
                raise RuntimeError(
                    "Nonlinear prestress for buckling did not converge (Newton). "
                    "Adjust [Newton] / [Nonlinear] settings, loads, or boundary conditions."
                )
            U_global = nl_runner.U_global
        total_dof = self._total_dof()
        Kg_global = assemble_global_geometric_stiffness(
            list(self.elements), U_global, total_dof
        )
        prescribed = self.settings.get("prescribed_displacement_dict")
        K_mod, Kg_mod, bc_dofs = apply_buckling_boundary_conditions(
            K_global, Kg_global, prescribed_displacements=prescribed, fixed_dofs=None
        )
        log_modal_diagnostics(K_mod, Kg_mod, job_results_dir)
        lambdas, modes = solve_linear_buckling_eigenpairs(
            K_mod, Kg_mod, num_modes, constrained_dofs=bc_dofs
        )
        return lambdas, modes

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
        """Execute modal vibration or linear buckling (``modal.analysis`` in simulation settings)."""
        try:
            self.setup_simulation()
            modal_cfg = self.simulation_settings.get("modal", {})
            num_modes = int(modal_cfg.get("num_modes", 10))
            analysis = str(modal_cfg.get("analysis", "vibration")).lower()
            job_results_dir = self.primary_results_dir

            K_global, M_global = self._assemble_global_matrices(job_results_dir)

            if analysis == "buckling":
                lambdas, mode_shapes = self._run_linear_buckling(
                    K_global, job_results_dir, num_modes, modal_cfg
                )
                self._save_buckling_results(lambdas, mode_shapes)
                self.secondary_results["global"]["buckling_load_factors"] = lambdas
                logger.info("Modal buckling simulation completed -> %s", self.results_root)
                return

            K_mod, M_mod, _ = self._modify_global_matrices(K_global, M_global, job_results_dir)
            frequencies, mode_shapes = self.solve_modal_vibration(
                K_mod, M_mod, num_modes, job_results_dir
            )

            self._save_primary_results(frequencies, mode_shapes)
            self._compute_secondary_results(frequencies, mode_shapes)
            self._save_secondary_results()

            logger.info("Modal simulation completed successfully -> %s", self.results_root)
        except Exception as exc:
            logger.exception("Modal simulation failed")
            raise RuntimeError("Modal simulation aborted") from exc