from __future__ import annotations

import copy
from contextlib import nullcontext
import logging
import numpy as np
import os
import datetime
from scipy.sparse import coo_matrix
from scipy.sparse import spmatrix

from processing.static.diagnostics.runtime_monitor_telemetry import RuntimeMonitorTelemetry
from processing.spectral.operations import (
    AssembleBucklingGeometricStiffness,
    AssembleSpectralGlobalSystem,
    ModifyBucklingGlobalMatrices,
    ModifySpectralGlobalSystem,
    PrepareSpectralLocalMatrices,
    SolveGeneralizedEigenproblem,
    SolveLinearBucklingEigenpairs,
)
from pre_processing.parsing.simulation_settings_resolution import (
    effective_buckling_config,
    effective_eigen_config,
)
from processing.boundary_supports import resolve_penalty_fixed_dofs
from processing.common.primary_artifact_manifest import write_primary_artifact_manifest
from processing.dynamic.assembly import assemble_global_force_vector
from processing.eigen.metrics.directional_effective_mass import modal_effective_mass_fraction_z
from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

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


class VibrationBucklingBackend:
    """
    Shared implementation for §2 eigen vibration and §5 linear buckling (undamped
    :math:`K` / :math:`M` and buckling pencil). Subclass with a concrete ``run()`` in
    :class:`~simulation_runner.eigen.eigen_simulation.EigenSimulationRunner` or
    :class:`~simulation_runner.buckling.buckling_simulation.BucklingSimulationRunner`.
    """

    def __init__(self, settings, job_name):
        self.settings = settings
        self.job_name = job_name
        self.start_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        self.primary_results = {"global": {}, "element": {"data": []}}
        self.secondary_results = {"global": {}, "element": {"data": []}}

        raw_el = self.settings.get("elements", np.array([]))
        self.elements = (
            raw_el if isinstance(raw_el, np.ndarray) else np.asarray(raw_el, dtype=object)
        )
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
            self.secondary_results_dir = os.path.join(job_results_dir, "secondary_results")
            self.tertiary_results_dir = os.path.join(job_results_dir, "tertiary_results")
            self.diagnostics_dir = os.path.join(job_results_dir, "diagnostics")
            self.logs_dir = os.path.join(job_results_dir, "logs")
        else:
            self.results_root = os.path.join(
                "post_processing", "results", f"{self.job_name}_{self.start_time}"
            )
            self.primary_results_dir = os.path.join(self.results_root, "primary_results")
            self.secondary_results_dir = os.path.join(self.results_root, "secondary_results")
            self.tertiary_results_dir = os.path.join(self.results_root, "tertiary_results")
            self.diagnostics_dir = os.path.join(self.results_root, "diagnostics")
            self.logs_dir = os.path.join(self.results_root, "logs")

        self.simulation_settings = self.settings.get("simulation_settings", {})

    def _ensure_sparse_format(self, matrices):
        """Converts matrices to sparse COO format if needed."""
        if matrices is None:
            return None
        out = []
        for matrix in matrices:
            if isinstance(matrix, spmatrix):
                out.append(matrix.tocoo() if not isinstance(matrix, coo_matrix) else matrix)
            else:
                dense = np.asarray(matrix, dtype=np.float64)
                out.append(coo_matrix(dense))
        return np.array(out, dtype=object)

    def setup_simulation(self):
        """Create output directory structure under results root."""
        logger.info("Setting up modal simulation for job: %s", self.job_name)
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.secondary_results_dir, exist_ok=True)
        os.makedirs(self.tertiary_results_dir, exist_ok=True)
        os.makedirs(self.diagnostics_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        logger.info("Results will be saved to: %s", self.results_root)

    # -------------------------------------------------------------------------
    # 1) ASSEMBLE GLOBAL MATRICES
    # -------------------------------------------------------------------------

    def _dof_per_node(self) -> int:
        ed = self.settings.get("element_dictionary")
        return 7 if ed is not None and mesh_uses_warping_dof(ed) else 6

    def _total_dof(self) -> int:
        num_nodes = len(self.mesh_dictionary["node_ids"])
        return int(num_nodes * self._dof_per_node())

    def _resolved_penalty_fixed_dofs(self) -> np.ndarray | None:
        """Penalty BC ``fixed_dofs`` from prescribed zeros + optional ``[Eigen] fixed_node_id``."""
        return resolve_penalty_fixed_dofs(
            total_dof=self._total_dof(),
            dof_per_node=self._dof_per_node(),
            prescribed_displacement_dict=self.settings.get("prescribed_displacement_dict"),
            section_settings=self.simulation_settings.get("eigen") or {},
            grid_node_ids=self.mesh_dictionary.get("node_ids"),
        )

    def _assemble_global_matrices(self, job_results_dir, ke_list=None, me_list=None):
        """Assemble global K and M via :class:`AssembleSpectralGlobalSystem` (optional pre-prepared lists)."""
        logger.info("Assembling global stiffness and mass matrices...")
        total_dof = self._total_dof()
        if ke_list is None or me_list is None:
            ke_list, me_list = PrepareSpectralLocalMatrices(
                element_stiffness_matrices=self.element_stiffness_matrices,
                element_mass_matrices=self.element_mass_matrices,
                job_results_dir=job_results_dir,
            ).run()
        K_global, M_global, _ = AssembleSpectralGlobalSystem(
            elements=list(self.elements),
            element_stiffness_matrices=ke_list,
            element_mass_matrices=me_list,
            total_dof=total_dof,
            job_results_dir=job_results_dir,
        ).run()
        logger.info("Global stiffness and mass matrices assembled.")
        return K_global, M_global

    # -------------------------------------------------------------------------
    # 2) APPLY BOUNDARY CONDITIONS
    # -------------------------------------------------------------------------

    def _modify_global_matrices(self, K_global, M_global, job_results_dir):
        """Apply boundary conditions via :class:`ModifySpectralGlobalSystem`."""
        logger.info("Applying boundary conditions to global matrices...")
        prescribed = self.settings.get("prescribed_displacement_dict")
        fixed_dofs = self._resolved_penalty_fixed_dofs()
        K_mod, M_mod, bc_dofs = ModifySpectralGlobalSystem(
            job_results_dir=job_results_dir,
            fixed_dofs=fixed_dofs,
            prescribed_displacements=prescribed,
        ).run(K_global, M_global)
        logger.info("Boundary conditions applied.")
        return K_mod, M_mod, bc_dofs

    # -------------------------------------------------------------------------
    # 3) SOLVE MODAL SYSTEM
    # -------------------------------------------------------------------------

    def solve_modal_vibration(self, K_mod, M_mod, num_modes, job_results_dir, *, dense_threshold=512):
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
        logger.info("Solving for %s natural frequencies and mode shapes...", num_modes)
        try:
            _ev, mode_shapes, frequencies = SolveGeneralizedEigenproblem(
                num_modes=num_modes,
                context="eigen vibration",
                dense_threshold=int(dense_threshold),
                job_results_dir=job_results_dir,
            ).run(K_mod, M_mod)
            logger.info("Computed %s natural frequencies.", num_modes)
            return frequencies, mode_shapes
        except Exception as e:
            logger.error("Modal solver failure: %s", e)
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
        mr = "primary_results/modal_results"
        write_primary_artifact_manifest(
            self.results_root,
            family="eigen_vibration",
            job_name=self.job_name,
            artifacts={
                "frequencies_hz": f"{mr}/{self.job_name}_frequencies.txt",
                "mode_shapes": f"{mr}/{self.job_name}_mode_shapes.txt",
            },
        )

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
        mr = "primary_results/modal_results"
        write_primary_artifact_manifest(
            self.results_root,
            family="buckling",
            job_name=self.job_name,
            artifacts={
                "buckling_load_factors": f"{mr}/{self.job_name}_buckling_load_factors.txt",
                "buckling_mode_shapes": f"{mr}/{self.job_name}_buckling_mode_shapes.txt",
            },
        )

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
        buck_cfg: dict,
        monitor: RuntimeMonitorTelemetry | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Prestress via linear static solve, assemble :math:`\\mathbf{K}_g`, solve buckling eigenpairs."""
        prestress = str(buck_cfg.get("buckling_prestress", "linear_static")).lower()
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

        def _stage(name: str):
            return monitor.stage(name) if monitor is not None else nullcontext()

        eo = self.settings.get("element_objects")
        fo = self.settings.get("force_objects")
        if eo is None or fo is None:
            raise ValueError(
                "Buckling requires element_objects and force_objects in modal settings (run_job wiring)."
            )
        element_objects_list = list(np.asarray(eo, dtype=object).ravel())
        force_objects_list = list(np.asarray(fo, dtype=object).ravel())
        load_scale = float(buck_cfg.get("buckling_load_factor", 1.0))

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

            with _stage("LinearStaticPrestress"):
                static_runner.solve_linear_system_only(force_scale=load_scale)
                U_global = static_runner.U_global
        else:
            use_twins = bool(buck_cfg.get("buckling_nonlinear_prestress_twins", False))
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

            with _stage("NonlinearStaticPrestress"):
                nl_runner.run()
                if not getattr(nl_runner, "newton_converged", False):
                    raise RuntimeError(
                        "Nonlinear prestress for buckling did not converge (Newton). "
                        "Adjust [Newton] / [Nonlinear] settings, loads, or boundary conditions."
                    )
                U_global = nl_runner.U_global
        self._buckling_prestress_u = np.asarray(U_global, dtype=np.float64).copy()
        total_dof = self._total_dof()

        with _stage("AssembleBucklingGeometricStiffness"):
            Kg_global = AssembleBucklingGeometricStiffness(
                elements=list(self.elements),
                U_global=U_global,
                total_dof=total_dof,
                job_results_dir=job_results_dir,
            ).run()
        prescribed = self.settings.get("prescribed_displacement_dict")
        buck_fixed = self._resolved_penalty_fixed_dofs()
        with _stage("ModifyBucklingGlobalMatrices"):
            K_mod, Kg_mod, bc_dofs = ModifyBucklingGlobalMatrices(
                job_results_dir=job_results_dir,
                prescribed_displacements=prescribed,
                fixed_dofs=buck_fixed,
            ).run(K_global, Kg_global)
        with _stage("SolveLinearBucklingEigenpairs"):
            lambdas, modes = SolveLinearBucklingEigenpairs(
                num_modes=num_modes,
                job_results_dir=job_results_dir,
            ).run(K_mod, Kg_mod, bc_dofs)
        return lambdas, modes

    # -------------------------------------------------------------------------
    # 5) OPTIONAL: formulation-cache secondary / tertiary (static-style pipeline)
    # -------------------------------------------------------------------------

    def _modal_post_processing_cfg(self) -> dict:
        return self.simulation_settings.get("post_processing") or {}

    def _modal_post_processing_enabled(self) -> bool:
        return bool(self._modal_post_processing_cfg().get("run_secondary_tertiary_modal", False))

    def _snapshot_u_vibration(self, mode_shapes: np.ndarray) -> np.ndarray:
        cfg = self._modal_post_processing_cfg()
        k = int(cfg.get("modal_mode_index", 0))
        amp = float(cfg.get("modal_amplitude", 1.0))
        ms = np.asarray(mode_shapes, dtype=np.float64)
        if ms.ndim != 2:
            raise ValueError("mode_shapes must be 2-D (n_dof × n_modes)")
        if k < 0 or k >= ms.shape[1]:
            raise IndexError(f"modal_mode_index={k} out of range for {ms.shape[1]} modes")
        return ms[:, k] * amp

    def _snapshot_u_buckling(self, buckling_modes: np.ndarray) -> np.ndarray:
        cfg = self._modal_post_processing_cfg()
        policy = str(cfg.get("buckling_displacement", "mode")).strip().lower()
        amp = float(cfg.get("modal_amplitude", 1.0))
        if policy == "prestress":
            u0 = getattr(self, "_buckling_prestress_u", None)
            if u0 is None:
                raise RuntimeError("buckling prestress displacement not stored")
            return np.asarray(u0, dtype=np.float64).ravel() * amp
        if policy != "mode":
            raise ValueError("post_processing.buckling_displacement must be 'mode' or 'prestress'")
        k = int(cfg.get("buckling_mode_index", cfg.get("modal_mode_index", 0)))
        bm = np.asarray(buckling_modes, dtype=np.float64)
        if bm.ndim != 2:
            raise ValueError("buckling mode matrix must be 2-D")
        if k < 0 or k >= bm.shape[1]:
            raise IndexError(f"buckling_mode_index={k} out of range for {bm.shape[1]} modes")
        return bm[:, k] * amp

    def _run_secondary_tertiary_from_formulation_cache(self, U_global: np.ndarray) -> None:
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
                "post_processing.run_secondary_tertiary_modal requires element_objects and force_objects "
                "(run_job wiring)."
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
        run_secondary_tertiary_from_formulation_cache(
            elements=list(self.elements),
            grid_dictionary=self.settings["grid_dictionary"],
            element_dictionary=self.settings["element_dictionary"],
            material_dictionary=self.settings["material_dictionary"],
            section_dictionary=self.settings["section_dictionary"],
            U_global=U,
            formulation_cache=cache,
            results_root=self.results_root,
            secondary_results_dir=self.secondary_results_dir,
            tertiary_results_dir=self.tertiary_results_dir,
        )

    # -------------------------------------------------------------------------
    # 6) COMPUTE SECONDARY RESULTS (PLACEHOLDERS)
    # -------------------------------------------------------------------------

    def _assemble_modal_reference_nodal_load(self) -> np.ndarray | None:
        """Global reference load from element ``F_e`` (same scatter as transient); no static imports."""
        fo = self.settings.get("force_objects")
        if fo is None:
            return None
        el_list = list(np.asarray(self.elements, dtype=object).ravel())
        fos = list(np.asarray(fo, dtype=object).ravel())
        if len(fos) != len(el_list):
            logger.warning("force_objects length does not match elements; skipping modal load assembly")
            return None
        total = self._total_dof()
        Fe_list = [np.asarray(obj.F_e).ravel() for obj in fos]
        jrd = os.path.join(self.primary_results_dir, "modal_force_assembly")
        return assemble_global_force_vector(el_list, Fe_list, total, job_results_dir=jrd)

    def _compute_secondary_results(self, frequencies, mode_shapes, M_mod=None):
        """
        Modal secondary metrics: generalized mass ``φ_jᵀ M_mod φ_j`` per mode.

        These are not modal participation factors (which require a chosen static load
        pattern); they quantify the norm of each mode in the constrained mass matrix.
        """
        del frequencies  # reserved for frequency-weighted extensions
        phi = np.asarray(mode_shapes, dtype=np.float64)
        if phi.ndim != 2 or M_mod is None:
            self.secondary_results["global"]["modal_generalized_mass"] = np.array([0.0])
            self.secondary_results["global"]["modal_load_participation"] = np.array([0.0])
            logger.info("Secondary modal metrics skipped (missing M_mod or invalid mode_shapes).")
            return
        n_modes = phi.shape[1]
        masses = np.empty(n_modes, dtype=np.float64)
        for j in range(n_modes):
            v = phi[:, j]
            masses[j] = float(np.dot(v, M_mod @ v))

        F_ref = self._assemble_modal_reference_nodal_load()
        participation = np.zeros(n_modes, dtype=np.float64)
        if F_ref is not None and np.max(np.abs(F_ref)) > 1e-18:
            Fv = np.asarray(F_ref, dtype=np.float64).ravel()
            for j in range(n_modes):
                mjj = max(masses[j], 1e-30)
                ph = phi[:, j] / np.sqrt(mjj)
                participation[j] = abs(float(np.dot(ph, Fv)))

        self.secondary_results["global"]["modal_generalized_mass"] = masses
        self.secondary_results["global"]["modal_load_participation"] = participation
        ed = self.settings.get("element_dictionary")
        dof_pn = 7 if ed is not None and mesh_uses_warping_dof(ed) else 6
        mez = modal_effective_mass_fraction_z(M_mod, phi, dof_per_node=dof_pn)
        if mez is not None:
            self.secondary_results["global"]["modal_effective_mass_fraction_z"] = mez
        else:
            self.secondary_results["global"].pop("modal_effective_mass_fraction_z", None)
        logger.info(
            "Computed modal generalized mass and load participation for %s mode(s).", n_modes
        )

    def _save_secondary_results(self):
        """Save secondary modal results."""
        results_dir = os.path.join(self.primary_results_dir, "modal_results")
        os.makedirs(results_dir, exist_ok=True)
        masses = self.secondary_results["global"].get(
            "modal_generalized_mass", np.array([0.0])
        )
        np.savetxt(
            os.path.join(results_dir, f"{self.job_name}_modal_generalized_mass.txt"),
            masses,
            fmt="%.6e",
        )
        part = self.secondary_results["global"].get(
            "modal_load_participation", np.zeros_like(masses)
        )
        np.savetxt(
            os.path.join(results_dir, f"{self.job_name}_modal_load_participation.txt"),
            np.asarray(part, dtype=np.float64).ravel(),
            fmt="%.6e",
        )
        mez = self.secondary_results["global"].get("modal_effective_mass_fraction_z")
        if mez is not None:
            np.savetxt(
                os.path.join(results_dir, f"{self.job_name}_modal_effective_mass_fraction_z.txt"),
                np.asarray(mez, dtype=np.float64).ravel(),
                fmt="%.6e",
            )
        logger.info("Saved secondary modal results (generalized mass and load participation).")

    # -------------------------------------------------------------------------
    # Vibration / buckling analysis bodies (shared with Eigen/Buckling runners)
    # -------------------------------------------------------------------------

    def _run_vibration_analysis(self) -> None:
        eigen_cfg = effective_eigen_config(self.simulation_settings)
        num_modes = int(eigen_cfg["num_modes"])
        dense_th = int(eigen_cfg.get("dense_threshold", 512))
        job_results_dir = self.primary_results_dir
        monitor = RuntimeMonitorTelemetry(job_results_dir=self.diagnostics_dir)
        with monitor.stage("PrepareSpectralLocalMatrices"):
            ke_list, me_list = PrepareSpectralLocalMatrices(
                element_stiffness_matrices=self.element_stiffness_matrices,
                element_mass_matrices=self.element_mass_matrices,
                job_results_dir=job_results_dir,
            ).run()
        with monitor.stage("AssembleSpectralGlobalSystem"):
            K_global, M_global = self._assemble_global_matrices(
                job_results_dir, ke_list=ke_list, me_list=me_list
            )
        with monitor.stage("ModifySpectralGlobalSystem"):
            K_mod, M_mod, _ = self._modify_global_matrices(K_global, M_global, job_results_dir)
        with monitor.stage("SolveGeneralizedEigenproblem"):
            frequencies, mode_shapes = self.solve_modal_vibration(
                K_mod, M_mod, num_modes, job_results_dir, dense_threshold=dense_th
            )
        self._save_primary_results(frequencies, mode_shapes)
        if self._modal_post_processing_enabled():
            U_snap = self._snapshot_u_vibration(mode_shapes)
            self._run_secondary_tertiary_from_formulation_cache(U_snap)
        else:
            self._compute_secondary_results(frequencies, mode_shapes, M_mod=M_mod)
            self._save_secondary_results()
        logger.info("Modal simulation completed successfully -> %s", self.results_root)

    def _run_buckling_analysis(self) -> None:
        buck_cfg = effective_buckling_config(self.simulation_settings)
        num_modes = int(buck_cfg["num_modes"])
        job_results_dir = self.primary_results_dir
        monitor = RuntimeMonitorTelemetry(job_results_dir=self.diagnostics_dir)
        with monitor.stage("PrepareSpectralLocalMatrices"):
            b_ke, b_me = PrepareSpectralLocalMatrices(
                element_stiffness_matrices=self.element_stiffness_matrices,
                element_mass_matrices=self.element_mass_matrices,
                job_results_dir=job_results_dir,
            ).run()
        with monitor.stage("AssembleSpectralGlobalSystem"):
            K_global, M_global = self._assemble_global_matrices(
                job_results_dir, ke_list=b_ke, me_list=b_me
            )
        with monitor.stage("BucklingPrestress"):
            lambdas, mode_shapes = self._run_linear_buckling(
                K_global, job_results_dir, num_modes, buck_cfg, monitor=monitor
            )
        self._save_buckling_results(lambdas, mode_shapes)
        self.secondary_results["global"]["buckling_load_factors"] = lambdas
        if self._modal_post_processing_enabled():
            U_snap = self._snapshot_u_buckling(mode_shapes)
            self._run_secondary_tertiary_from_formulation_cache(U_snap)
        logger.info("Modal buckling simulation completed -> %s", self.results_root)