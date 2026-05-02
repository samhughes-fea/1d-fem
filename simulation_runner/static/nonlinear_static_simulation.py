# simulation_runner/static/nonlinear_static_simulation.py
# Nonlinear (geometric) static finite element analysis (Newton–Raphson iteration).

import logging
import numpy as np
import os
from pathlib import Path
import datetime
from scipy.sparse import coo_matrix

# Diagnostics
from processing.static.diagnostics.runtime_monitor_telemetry import RuntimeMonitorTelemetry

# Operations
from processing.static.operations.preparation import PrepareLocalSystem
from processing.static.operations.assembly import AssembleGlobalSystem
from processing.static.operations.modification import ModifyGlobalSystem
from processing.static.operations.condensation import CondenseModifiedSystem
from processing.static.operations.solver import SolveCondensedSystem
from processing.static.operations.reconstruction import ReconstructGlobalSystem
from processing.static.operations.disassembly import DisassembleGlobalSystem

# Results: containers
from processing.static.results.containers.global_results import GlobalResults
from processing.static.results.containers.elemental_results import ElementalResults
from processing.static.results.containers.nodal_results import NodalResults
from processing.static.results.containers.gaussian_results import GaussianResults
from processing.static.results.containers.container_hopper import PrimaryResultSet, SecondaryResultSet, IndexMapSet

# Results: compute
from processing.static.results.compute_primary.element_formulation_processor import ElementFormulationProcessor
from processing.static.results.compute_primary.primary_results_orchestrator import PrimaryResultsOrchestrator
from processing.static.results.compute_secondary.secondary_results_orchestrator import SecondaryResultsOrchestrator
from processing.static.results.compute_tertiary.tertiary_results_orchestrator import TertiaryResultsOrchestrator

# Results: save
from processing.static.results.save_primary_container import SavePrimaryResults, SavePrimaryResultsSummary
from processing.static.results.save_index_map_container import SaveIndexMaps
from processing.static.results.save_secondary_container import SaveSecondaryResults, SaveSecondaryResultsSummary
from processing.static.results.save_tertiary_container import SaveTertiaryResults, SaveTertiaryResultsSummary
from processing.static.results.build_converged_formulation_cache import build_converged_formulation_cache

logger = logging.getLogger(__name__)


def newton_condensed_residual_converged(
    norm_r_cond: float,
    atol: float,
    rtol: float | None,
    ref_scale: float,
) -> bool:
    """
    Compare ‖F_cond‖₂ (condensed Newton RHS) to atol + rtol·ref_scale.

    ``ref_scale`` is set once from the first Newton iteration (first_residual)
    or from ‖F_ext‖ on condensed DOFs (external_force), depending on settings.

    Parameters
    ----------
    norm_r_cond
        Euclidean norm of the condensed residual vector ``F_cond``.
    atol
        Absolute tolerance on forces (same units as ``F_cond`` entries).
    rtol
        Relative multiplier; if ``None``, only ``atol`` applies via ``atol + 0``.
    ref_scale
        Positive reference scale (first residual or external-load norm).
    """
    rt = 0.0 if rtol is None else float(rtol)
    ref = max(float(ref_scale), np.finfo(float).eps)
    return bool(float(norm_r_cond) <= float(atol) + rt * ref)


def _to_coo(Ke):
    """Convert element stiffness to COO for assembly."""
    if hasattr(Ke, "tocoo"):
        return Ke.tocoo() if not isinstance(Ke, coo_matrix) else Ke
    return coo_matrix(np.asarray(Ke, dtype=np.float64))


class NonlinearStaticSimulationRunner:
    """
    Nonlinear static finite element analysis via Newton–Raphson:
    R = F_ext - F_int(U), solve K_T(U) @ delta_U = R, U += delta_U until convergence.
    Supports mixed meshes (linear and nonlinear elements).
    """

    def __init__(
        self,
        elements,
        grid_dictionary,
        element_dictionary,
        material_dictionary,
        section_dictionary,
        point_load_array,
        distributed_load_array,
        element_objects,
        force_objects,
        job_name,
        job_results_dir,
        simulation_settings=None,
    ):
        from processing.static.results.containers import (
            FormulationResultSet,
            validate_shape_functions_populated,
        )

        self.elements = elements
        self.grid_dictionary = grid_dictionary
        self.element_dictionary = element_dictionary
        self.material_dictionary = material_dictionary
        self.section_dictionary = section_dictionary
        self.point_load_array = point_load_array
        self.distributed_load_array = distributed_load_array
        self.job_name = job_name

        self.formulation_cache = FormulationResultSet(
            element_objects=list(element_objects),
            force_objects=list(force_objects),
        )
        validate_shape_functions_populated(
            self.formulation_cache.element_objects,
            self.formulation_cache.force_objects,
            strict=False,
        )

        # Initial K_e and F_e (nonlinear elements return tangent at U=0 and F_e from loads)
        self.element_stiffness_matrices = np.array([obj.K_e for obj in element_objects], dtype=object)
        self.element_force_vectors = np.array([obj.F_e for obj in force_objects], dtype=object)

        self.results_root = job_results_dir
        self.primary_results_dir = os.path.join(self.results_root, "primary_results")
        self.secondary_results_dir = os.path.join(self.results_root, "secondary_results")
        self.tertiary_results_dir = os.path.join(self.results_root, "tertiary_results")
        self.maps_dir = os.path.join(self.results_root, "maps")
        self.logs_dir = os.path.join(self.results_root, "logs")
        self.diagnostics_dir = os.path.join(self.results_root, "diagnostics")

        self.monitor = RuntimeMonitorTelemetry(job_results_dir=self.primary_results_dir)

        self.primary_results_set = PrimaryResultSet(
            gaussian_results=GaussianResults(),
            nodal_results=NodalResults(),
            elemental_results=ElementalResults(),
            global_results=GlobalResults(),
        )
        self.secondary_results_set = SecondaryResultSet(
            gaussian_results=GaussianResults(),
            nodal_results=NodalResults(),
            elemental_results=ElementalResults(),
            global_results=GlobalResults(),
        )
        self.maps = IndexMapSet()

        # Validate mesh and element data
        if len(self.elements) == 0 or not self.element_dictionary or not self.grid_dictionary:
            logger.error("❌ Error: Missing elements or mesh data in constructor!")
            raise ValueError("❌ Error: Missing elements or mesh data!")

        if "ids" in self.grid_dictionary:
            self.node_ids = self.grid_dictionary["ids"]
        elif "nodes" in self.grid_dictionary and "ids" in self.grid_dictionary["nodes"]:
            self.node_ids = self.grid_dictionary["nodes"]["ids"]
        else:
            raise KeyError("grid_dictionary must contain an 'ids' array of node IDs")

        self.total_dof = len(self.node_ids) * 6
        logger.debug("Total DOFs initialised: %d (nodes %d × 6)", self.total_dof, len(self.node_ids))
# -----------------------------------------------------------------------

        simulation_settings = simulation_settings or {}
        self.solver_config = simulation_settings.get("solver", {
            "type": "cg",
            "tolerance": 1e-6,
            "max_iterations": 1000,
            "restart": 20,
            "ilu_drop_tol": 1e-6,
            "ilu_fill_factor": 1.0,
            "disable_scaling": False,
        })
        self.condensation_config = simulation_settings.get("condensation", {"base_tol": 1e-12})
        self.integration_config = simulation_settings.get("integration", {})
        self.output_config = simulation_settings.get("output", {})

        # Newton–Raphson settings
        newton = simulation_settings.get("newton", {})
        self.newton_tol = float(newton.get("tolerance", 1e-8))
        self.newton_max_iter = int(newton.get("max_iterations", 50))
        self.newton_tol_delta_u = float(newton.get("tolerance_delta_u", 1e-10))
        rel = newton.get("relative_tolerance", None)
        self.newton_relative_tol = None if rel is None else float(rel)
        ref_mode = str(newton.get("relative_reference", "first_residual")).lower().strip()
        if ref_mode not in ("first_residual", "external_force"):
            raise ValueError(
                "newton.relative_reference must be 'first_residual' or 'external_force', "
                f"got {ref_mode!r}",
            )
        self.newton_relative_reference = ref_mode

        self.prescribed_displacements = None

        # Intermediate system storage
        self.K_global = None
        self.F_global = None
        self.K_mod = None
        self.F_mod = None
        self.fixed_dofs = None
        self.condensed_dofs = None
        self.inactive_dofs = None
        self.K_cond = None
        self.F_cond = None
        self.U_cond = None
        self.U_global = None
        self.local_global_dof_map = None

    def setup_simulation(self):
        """Create the necessary output directory structure."""
        logger.info("✅ Setting up nonlinear static simulation for job: %s...", self.job_name)
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.secondary_results_dir, exist_ok=True)
        os.makedirs(self.tertiary_results_dir, exist_ok=True)
        os.makedirs(self.maps_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        os.makedirs(self.diagnostics_dir, exist_ok=True)
        logger.info("📁 Results will be saved to: %s", self.results_root)

    # -------------------------------------------------------------------------
    # 0) PREPARE LOCAL SYSTEM
    # -------------------------------------------------------------------------

    def prepare_local_system(self, job_results_dir: str):
        """
        Validate and format local element stiffness matrices (Ke) and force vectors (Fe)
        using the PrepareLocalSystem helper for assembly.

        Parameters
        ----------
        job_results_dir : str
            Directory for logging local system processing.

        Returns
        -------
        element_stiffness_matrices : List[scipy.sparse.coo_matrix]
        element_force_vectors : List[numpy.ndarray]
        """
        logger.info("📦 Preparing local element stiffness matrices and force vectors...")

        preparer = PrepareLocalSystem(
            Ke_raw=self.element_stiffness_matrices,
            Fe_raw=self.element_force_vectors,
            job_results_dir=job_results_dir,
        )

        try:
            Ke_formatted, Fe_formatted = preparer.validate_and_format()
        except Exception as exc:
            logger.error("❌ Local system preparation failed – see prepare_local_system.log for details")
            raise

        self.element_stiffness_matrices = Ke_formatted
        self.element_force_vectors = Fe_formatted

        logger.info("✅ Local element systems validated and formatted for assembly")
        return Ke_formatted, Fe_formatted

    # -------------------------------------------------------------------------
    # 1) ASSEMBLE GLOBAL SYSTEM
    # -------------------------------------------------------------------------

    def assemble_global_system(
        self,
        job_results_dir: str,
        *,
        element_stiffness_matrices=None,
        element_force_vectors=None,
    ):
        """
        Assemble the global stiffness matrix and force vector from element data.

        Parameters
        ----------
        job_results_dir : str
            Directory for assembly logs.
        element_stiffness_matrices : array-like, optional
            Per-element stiffness (default: self.element_stiffness_matrices).
        element_force_vectors : array-like, optional
            Per-element force vectors (default: self.element_force_vectors).

        Returns
        -------
        K_global : scipy.sparse matrix
        F_global : np.ndarray
        local_global_dof_map : list
        assembly_map : dict
        """
        logger.info("🔧 Assembling global stiffness and force matrices...")

        Ke = self.element_stiffness_matrices if element_stiffness_matrices is None else element_stiffness_matrices
        Fe = self.element_force_vectors if element_force_vectors is None else element_force_vectors

        assembler = AssembleGlobalSystem(
            elements=self.elements,
            element_stiffness_matrices=Ke,
            element_force_vectors=Fe,
            total_dof=self.total_dof,
            job_results_dir=job_results_dir,
        )

        try:
            K_global, F_global, local_global_dof_map, assembly_map = assembler.assemble()
        except Exception as exc:
            logger.error("❌ Assembly failed – see assembly.log for details")
            raise

        F_global = np.asarray(F_global, dtype=np.float64).ravel()

        self.K_global = K_global
        self.F_global = F_global
        self.local_global_dof_map = local_global_dof_map
        self.maps.assembly_map = assembly_map
        self.primary_results_set.global_results.F_global = self.F_global
        self.primary_results_set.global_results.K_global = self.K_global

        logger.info("✅ Global stiffness matrix and force vector successfully assembled")
        return K_global, F_global, local_global_dof_map, assembly_map

    # -------------------------------------------------------------------------
    # 2) MODIFY GLOBAL SYSTEM (BOUNDARY CONDITIONS)
    # -------------------------------------------------------------------------

    def modify_global_system(
        self,
        K_global,
        F_global,
        local_global_dof_map,
        job_results_dir: str,
        *,
        fixed_dofs: np.ndarray | None = None,
        prescribed_displacements: dict | None = None,
    ):
        """
        Apply boundary conditions to the global system.

        Parameters
        ----------
        K_global : scipy.sparse matrix
            Global stiffness matrix.
        F_global : np.ndarray
            Global force vector.
        local_global_dof_map : list
            Local to global DOF mapping.
        job_results_dir : str
            Results directory.
        fixed_dofs : np.ndarray | None, optional
            Fixed DOF indices.
        prescribed_displacements : dict | None, optional
            Dictionary with 'global_dof' and 'value' arrays.

        Returns
        -------
        K_mod : scipy.sparse matrix
        F_mod : np.ndarray
        fixed_dofs : np.ndarray
        modification_map : dict
        """
        logger.info("🔒 Applying boundary conditions to global matrices…")

        modifier = ModifyGlobalSystem(
            K_global=K_global,
            F_global=np.asarray(F_global).ravel(),
            job_results_dir=job_results_dir,
            fixed_dofs=fixed_dofs,
            prescribed_displacements=prescribed_displacements,
            local_global_dof_map=local_global_dof_map,
        )

        try:
            K_mod, F_mod, fixed_dofs, modification_map = modifier.apply_boundary_conditions()
        except Exception:
            logger.error("❌ Boundary-condition step failed – see ModifyGlobalSystem.log")
            raise

        self.K_mod = K_mod
        self.F_mod = F_mod
        self.fixed_dofs = fixed_dofs
        self.maps.modification_map = modification_map
        self.primary_results_set.global_results.K_mod = K_mod
        self.primary_results_set.global_results.F_mod = F_mod

        logger.info("✅ Boundary conditions successfully applied")
        return K_mod, F_mod, fixed_dofs, modification_map

    # -------------------------------------------------------------------------
    # 3) CONDENSE MODIFIED SYSTEM
    # -------------------------------------------------------------------------

    def condense_modified_system(
        self,
        K_mod,
        F_mod,
        fixed_dofs,
        local_global_dof_map,
        job_results_dir: str,
        *,
        base_tol: float = 1e-12,
    ):
        """
        Perform static condensation on the boundary-conditioned system.

        Parameters
        ----------
        K_mod : scipy.sparse matrix
            Stiffness matrix after boundary conditions.
        F_mod : np.ndarray
            Force vector (flattened).
        fixed_dofs : np.ndarray
            Fixed DOF indices.
        local_global_dof_map : list
            Local-to-global DOF map.
        job_results_dir : str
            Directory for condensation logs.
        base_tol : float, optional
            Baseline tolerance for condensation.

        Returns
        -------
        condensed_dofs : np.ndarray
        inactive_dofs : np.ndarray
        K_cond : scipy.sparse matrix
        F_cond : np.ndarray
        condensation_map : dict
        """
        logger.info("📉 Condensing modified system…")

        condenser = CondenseModifiedSystem(
            K_mod=K_mod,
            F_mod=np.asarray(F_mod).ravel(),
            fixed_dofs=fixed_dofs,
            local_global_dof_map=local_global_dof_map,
            job_results_dir=job_results_dir,
            base_tol=base_tol,
        )

        try:
            condensed_dofs, inactive_dofs, K_cond, F_cond, cond_map = condenser.apply_condensation()
        except Exception as exc:
            logger.error("❌ Condensation failed – see condensation.log for details")
            raise

        self.K_cond = K_cond
        self.F_cond = F_cond
        self.inactive_dofs = inactive_dofs
        self.condensed_dofs = condensed_dofs
        self.maps.condensation_map = cond_map
        self.primary_results_set.global_results.K_cond = K_cond
        self.primary_results_set.global_results.F_cond = F_cond

        logger.info(
            "✅ Static condensation complete "
            "(%d DOFs retained, %d pruned)",
            len(condensed_dofs), len(inactive_dofs),
        )
        return condensed_dofs, inactive_dofs, K_cond, F_cond, cond_map

    # -------------------------------------------------------------------------
    # 4) SOLVE CONDENSED SYSTEM
    # -------------------------------------------------------------------------

    def solve_condensed_system(
        self,
        K_cond,
        F_cond,
        job_results_dir: str,
        *,
        solver_name: str = "cg",
        preconditioner: str | None = "auto",
        tolerance: float = 1e-6,
        max_iterations: int = 1000,
        restart: int = 20,
        ilu_drop_tol: float = 1e-6,
        ilu_fill_factor: float = 1.0,
        disable_scaling: bool = False,
    ):
        """
        Solve the condensed linear system K_cond · U_cond = F_cond.

        Parameters
        ----------
        K_cond, F_cond : scipy.sparse matrix, np.ndarray
            Condensed stiffness and force vector.
        job_results_dir : str
            Folder for solver logs and output.
        solver_name : str, optional
            Solver type (default "cg").
        preconditioner : str | None, optional
            Preconditioning strategy (default "auto").
        tolerance : float, optional
            Solver tolerance.
        max_iterations : int, optional
            Maximum iterations.
        restart : int, optional
            GMRES restart.
        ilu_drop_tol, ilu_fill_factor : float, optional
            ILU options.
        disable_scaling : bool, optional
            Disable scaling.

        Returns
        -------
        U_cond : np.ndarray
        """
        logger.info("🧮 Solving condensed system…")

        solver = SolveCondensedSystem(
            K_cond=K_cond,
            F_cond=F_cond,
            solver_name=solver_name,
            job_results_dir=job_results_dir,
            preconditioner=preconditioner,
            tolerance=tolerance,
            max_iterations=max_iterations,
            restart=restart,
            ilu_drop_tol=ilu_drop_tol,
            ilu_fill_factor=ilu_fill_factor,
            disable_scaling=disable_scaling,
        )

        U_cond = solver.solve()
        if U_cond is None:
            raise RuntimeError("Condensed solver failed – see SolveCondensedSystem.log")

        self.U_cond = U_cond
        self.primary_results_set.global_results.U_cond = U_cond

        logger.info("✅ Condensed system successfully solved")
        return U_cond

    # -------------------------------------------------------------------------
    # 5) RECONSTRUCT GLOBAL SYSTEM
    # -------------------------------------------------------------------------

    def reconstruct_global_system(
        self,
        condensed_dofs: np.ndarray,
        U_cond: np.ndarray,
        total_dof: int,
        job_results_dir: str | Path,
        *,
        fixed_dofs: np.ndarray | None = None,
        inactive_dofs: np.ndarray | None = None,
        local_global_dof_map=None,
    ) -> np.ndarray:
        """
        Reconstruct the full displacement vector from the condensed solution.

        Parameters
        ----------
        condensed_dofs : np.ndarray
            Active DOF indices in original numbering.
        U_cond : np.ndarray
            Displacements in condensed ordering.
        total_dof : int
            Total DOFs in the original system.
        job_results_dir : str | Path
            Directory for reconstruction logs.
        fixed_dofs : np.ndarray | None, optional
            Fixed DOF indices.
        inactive_dofs : np.ndarray | None, optional
            Inactive DOFs removed during condensation.
        local_global_dof_map : list, optional
            Local-to-global DOF map (default: self.local_global_dof_map).

        Returns
        -------
        U_global : np.ndarray
        """
        logger.info("🔄 Reconstructing global displacement vector…")

        dof_map = local_global_dof_map if local_global_dof_map is not None else self.local_global_dof_map
        reconstructor = ReconstructGlobalSystem(
            active_dofs=condensed_dofs,
            U_cond=U_cond,
            total_dofs=total_dof,
            job_results_dir=Path(job_results_dir),
            fixed_dofs=fixed_dofs,
            inactive_dofs=inactive_dofs,
            local_global_dof_map=dof_map,
        )

        try:
            U_global, reconstruction_map = reconstructor.reconstruct()
        except Exception:
            logger.error("❌ Reconstruction failed – see ReconstructGlobalSystem.log")
            raise

        self.U_global = U_global
        self.primary_results_set.global_results.U_global = U_global
        self.maps.reconstruction_map = reconstruction_map

        logger.info(
            "✅ Displacement reconstruction complete "
            "(min=%.3e, max=%.3e)", U_global.min(), U_global.max(),
        )
        return U_global

    # -------------------------------------------------------------------------
    # 6) PRIMARY RESULTS PIPELINE
    # -------------------------------------------------------------------------

    def compute_primary_results(self) -> None:
        """
        Build global-level primary outputs (reactions), disassemble U_global and R_global
        onto elements, cache results, and save CSVs.

        For nonlinear runs, element_stiffness_matrices and element_force_vectors are
        expected to be the converged formulation (K_T, F_int) after BuildConvergedFormulationCache.
        """
        logger.info("📊 Computing primary results…")

        # Step 1: Validate and standardise element matrices (COO → CSR)
        self.element_stiffness_matrices, self.element_force_vectors = ElementFormulationProcessor(
            F_e=self.element_force_vectors,
            K_e=self.element_stiffness_matrices,
        ).process()

        # Step 2: Compute reactions and residuals
        self.R_global, self.R_residual = PrimaryResultsOrchestrator(
            K_global=self.K_global,
            F_global=self.F_global,
            U_global=self.U_global,
            fixed_dofs=self.fixed_dofs,
            job_results_dir=self.primary_results_dir,
        ).compute()

        # Step 3: Disassemble element-level primary results
        self.R_e, self.R_residual_e, self.U_e = DisassembleGlobalSystem(
            U_global=self.U_global,
            R_global=self.R_global,
            R_residual=self.R_residual,
            F_e=self.element_force_vectors,
            local_global_dof_map=self.local_global_dof_map,
            job_results_dir=self.primary_results_dir,
        ).disassemble()

        # Step 4: Construct global results
        global_results = GlobalResults(
            F_global=self.F_global,
            K_global=self.K_global,
            F_mod=self.F_mod,
            K_mod=self.K_mod,
            F_cond=self.F_cond,
            K_cond=self.K_cond,
            U_cond=self.U_cond,
            U_global=self.U_global,
            R_global=self.R_global,
            R_residual=self.R_residual,
        )

        # Step 5: Construct element results
        elemental_results = ElementalResults(
            K_e=self.element_stiffness_matrices,
            F_e=self.element_force_vectors,
            U_e=self.U_e,
            R_e=self.R_e,
            R_residual_e=self.R_residual_e,
        )

        # Step 6: Cache results
        self.primary_results = PrimaryResultSet(
            global_results=global_results,
            elemental_results=elemental_results,
        )

        # Step 7: Save CSVs
        SavePrimaryResults(
            primary_results_set=self.primary_results,
            index_map_set=self.maps,
            save_dir=self.primary_results_dir,
        ).run()
        SavePrimaryResultsSummary(primary_results_set=self.primary_results, save_dir=self.primary_results_dir).save()
        SaveIndexMaps(index_map_set=self.maps, save_dir=self.maps_dir).run()

        logger.info("✅ Primary results computed, written and saved")

    # -------------------------------------------------------------------------
    # 7) SECONDARY RESULTS PIPELINE
    # -------------------------------------------------------------------------

    def compute_secondary_results(self) -> None:
        """Compute and save secondary results (stresses, strains, etc.)."""
        logger.info("📈 Computing secondary results using cached formulation data...")

        orchestrator = SecondaryResultsOrchestrator(
            elements=self.elements,
            grid_dictionary=self.grid_dictionary,
            element_dictionary=self.element_dictionary,
            material_dictionary=self.material_dictionary,
            section_dictionary=self.section_dictionary,
            global_displacement=self.U_global,
            formulation_cache=self.formulation_cache,
            job_results_dir=self.secondary_results_dir,
        )
        self.secondary_results_set = orchestrator.compute_all()

        SaveSecondaryResults(
            secondary_results=self.secondary_results_set,
            save_dir=self.secondary_results_dir,
        ).save_all()
        SaveSecondaryResultsSummary(
            secondary_results=self.secondary_results_set,
            save_dir=self.secondary_results_dir,
        ).save()

        logger.info("✅ Secondary results computed and saved")

    # -------------------------------------------------------------------------
    # 8) TERTIARY RESULTS PIPELINE
    # -------------------------------------------------------------------------

    def compute_tertiary_results(self) -> None:
        """Compute and save tertiary results (section forces, etc.)."""
        logger.info("📊 Computing tertiary results...")

        orchestrator = TertiaryResultsOrchestrator(
            secondary_results=self.secondary_results_set,
            formulation_cache=self.formulation_cache,
            element_dictionary=self.element_dictionary,
            grid_dictionary=self.grid_dictionary,
            job_results_dir=self.tertiary_results_dir,
        )
        self.tertiary_results = orchestrator.compute()

        SaveTertiaryResults(tertiary_results=self.tertiary_results, save_dir=self.results_root).save_all()
        SaveTertiaryResultsSummary(tertiary_results=self.tertiary_results, save_dir=self.results_root).save()

        logger.info("✅ Tertiary results computed and saved")

    # -------------------------------------------------------------------------
    # NEWTON HELPERS
    # -------------------------------------------------------------------------

    def _build_K_T_and_F_int(self, U_global: np.ndarray):
        """Build lists of tangent stiffness and internal force for current U (for Newton step)."""
        K_T_list = []
        F_int_list = []
        for i, elem in enumerate(self.elements):
            dof_map = elem.assemble_global_dof_indices()
            U_e = np.asarray(U_global[dof_map], dtype=np.float64).ravel()
            if hasattr(elem, "tangent_stiffness_matrix") and hasattr(elem, "internal_force_vector"):
                K_T_list.append(_to_coo(elem.tangent_stiffness_matrix(U_e)))
                F_int_list.append(np.asarray(elem.internal_force_vector(U_e), dtype=np.float64).ravel())
            else:
                Ke = self.element_stiffness_matrices[i]
                K_T_list.append(_to_coo(Ke))
                Ke_dense = Ke.toarray() if hasattr(Ke, "toarray") else np.asarray(Ke)
                F_int_list.append(Ke_dense @ U_e)
        return K_T_list, F_int_list

    # -------------------------------------------------------------------------
    # TOP-LEVEL SIMULATION FLOW
    # -------------------------------------------------------------------------

    def run(self) -> None:
        """Execute the complete nonlinear static workflow (Newton–Raphson then primary/secondary/tertiary)."""
        try:
            self.monitor = RuntimeMonitorTelemetry(job_results_dir=self.diagnostics_dir)

            # -----------------------------------------------------------------
            # 0  Housekeeping and local systems
            # -----------------------------------------------------------------
            self.setup_simulation()
            self.start_time = datetime.datetime.now()

            with self.monitor.stage("PrepareLocalSystem"):
                self.prepare_local_system(job_results_dir=self.primary_results_dir)

            # -----------------------------------------------------------------
            # Initial assembly (F_ext and DOF map; initial K not used in Newton loop)
            # -----------------------------------------------------------------
            with self.monitor.stage("AssembleInitial"):
                self.assemble_global_system(self.primary_results_dir)
                F_ext_global = np.asarray(self.F_global, dtype=np.float64).ravel()

            # -----------------------------------------------------------------
            # Newton–Raphson loop
            # -----------------------------------------------------------------
            U_global = np.zeros(self.total_dof, dtype=np.float64)
            if self.prescribed_displacements is not None:
                gd = np.asarray(self.prescribed_displacements["global_dof"], dtype=np.int32)
                val = np.asarray(self.prescribed_displacements["value"], dtype=np.float64)
                U_global[gd] = val

            newton_tol = self.newton_tol
            newton_max = self.newton_max_iter
            tol_du = self.newton_tol_delta_u
            delta_U_cond = None
            newton_ref_residual = None
            newton_ref_external = None

            for iteration in range(newton_max):
                with self.monitor.stage("NewtonBuildKT_Fint"):
                    K_T_list, F_int_list = self._build_K_T_and_F_int(U_global)

                with self.monitor.stage("NewtonAssemble"):
                    self.assemble_global_system(
                        self.primary_results_dir,
                        element_stiffness_matrices=K_T_list,
                        element_force_vectors=F_int_list,
                    )
                    F_int_global = np.asarray(self.F_global, dtype=np.float64).ravel()

                R_global = F_ext_global - F_int_global

                with self.monitor.stage("NewtonModify"):
                    self.modify_global_system(
                        self.K_global,
                        R_global,
                        self.local_global_dof_map,
                        self.primary_results_dir,
                        prescribed_displacements=getattr(self, "prescribed_displacements", None),
                    )
                    R_mod = self.F_mod

                with self.monitor.stage("NewtonCondense"):
                    self.condense_modified_system(
                        self.K_mod,
                        np.asarray(R_mod).ravel(),
                        self.fixed_dofs,
                        self.local_global_dof_map,
                        self.primary_results_dir,
                        base_tol=self.condensation_config.get("base_tol", 1e-12),
                    )

                F_cond_arr = np.asarray(self.F_cond, dtype=np.float64).ravel()
                norm_R_conv = float(np.linalg.norm(F_cond_arr))
                norm_R_full = float(np.linalg.norm(R_global))
                if iteration == 0:
                    newton_ref_residual = norm_R_conv
                    cd0 = np.asarray(self.condensed_dofs, dtype=np.intp)
                    newton_ref_external = float(np.linalg.norm(F_ext_global[cd0]))

                if self.newton_relative_reference == "external_force":
                    ref_scale = newton_ref_external
                else:
                    ref_scale = newton_ref_residual
                residual_ok = newton_condensed_residual_converged(
                    norm_R_conv,
                    newton_tol,
                    self.newton_relative_tol,
                    ref_scale,
                )

                with self.monitor.stage("NewtonSolve"):
                    delta_U_cond = self.solve_condensed_system(
                        self.K_cond,
                        self.F_cond,
                        self.primary_results_dir,
                        solver_name=self.solver_config.get("type", "cg"),
                        preconditioner="auto",
                        tolerance=self.solver_config.get("tolerance", 1e-6),
                        max_iterations=self.solver_config.get("max_iterations", 1000),
                        restart=self.solver_config.get("restart", 20),
                        ilu_drop_tol=self.solver_config.get("ilu_drop_tol", 1e-6),
                        ilu_fill_factor=self.solver_config.get("ilu_fill_factor", 1.0),
                        disable_scaling=self.solver_config.get("disable_scaling", False),
                    )
                    if delta_U_cond is None:
                        raise RuntimeError("Newton step solve failed")

                with self.monitor.stage("NewtonReconstruct"):
                    delta_U_global = self.reconstruct_global_system(
                        self.condensed_dofs,
                        delta_U_cond,
                        self.total_dof,
                        self.primary_results_dir,
                        fixed_dofs=self.fixed_dofs,
                        inactive_dofs=self.inactive_dofs,
                        local_global_dof_map=self.local_global_dof_map,
                    )
                U_global += delta_U_global

                norm_du = float(np.linalg.norm(delta_U_global))
                thr = newton_tol
                if self.newton_relative_tol is not None:
                    rs = max(float(ref_scale), np.finfo(float).eps)
                    thr = newton_tol + self.newton_relative_tol * rs
                logger.info(
                    "Newton iter %d: ||R||_full=%.4e ||F_cond||=%.4e (vs thr=%.4e) ||delta_U||=%.4e",
                    iteration + 1,
                    norm_R_full,
                    norm_R_conv,
                    thr,
                    norm_du,
                )
                if residual_ok and norm_du < tol_du:
                    logger.info("Newton converged at iteration %d", iteration + 1)
                    break
            else:
                logger.warning("Newton did not converge within %d iterations", newton_max)

            self.U_global = U_global
            self.F_global = F_ext_global
            self.primary_results_set.global_results.U_global = self.U_global
            self.primary_results_set.global_results.K_global = self.K_global
            self.primary_results_set.global_results.F_global = self.F_global
            self.primary_results_set.global_results.K_mod = self.K_mod
            self.primary_results_set.global_results.F_mod = self.F_mod
            self.primary_results_set.global_results.K_cond = self.K_cond
            self.primary_results_set.global_results.F_cond = self.F_cond
            self.primary_results_set.global_results.U_cond = delta_U_cond

            # -----------------------------------------------------------------
            # 6  Build converged formulation cache (Option B)
            # -----------------------------------------------------------------
            with self.monitor.stage("BuildConvergedFormulationCache"):
                converged_cache = build_converged_formulation_cache(
                    elements=self.elements,
                    reference_cache=self.formulation_cache,
                    U_global=self.U_global,
                )
                self.formulation_cache = converged_cache
                self.element_stiffness_matrices = np.array(
                    [obj.K_e for obj in converged_cache.element_objects], dtype=object
                )
                self.element_force_vectors = np.array(
                    [obj.F_e for obj in converged_cache.force_objects], dtype=object
                )

            # -----------------------------------------------------------------
            # 7  Primary results (disassembly and CSV export)
            # -----------------------------------------------------------------
            with self.monitor.stage("ComputePrimaryResults"):
                self.compute_primary_results()

            # -----------------------------------------------------------------
            # 8  Save formulation cache
            # -----------------------------------------------------------------
            with self.monitor.stage("SaveFormulationCache"):
                from processing.static.results.save_formulation_container import SaveFormulationData
                SaveFormulationData(
                    formulation_cache=self.formulation_cache,
                    save_dir=self.primary_results_dir,
                    save_gauss_data=False,
                ).save_all()

            # -----------------------------------------------------------------
            # 9  Secondary (derived) results
            # -----------------------------------------------------------------
            with self.monitor.stage("ComputeSecondaryResults"):
                self.compute_secondary_results()

            # -----------------------------------------------------------------
            # 10  Tertiary (design/verification) results
            # -----------------------------------------------------------------
            with self.monitor.stage("ComputeTertiaryResults"):
                self.compute_tertiary_results()

            logger.info("🏁 Nonlinear static simulation completed successfully → %s", self.results_root)

        except Exception as exc:
            logger.exception("💥 Simulation failed with critical error")
            raise RuntimeError("Simulation aborted") from exc