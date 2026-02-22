# simulation_runner/static_simulation.py

import logging
import numpy as np
import os
from pathlib import Path
import datetime

# Diagnostics 
from processing_OOP.static.diagnostics.linear_static_diagnostic import DiagnoseLinearStaticSystem
from processing_OOP.static.diagnostics.runtime_monitor_telemetry import RuntimeMonitorTelemetry

# Operations
from processing_OOP.static.operations.preparation import PrepareLocalSystem
from processing_OOP.static.operations.assembly import AssembleGlobalSystem
from processing_OOP.static.operations.modification import ModifyGlobalSystem
from processing_OOP.static.operations.condensation import CondenseModifiedSystem
from processing_OOP.static.operations.solver import SolveCondensedSystem
from processing_OOP.static.operations.reconstruction import ReconstructGlobalSystem
from processing_OOP.static.operations.disassembly import DisassembleGlobalSystem

# Results containers

from processing_OOP.static.results.containers.global_results import GlobalResults
from processing_OOP.static.results.containers.elemental_results import ElementalResults
from processing_OOP.static.results.containers.nodal_results import NodalResults
from processing_OOP.static.results.containers.gaussian_results import GaussianResults


from processing_OOP.static.results.compute_primary.element_formulation_processor import ElementFormulationProcessor
from processing_OOP.static.results.compute_primary.primary_results_orchestrator import PrimaryResultsOrchestrator
from processing_OOP.static.results.containers.container_hopper import PrimaryResultSet, SecondaryResultSet, IndexMapSet

from processing_OOP.static.results.save_primary_container import SavePrimaryResults, SavePrimaryResultsSummary
from processing_OOP.static.results.save_index_map_container import SaveIndexMaps

# Configure module-level logger
logger = logging.getLogger(__name__)


class StaticSimulationRunner:
    """
    Handles static finite element analysis.
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
        element_objects,           # NEW: ElementObject[]
        force_objects,             # NEW: ForceObject[]
        job_name,
        job_results_dir,
        simulation_settings=None   # NEW: Optional simulation settings dict
    ):
        from processing_OOP.static.results.containers import (
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

        # Create FormulationResultSet from objects
        self.formulation_cache = FormulationResultSet(
            element_objects=list(element_objects),
            force_objects=list(force_objects)
        )
        validate_shape_functions_populated(
            self.formulation_cache.element_objects,
            self.formulation_cache.force_objects,
            strict=False,
        )

        # Extract matrices/vectors for backward compatibility with assembly
        self.element_stiffness_matrices = np.array([obj.K_e for obj in element_objects], dtype=object)
        self.element_force_vectors = np.array([obj.F_e for obj in force_objects], dtype=object)

        self.results_root = job_results_dir  # Use exact path passed from run_job.py
        self.primary_results_dir = os.path.join(self.results_root, "primary_results")
        self.secondary_results_dir = os.path.join(self.results_root, "secondary_results")
        self.tertiary_results_dir = os.path.join(self.results_root, "tertiary_results")
        self.maps_dir = os.path.join(self.results_root, "maps")
        self.logs_dir = os.path.join(self.results_root, "logs")
        self.diagnostics_dir = os.path.join(self.results_root, "diagnostics")

        self.monitor = RuntimeMonitorTelemetry(job_results_dir=self.primary_results_dir)

        # Initialize results storage
        self.primary_results_set = PrimaryResultSet(
            gaussian_results = GaussianResults(),
            nodal_results    = NodalResults(),
            elemental_results  = ElementalResults(),
            global_results   = GlobalResults(),
        )

        # Initialize results storage
        self.secondary_results_set = SecondaryResultSet(
            gaussian_results = GaussianResults(),
            nodal_results    = NodalResults(),
            elemental_results  = ElementalResults(),
            global_results   = GlobalResults(),
        )

        # Initialise mapping storage
        self.maps = IndexMapSet()


        # Validate mesh and element data
        if len(self.elements) == 0 or not self.element_dictionary or not self.grid_dictionary:
            logger.error("❌ Error: Missing elements or mesh data in constructor!")
            raise ValueError("❌ Error: Missing elements or mesh data!")
        
        # ---------------- total DOF initialisation ------------------
        # GridParser outputs a flat key: grid_dictionary["ids"]
        if "ids" in self.grid_dictionary:
            self.node_ids = self.grid_dictionary["ids"]
        elif "nodes" in self.grid_dictionary and "ids" in self.grid_dictionary["nodes"]:
            # legacy nested layout still supported
            self.node_ids = self.grid_dictionary["nodes"]["ids"]
        else:
            raise KeyError("grid_dictionary must contain an 'ids' array of node IDs")

        self.total_dof = len(self.node_ids) * 6
        logger.debug("Total DOFs initialised: %d (nodes %d × 6)", self.total_dof, len(self.node_ids))
# -----------------------------------------------------------------------

        # Extract simulation settings with defaults
        if simulation_settings is None:
            simulation_settings = {}
        
        # Extract solver configuration with defaults
        self.solver_config = simulation_settings.get("solver", {
            "type": "cg",
            "tolerance": 1e-6,
            "max_iterations": 1000,
            "restart": 20,
            "ilu_drop_tol": 1e-6,
            "ilu_fill_factor": 1.0,
            "disable_scaling": False,
        })
        
        # Extract condensation configuration with defaults
        self.condensation_config = simulation_settings.get("condensation", {
            "base_tol": 1e-12,
        })
        
        # Extract integration configuration with defaults (for future use)
        self.integration_config = simulation_settings.get("integration", {})
        
        # Extract output configuration with defaults (for future use)
        self.output_config = simulation_settings.get("output", {})
        
        # General simulation settings (optional or future expansion)
        self.solver_name = self.solver_config.get("type", "cg")
        
        # Prescribed displacements (set externally via setter or passed in)
        self.prescribed_displacements = None

        # Initialize intermediate system storage
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

    def setup_simulation(self):
        """Creates the necessary output directory structure."""
        logger.info(f"✅ Setting up static simulation for job: {self.job_name}...")
        os.makedirs(self.primary_results_dir, exist_ok=True)
        os.makedirs(self.secondary_results_dir, exist_ok=True)
        os.makedirs(self.tertiary_results_dir, exist_ok=True)
        os.makedirs(self.maps_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        os.makedirs(self.diagnostics_dir, exist_ok=True)
        logger.info(f"📁 Results will be saved to: {self.results_root}")

    # -------------------------------------------------------------------------
    # 0) PREPARE LOCAL SYSTEM
    # -------------------------------------------------------------------------

    def prepare_local_system(self, job_results_dir: str):
        """
        Validate and format local element stiffness matrices (Ke) and force vectors (Fe)
        using the PrepareLocalSystem helper to format the data structures ready for assembly.

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
            job_results_dir=job_results_dir
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

    def assemble_global_system(self, job_results_dir: str):
        """
        Assemble the global stiffness matrix (K_global) and force vector (F_global),
        and cache them immediately in self.global_results.
        """

        logger.info("🔧 Assembling global stiffness and force matrices...")

        assembler = AssembleGlobalSystem(
            elements=self.elements,
            element_stiffness_matrices=self.element_stiffness_matrices,
            element_force_vectors=self.element_force_vectors,
            total_dof=self.total_dof,
            job_results_dir=job_results_dir
        )

        try:
            K_global, F_global, local_global_dof_map, assembly_map = assembler.assemble()
        except Exception as exc:
            logger.error("❌ Assembly failed – see assembly.log for details")
            raise

        F_global = np.asarray(F_global, dtype=np.float64).ravel()

        # Cache to self
        self.K_global = K_global
        self.F_global = F_global
        self.local_global_dof_map = local_global_dof_map

        # Cache map and results to simulation runner results container
        self.maps.assembly_map = assembly_map
        self.primary_results_set.global_results.F_global = self.F_global
        self.primary_results_set.global_results.K_global = self.K_global

        logger.info("✅ Global stiffness matrix and force vector successfully assembled")

        return K_global, F_global, local_global_dof_map

    # -------------------------------------------------------------------------
    # 2) MODIFY GLOBAL MATRICES (BOUNDARY CONDITIONS)
    # -------------------------------------------------------------------------
    def modify_global_system(
            self,
            K_global,
            F_global,
            local_global_dof_map,
            job_results_dir: str,
            *,
            fixed_dofs: np.ndarray | None = None,
            prescribed_displacements: dict | None = None):
        """
        Apply boundary conditions to the global system.

        Parameters
        ----------
        K_global : csr_matrix
            Global stiffness matrix
        F_global : np.ndarray
            Global force vector
        local_global_dof_map : list
            Local to global DOF mapping
        job_results_dir : str
            Results directory
        fixed_dofs : np.ndarray | None, optional
            Fixed DOF indices (for backward compatibility)
        prescribed_displacements : dict | None, optional
            Dictionary with 'global_dof' and 'value' arrays for prescribed displacements

        Returns
        -------
        K_mod : csr_matrix
        F_mod : np.ndarray
        fixed_dofs : np.ndarray
        """
        logger.info("🔒 Applying boundary conditions to global matrices…")

        # 1. Instantiate the boundary-condition helper
        modifier = ModifyGlobalSystem(
            K_global=K_global,
            F_global=np.asarray(F_global).ravel(),
            job_results_dir=job_results_dir,
            fixed_dofs=fixed_dofs,
            prescribed_displacements=prescribed_displacements,
            local_global_dof_map=local_global_dof_map  # ← pass map once
        )

        # 2. Run the BC pipeline
        try:
            K_mod, F_mod, fixed_dofs, modification_map = modifier.apply_boundary_conditions()
        except Exception:
            logger.error("❌ Boundary-condition step failed – see ModifyGlobalSystem.log")
            raise

        # Cache into self
        self.K_mod = K_mod
        self.F_mod = F_mod
        self.fixed_dofs = fixed_dofs

        # Cache map and results to simulation runner results container
        self.maps.modification_map = modification_map
        self.primary_results_set.global_results.K_mod = K_mod
        self.primary_results_set.global_results.F_mod = F_mod

        logger.info("✅ Boundary conditions successfully applied")
        return K_mod, F_mod, fixed_dofs

    # ----------------------------------------=---------------------------------
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
            base_tol: float = 1e-12
        ):
        """
        Perform static condensation on the boundary-conditioned system using the
        high-performance `CondenseModifiedSystem` helper.

        Parameters
        ----------
        K_mod : scipy.sparse.csr_matrix
            Global stiffness matrix *after* boundary conditions have been applied.
        F_mod : np.ndarray
            Corresponding force vector (flattened to 1-D internally).
        fixed_dofs : np.ndarray
            Array of fixed DOF indices (passed straight through to the helper).
        job_results_dir : str
            Directory where condensation logs/diagnostics should be written.
        base_tol : float, optional
            Baseline tolerance; the helper auto-scales from this value.

        Returns
        -------
        condensed_dofs : np.ndarray
            DOF indices retained in the condensed system (w.r.t. the **original**
            global numbering).
        inactive_dofs : np.ndarray
            Active DOFs that were fully eliminated by the secondary condensation.
        K_cond : scipy.sparse.csr_matrix
            Condensed stiffness matrix.
        F_cond : np.ndarray
            Condensed force vector (1-D).
        """
        logger.info("📉 Condensing modified system…")

        
        # 1.  Instantiate the condensation helper
        
        condenser = CondenseModifiedSystem(
            K_mod        = K_mod,
            F_mod        = np.asarray(F_mod).ravel(),   # ensure 1-D
            fixed_dofs   = fixed_dofs,
            local_global_dof_map=local_global_dof_map,
            job_results_dir = job_results_dir,
            base_tol        = base_tol
        )

        # 2.  Execute the condensation pipeline
        
        try:
            condensed_dofs, inactive_dofs, K_cond, F_cond, condensation_map = (
                condenser.apply_condensation()
            )
        except Exception as exc:                       # helper already logs detail
            logger.error("❌ Condensation failed – see condensation.log for details")
            raise                                       # propagate to caller

        
        # Cache into self
        self.K_cond = K_cond
        self.F_cond = F_cond
        self.inactive_dofs = inactive_dofs
        self.condensed_dofs = condensed_dofs

        # Cache map and results to simulation runner results container
        self.maps.condensation_map = condensation_map
        self.primary_results_set.global_results.K_cond = K_cond
        self.primary_results_set.global_results.F_cond = F_cond

        # 3.  High-level diagnostics (optional, keeps your existing style)

        logger.info(
            "✅ Static condensation complete "
            f"({len(condensed_dofs)} DOFs retained, "
            f"{len(inactive_dofs)} pruned)"
        )
        return condensed_dofs, inactive_dofs, K_cond, F_cond
    
    # -------------------------------------------------------------------------
    # 4) SOLVE CONDENSED SYSTEM
    # -------------------------------------------------------------------------

    def solve_condensed_system(self,
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
                            disable_scaling: bool = False):
        """
        Solve the condensed linear system K_cond · U_cond = F_cond.

        Parameters
        ----------
        K_cond, F_cond : csr_matrix, np.ndarray
            Condensed stiffness matrix and force vector.
        job_results_dir : str
            Folder for logs and CSV output.
        solver_name : {"cg","gmres","bicgstab","direct",...}, optional
            Which registered solver to use (default "cg").
        preconditioner : {"auto","ilu","jacobi",None}, optional
            Preconditioning strategy (default "auto").
        tolerance : float, optional
            Solver tolerance (default 1e-6).
        max_iterations : int, optional
            Maximum solver iterations (default 1000).
        restart : int, optional
            Restart parameter for GMRES (default 20).
        ilu_drop_tol : float, optional
            ILU drop tolerance (default 1e-6).
        ilu_fill_factor : float, optional
            ILU fill factor (default 1.0).
        disable_scaling : bool, optional
            Disable row/column scaling (default False).

        Returns
        -------
        U_cond : np.ndarray
            Displacements in condensed system ordering.
        """
        logger.info("🧮 Solving condensed system…")

        cond_solver = SolveCondensedSystem(
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
            disable_scaling=disable_scaling
        )

        U_cond = cond_solver.solve()
        if U_cond is None:
            raise RuntimeError(
                "Condensed solver failed — see SolveCondensedSystem.log")
        
        # Cache into self
        self.U_cond = U_cond

        # Cache map and results to simulation runner results container
        self.primary_results_set.global_results.U_cond = U_cond

        logger.info("✅ Condensed system successfully solved")
        return U_cond

    # ---------------------------------------------------------------------------------
    # 5) RECONSTRUCT GLOBAL SYSTEM
    # ---------------------------------------------------------------------------------

    def reconstruct_global_system(self,
                                    condensed_dofs: np.ndarray,
                                    U_cond: np.ndarray,
                                    total_dof: int,
                                    job_results_dir: str | Path,
                                    *,
                                    fixed_dofs: np.ndarray | None = None,
                                    inactive_dofs: np.ndarray | None = None,
                                    local_global_dof_map,) -> np.ndarray:
        """
        Reconstruct the full displacement vector from the condensed solution.

        Parameters
        ----------
        condensed_dofs : np.ndarray
            Active DOF indices in original numbering.
        U_cond : np.ndarray
            Displacements solved in condensed ordering.
        total_dof : int
            Total DOFs in the original system.
        job_results_dir : str | Path
            Directory for reconstruction logs and CSV output.
        fixed_dofs : np.ndarray | None, optional
            Fixed DOF indices for validation.
        inactive_dofs : np.ndarray | None, optional
            Inactive DOFs removed during condensation.

        Returns
        -------
        np.ndarray
            Reconstructed global displacement vector (
            ).
        """
        logger.info("🔄 Reconstructing global displacement vector…")

        reconstructor = ReconstructGlobalSystem(
            active_dofs=condensed_dofs,
            U_cond=U_cond,
            total_dofs=total_dof,
            job_results_dir=Path(job_results_dir),
            fixed_dofs=fixed_dofs,
            inactive_dofs=inactive_dofs,
            local_global_dof_map = local_global_dof_map
        )

        try:
            U_global, reconstruction_map = reconstructor.reconstruct()
        except Exception:
            logger.error("❌ Reconstruction failed – see ReconstructGlobalSystem.log")
            raise

        # Cache into self
        self.U_global = U_global

        # Cache map and results to simulation runner results container
        self.primary_results_set.global_results.U_global = U_global
        self.maps.reconstruction_map = reconstruction_map

        logger.info(
            f"✅ Displacement reconstruction complete "
            f"(min={U_global.min():.3e}, max={U_global.max():.3e})"
        )
        return U_global

    # ─────────────────────────────────────────────────────────────────────────────
    # 6) PRIMARY RESULTS PIPELINE
    # ─────────────────────────────────────────────────────────────────────────────
    def compute_primary_results(self) -> None:
        """
        1. Build all global-level first-order outputs (incl. reactions).
        2. Disassemble U_global & R_global back onto each element.
        3. Cache the two result sets in-memory.
        (All CSVs are persisted by the helpers themselves.)
        """
        logger.info("📊 Computing primary results …")

        # Step 1: Validate and standardize element matrices from COO to CSR (in-place)
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
            job_results_dir=self.primary_results_dir
        ).compute()

        # Step 3: Disassemble element-level primary results
        self.R_e, self.R_residual_e, self.U_e = DisassembleGlobalSystem(
            U_global=self.U_global,
            R_global=self.R_global,
            R_residual=self.R_residual,
            F_e = self.element_force_vectors,
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
            R_residual=self.R_residual
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
            elemental_results=elemental_results
        )

        
        # Step 7: Save .csv files (all resolutions: global + elemental)
        SavePrimaryResults(
            primary_results_set=self.primary_results,
            index_map_set=self.maps,
            save_dir=self.primary_results_dir,
        ).run()

        SavePrimaryResultsSummary(
            primary_results_set=self.primary_results,
            save_dir=self.primary_results_dir,
        ).save()

        SaveIndexMaps(
            index_map_set = self.maps,
            save_dir = self.maps_dir                                      #os.path.join(self.maps_dir, "NEW")
        ).run()

        logger.info("✅ Primary results computed, written and saved")

    # -------------------------------------------------------------------------
    # 7) SECONDARY RESULTS PIPELINE
    # -------------------------------------------------------------------------
    def compute_secondary_results(self) -> None:
        """Compute and save secondary results (stresses, strains, etc.)."""
        logger.info("📈 Computing secondary results using cached formulation data...")
        
        from processing_OOP.static.results.compute_secondary.compute_secondary_results import SecondaryResultsOrchestrator
        from processing_OOP.static.results.save_secondary_container import SaveSecondaryResults, SaveSecondaryResultsSummary

        # Compute secondary results using cached Gauss point data
        orchestrator = SecondaryResultsOrchestrator(
            elements=self.elements,
            grid_dictionary=self.grid_dictionary,
            element_dictionary=self.element_dictionary,
            material_dictionary=self.material_dictionary,
            section_dictionary=self.section_dictionary,
            global_displacement=self.U_global,
            formulation_cache=self.formulation_cache,
            job_results_dir=self.secondary_results_dir
        )
        
        self.secondary_results_set = orchestrator.compute_all()
        
        # Save secondary results (all resolutions)
        saver = SaveSecondaryResults(
            secondary_results=self.secondary_results_set,
            save_dir=self.secondary_results_dir,
        )
        saver.save_all()

        SaveSecondaryResultsSummary(
            secondary_results=self.secondary_results_set,
            save_dir=self.secondary_results_dir,
        ).save()

        logger.info("✅ Secondary results computed and saved")

    # -------------------------------------------------------------------------
    # 8) TERTIARY RESULTS PIPELINE
    # -------------------------------------------------------------------------
    def compute_tertiary_results(self) -> None:
        """Compute and save tertiary results (section forces, principal stresses, etc.)."""
        logger.info("📊 Computing tertiary results...")

        from processing_OOP.static.results.compute_tertiary.compute_tertiary_results import TertiaryResultsOrchestrator
        from processing_OOP.static.results.save_tertiary_container import SaveTertiaryResults, SaveTertiaryResultsSummary

        orchestrator = TertiaryResultsOrchestrator(
            secondary_results=self.secondary_results_set,
            formulation_cache=self.formulation_cache,
            job_results_dir=self.tertiary_results_dir,
        )
        self.tertiary_results = orchestrator.compute()

        # Save tertiary results (all resolutions: gaussian + elemental)
        SaveTertiaryResults(
            tertiary_results=self.tertiary_results,
            save_dir=self.results_root,
        ).save_all()

        SaveTertiaryResultsSummary(
            tertiary_results=self.tertiary_results,
            save_dir=self.results_root,
        ).save()

        logger.info("✅ Tertiary results computed and saved")

    # -------------------------------------------------------------------------
    # TOP-LEVEL SIMULATION FLOW
    # -------------------------------------------------------------------------
    def run(self) -> None:
        """Execute the complete linear-static workflow."""
        try:
            # initialise once – writes machine header
            self.monitor = RuntimeMonitorTelemetry(job_results_dir=self.diagnostics_dir)

            # 0) housekeeping & local systems
            self.setup_simulation()                       # create folders
            self.start_time = datetime.datetime.now()

            # -----------------------------------------------------------------
            # 0  Prepare local element data
            # -----------------------------------------------------------------
            with self.monitor.stage("PrepareLocalSystem"):
                self.prepare_local_system(job_results_dir=self.primary_results_dir)

            # -----------------------------------------------------------------
            # 1  Global assembly
            # -----------------------------------------------------------------
            with self.monitor.stage("AssembleGlobalSystem"):
                (self.K_global,
                 self.F_global,
                 self.local_global_dof_map) = self.assemble_global_system(self.primary_results_dir)

            DiagnoseLinearStaticSystem(stage="Global System",
                                       A_current       = self.K_global,  
                                       b_current       = self.F_global,
                                       A_full          = self.K_global,  
                                       b_full          = self.F_global,
                                       fixed_dofs      = [], 
                                       condensed_dofs  = [],
                                       job_results_dir = self.primary_results_dir,)

            # -----------------------------------------------------------------
            # 2  Boundary conditions
            # -----------------------------------------------------------------
            with self.monitor.stage("ModifyGlobalSystem"):
                (self.K_mod,
                 self.F_mod,
                 self.fixed_dofs) = self.modify_global_system(
                    self.K_global,
                    self.F_global,
                    self.local_global_dof_map,
                    self.primary_results_dir,
                    prescribed_displacements=getattr(self, 'prescribed_displacements', None)
                )

            DiagnoseLinearStaticSystem(stage="Modified System",
                                       A_current       = self.K_mod,
                                       b_current       = self.F_mod,
                                       A_full          = self.K_global,
                                       b_full          = self.F_global,
                                       fixed_dofs      = self.fixed_dofs,
                                       condensed_dofs  = [],
                                       job_results_dir = self.primary_results_dir,)

            # -----------------------------------------------------------------
            # 3  Static condensation
            # -----------------------------------------------------------------
            with self.monitor.stage("CondenseModifiedSystem"):
                (self.condensed_dofs,
                 self.inactive_dofs,
                 self.K_cond,
                 self.F_cond) = self.condense_modified_system(
                    self.K_mod, 
                    self.F_mod,
                    self.fixed_dofs,
                    self.local_global_dof_map,
                    self.primary_results_dir,
                    base_tol=self.condensation_config.get("base_tol", 1e-12)
                )

            DiagnoseLinearStaticSystem(stage="Condensed System",
                                       A_current       = self.K_cond,
                                       b_current       = self.F_cond,
                                       A_full          = self.K_global,
                                       b_full          = self.F_global,
                                       fixed_dofs      = self.fixed_dofs,
                                       condensed_dofs  = self.condensed_dofs,
                                       job_results_dir = self.primary_results_dir,)

            # -----------------------------------------------------------------
            # 4  Solve condensed system
            # -----------------------------------------------------------------
            with self.monitor.stage("SolveCondensedSystem"):
                self.U_cond = self.solve_condensed_system(
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
                    disable_scaling=self.solver_config.get("disable_scaling", False)
                )

            # -----------------------------------------------------------------
            # 5  Reconstruct full-length displacement vector
            # -----------------------------------------------------------------
            with self.monitor.stage("ReconstructGlobalSystem"):
                total_dofs = self.total_dof
                self.U_global = self.reconstruct_global_system(self.condensed_dofs,
                                                               self.U_cond,total_dofs,
                                                               self.primary_results_dir,
                                                               fixed_dofs=self.fixed_dofs,
                                                               inactive_dofs=self.inactive_dofs,
                                                               local_global_dof_map=self.local_global_dof_map)

            # -----------------------------------------------------------------
            # 6  Primary results (includes disassembly & CSV export)
            # -----------------------------------------------------------------
            with self.monitor.stage("ComputePrimaryResults"):
                self.compute_primary_results()
            
            # -----------------------------------------------------------------
            # 6.5 Save formulation cache
            # -----------------------------------------------------------------
            with self.monitor.stage("SaveFormulationCache"):
                from processing_OOP.static.results.save_formulation_container import SaveFormulationData
                saver = SaveFormulationData(
                    formulation_cache=self.formulation_cache,
                    save_dir=self.primary_results_dir,
                    save_gauss_data=False  # Can be enabled for debugging
                )
                saver.save_all()

            # -----------------------------------------------------------------
            # 7  Secondary (derived) results
            # -----------------------------------------------------------------
            with self.monitor.stage("ComputeSecondaryResults"):
                self.compute_secondary_results()

            # -----------------------------------------------------------------
            # 8  Tertiary (design/verification) results
            # -----------------------------------------------------------------
            with self.monitor.stage("ComputeTertiaryResults"):
                self.compute_tertiary_results()

            logger.info("🏁 Simulation completed successfully → %s", self.results_root)

        except Exception as exc:
            logger.exception("💥 Simulation failed with critical error")
            raise RuntimeError("Simulation aborted") from exc