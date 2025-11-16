# processing_OOP\static\results\compute_secondary\compute_secondary_results.py

"""
Secondary Results Orchestrator

Computes derived field quantities (strain, stress, energy density) using
cached Gauss point formulation data from element computation phase.

This eliminates the need to recompute shape functions and B matrices,
improving performance and ensuring exact consistency with stiffness computation.
"""

import numpy as np
from pathlib import Path
from typing import Optional
import logging

from processing_OOP.static.results.containers import (
    SecondaryResultSet,
    GaussianResults,
    NodalResults,
    FormulationResultSet
)


class SecondaryResultsOrchestrator:
    """
    Orchestrates computation of secondary results using cached formulation data.
    
    Leverages Gauss-level caching from FormulationResultSet to compute:
    - Strains at Gauss points
    - Stresses at Gauss points
    - Strain energy density at Gauss points
    - Nodal interpolations (optional)
    """
    
    def __init__(
        self,
        *,
        elements,
        grid_dictionary,
        element_dictionary,
        material_dictionary,
        section_dictionary,
        global_displacement: np.ndarray,
        formulation_cache: FormulationResultSet,
        job_results_dir: Optional[Path] = None,
    ):
        """
        Parameters
        ----------
        elements : list
            List of element objects
        grid_dictionary : dict
            Node/grid data
        element_dictionary : dict
            Element connectivity data
        material_dictionary : dict
            Material properties
        section_dictionary : dict
            Cross-section properties
        global_displacement : np.ndarray
            Global displacement vector
        formulation_cache : FormulationResultSet
            Cached element formulation data with Gauss point information
        job_results_dir : Path, optional
            Directory for logging
        """
        self.elements = elements
        self.grid_dictionary = grid_dictionary
        self.element_dictionary = element_dictionary
        self.material_dictionary = material_dictionary
        self.section_dictionary = section_dictionary
        self.U_global = global_displacement.reshape(-1)
        self.formulation_cache = formulation_cache
        self.job_results_dir = Path(job_results_dir) if job_results_dir else None
        self.logger = self._init_logging()
    
    def _init_logging(self):
        """Initialize logger with optional file output."""
        logger = logging.getLogger(f"SecondaryResultsOrchestrator.{id(self)}")
        logger.handlers.clear()
        logger.setLevel(logging.DEBUG)
        logger.propagate = False
        
        if self.job_results_dir:
            logs_dir = self.job_results_dir.parent / "logs"
            logs_dir.mkdir(parents=True, exist_ok=True)
            log_path = logs_dir / "SecondaryResults.log"
            
            try:
                file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
                file_handler.setFormatter(logging.Formatter(
                    "%(asctime)s [%(levelname)s] %(message)s"
                ))
                logger.addHandler(file_handler)
            except Exception as e:
                print(f"⚠️ Failed to create log file for SecondaryResults: {e}")
        
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        logger.addHandler(stream_handler)
        
        return logger
    
    def compute_all(self) -> SecondaryResultSet:
        """
        Compute all secondary results using cached Gauss data.
        
        Returns
        -------
        SecondaryResultSet
            Container with Gaussian and nodal results
        """
        self.logger.info("🔬 Computing secondary results from cached formulation data...")
        
        gaussian_results = self.compute_gaussian_results()
        nodal_results = self.compute_nodal_results(gaussian_results)
        
        return SecondaryResultSet(
            gaussian_results=gaussian_results,
            nodal_results=nodal_results
        )
    
    def compute_gaussian_results(self) -> GaussianResults:
        """
        Compute strain, stress, and energy density at Gauss points using cached B and D matrices.
        
        Returns
        -------
        GaussianResults
            Field quantities at Gauss integration points
        """
        self.logger.info("  Computing Gaussian field quantities...")
        
        strains_all = []
        stresses_all = []
        energy_all = []
        
        for elem_obj in self.formulation_cache.element_objects:
            elem_id = elem_obj.element_id
            element = self.elements[elem_id]
            
            # Extract element displacement from global
            dof_indices = element.assemble_global_dof_indices()
            U_e = self.U_global[dof_indices]
            
            # Compute at each Gauss point using cached data
            elem_strains = []
            elem_stresses = []
            elem_energy = []
            
            for gp in elem_obj.gauss_data:
                # Use cached B and D matrices - no shape function recomputation!
                strain = gp.B_matrix @ U_e
                stress = gp.D_matrix @ strain
                energy_density = 0.5 * (strain.T @ stress)
                
                elem_strains.append(strain)
                elem_stresses.append(stress)
                elem_energy.append(float(energy_density))
            
            strains_all.append(elem_strains)
            stresses_all.append(elem_stresses)
            energy_all.append(elem_energy)
        
        self.logger.info(f"  ✅ Computed results for {len(strains_all)} elements")
        
        return GaussianResults(
            strain=strains_all,
            stress=stresses_all,
            internal_energy_density=energy_all
        )
    
    def compute_nodal_results(self, gaussian_results: GaussianResults) -> NodalResults:
        """
        Interpolate Gaussian results to nodes (if needed).
        
        Parameters
        ----------
        gaussian_results : GaussianResults
            Field quantities at Gauss points
        
        Returns
        -------
        NodalResults
            Field quantities interpolated to nodes
        """
        self.logger.info("  Computing nodal field quantities...")
        
        # For now, return empty nodal results - can implement extrapolation later
        # This would involve extrapolating from Gauss points to element nodes
        
        return NodalResults(
            strain=None,
            stress=None,
            strain_energy_density=None
        )
