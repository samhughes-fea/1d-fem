# processing\static\results\compute_tertiary\integrated_elemental_results.py

"""
Integrated Elemental Results

Computes integrated quantities per element from Gaussian-level results:
- Total strain energy per element
- Integrated section forces per element

These are obtained by integrating field quantities over the element domain
using Gauss quadrature.
"""

import numpy as np
from typing import List, Optional
import logging

logger = logging.getLogger(__name__)


class ComputeIntegratedElementalResults:
    """
    Computes integrated elemental results from Gaussian-level field quantities.
    
    Integrates:
    1. Strain energy density → Total strain energy per element
    2. Section forces → Integrated/average section forces per element
    """
    
    def __init__(
        self,
        secondary_results,
        formulation_cache
    ):
        """
        Parameters
        ----------
        secondary_results
            SecondaryResultSet with Gaussian results
        formulation_cache
            FormulationResultSet with element objects containing Gauss point data
        """
        self.secondary_results = secondary_results
        self.formulation_cache = formulation_cache
        
    def compute_total_strain_energy(self) -> List[float]:
        """
        Compute total strain energy per element by integrating energy density.
        
        Formula: U_e = ∫_Ω w(x) dΩ ≈ ∑_g w_g ⋅ w(x_g) ⋅ |J(x_g)|
        
        Returns
        -------
        List[float]
            Total strain energy per element (units: J)
            Shape: List[element] -> float
        """
        total_energies = []
        
        gaussian_results = self.secondary_results.gaussian_results
        if gaussian_results is None or gaussian_results.internal_energy_density is None:
            logger.warning("No energy density data available for integration")
            return []
        
        # Build mapping from element_id to index in GaussianResults
        elem_id_to_idx = {}
        for idx, elem_obj in enumerate(self.formulation_cache.element_objects):
            elem_id_to_idx[elem_obj.element_id] = idx
        
        for elem_obj in self.formulation_cache.element_objects:
            elem_id = elem_obj.element_id
            gauss_idx = elem_id_to_idx.get(elem_id)
            
            if gauss_idx is None:
                logger.warning(f"Element {elem_id} not found in GaussianResults, skipping")
                total_energies.append(0.0)
                continue
            
            # Get energy density at Gauss points for this element
            energy_densities = gaussian_results.internal_energy_density[gauss_idx]
            
            # Integrate using Gauss quadrature
            total_energy = 0.0
            for gp, energy_density in zip(elem_obj.gauss_data, energy_densities):
                # U_e = ∫ w dΩ ≈ ∑ w_g ⋅ w(x_g) ⋅ |J(x_g)|
                total_energy += gp.weight * float(energy_density) * gp.jacobian
            
            total_energies.append(total_energy)
        
        logger.info(f"✅ Computed total strain energy for {len(total_energies)} elements")
        return total_energies
    
    def compute_integrated_section_forces(
        self,
        section_forces_gauss: List[List[np.ndarray]]
    ) -> List[np.ndarray]:
        """
        Compute integrated/average section forces per element.
        
        For beam elements, this computes the average section force over the
        element length, weighted by the element geometry.
        
        Formula: F_avg = (1/L) ∫_Ω F(x) dΩ ≈ (1/L) ∑_g w_g ⋅ F(x_g) ⋅ |J(x_g)|
        
        Parameters
        ----------
        section_forces_gauss : List[List[np.ndarray]]
            Section forces at Gauss points
            Shape: List[element] -> List[gauss_point] -> np.ndarray(6,)
        
        Returns
        -------
        List[np.ndarray]
            Integrated section forces per element [N, Vy, Vz, T, My, Mz]
            Shape: List[element] -> np.ndarray(6,)
        """
        integrated_forces = []
        
        # Build mapping from element_id to index
        elem_id_to_idx = {}
        for idx, elem_obj in enumerate(self.formulation_cache.element_objects):
            elem_id_to_idx[elem_obj.element_id] = idx
        
        for elem_idx, elem_obj in enumerate(self.formulation_cache.element_objects):
            elem_id = elem_obj.element_id
            
            if elem_idx >= len(section_forces_gauss):
                logger.warning(f"Element {elem_id} index out of range, skipping")
                integrated_forces.append(np.zeros(6))
                continue
            
            # Get section forces at Gauss points for this element
            elem_section_forces = section_forces_gauss[elem_idx]
            
            if len(elem_section_forces) != len(elem_obj.gauss_data):
                logger.warning(
                    f"Element {elem_id}: Mismatch between Gauss points "
                    f"({len(elem_obj.gauss_data)}) and section forces ({len(elem_section_forces)})"
                )
                integrated_forces.append(np.zeros(6))
                continue
            
            # Integrate section forces using Gauss quadrature
            # For beam elements, compute weighted average over element length
            integrated_force = np.zeros(6)
            total_weight = 0.0
            
            for gp, section_force in zip(elem_obj.gauss_data, elem_section_forces):
                # Weighted contribution: w_g ⋅ F(x_g) ⋅ |J(x_g)|
                weight_contrib = gp.weight * gp.jacobian
                integrated_force += weight_contrib * np.array(section_force)
                total_weight += weight_contrib
            
            # Normalize by total weight to get average
            if total_weight > 1e-10:
                integrated_force /= total_weight
            else:
                logger.warning(f"Element {elem_id}: Zero total weight in integration")
            
            integrated_forces.append(integrated_force)
        
        logger.info(f"✅ Computed integrated section forces for {len(integrated_forces)} elements")
        return integrated_forces

