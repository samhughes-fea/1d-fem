# processing_OOP\static\results\compute_secondary\nodal_result_projector.py

"""
Nodal Result Projector

Projects Gaussian-level field quantities (strain, stress, energy density) 
to nodal resolution using shape function extrapolation.

For each element:
1. Evaluates shape functions at node locations (ξ = -1, +1 for 2-node elements)
2. Extrapolates Gauss point values to element nodes
3. For shared nodes, accumulates and averages contributions from all connected elements
"""

import numpy as np
from typing import List, Dict
import logging

from processing_OOP.static.results.containers import (
    GaussianResults,
    NodalResults
)


class NodalResultProjector:
    """
    Projects Gaussian-level field quantities to nodal resolution.
    
    Uses shape function extrapolation to project strain, stress, and energy
    density from Gauss integration points to element nodes, then averages
    contributions at shared nodes.
    """
    
    def __init__(
        self,
        elements: List,
        element_dictionary: Dict,
        grid_dictionary: Dict,
        gaussian_results: GaussianResults,
        formulation_cache
    ):
        """
        Parameters
        ----------
        elements : List
            List of element objects
        element_dictionary : Dict
            Element connectivity and metadata
        grid_dictionary : Dict
            Node/grid data
        gaussian_results : GaussianResults
            Field quantities at Gauss points
        formulation_cache
            Cached formulation data with element objects
        """
        self.elements = elements
        self.element_dictionary = element_dictionary
        self.grid_dictionary = grid_dictionary
        self.gaussian_results = gaussian_results
        self.formulation_cache = formulation_cache
        
        # Build mapping from element_id to index in GaussianResults
        # This ensures we can reliably access results by element_id rather than relying on order
        self.elem_id_to_gauss_idx = {}
        for idx, elem_obj in enumerate(self.formulation_cache.element_objects):
            self.elem_id_to_gauss_idx[elem_obj.element_id] = idx
        
        # Get number of nodes from coordinates array
        if "coordinates" in grid_dictionary:
            self.n_nodes = len(grid_dictionary["coordinates"])
        elif "grid_dictionary" in grid_dictionary and "coordinates" in grid_dictionary["grid_dictionary"]:
            self.n_nodes = len(grid_dictionary["grid_dictionary"]["coordinates"])
        else:
            raise ValueError("Cannot determine number of nodes from grid_dictionary")
        
        self.logger = logging.getLogger(f"{__name__}.{id(self)}")
        
    def project(self) -> NodalResults:
        """
        Project all Gaussian results to nodal resolution.
        
        Returns
        -------
        NodalResults
            Field quantities interpolated to nodes
        """
        self.logger.info("  Projecting Gaussian results to nodes...")
        
        # Initialize nodal arrays
        nodal_strain = np.zeros((self.n_nodes, 6))
        nodal_stress = np.zeros((self.n_nodes, 6))
        nodal_energy = np.zeros(self.n_nodes)
        nodal_weight = np.zeros(self.n_nodes)  # For averaging shared nodes
        
        # Process each element
        # Use element_id mapping to reliably access GaussianResults
        for elem_obj in self.formulation_cache.element_objects:
            elem_id = elem_obj.element_id
            
            # Get the index in GaussianResults for this element_id
            if elem_id not in self.elem_id_to_gauss_idx:
                self.logger.warning(f"Element {elem_id} not found in GaussianResults mapping, skipping")
                continue
            
            gauss_idx = self.elem_id_to_gauss_idx[elem_id]
            element = self.elements[elem_id]
            
            # Get element node IDs
            node_ids = self._get_element_node_ids(elem_id)
            n_nodes_elem = len(node_ids)
            
            # Natural coordinates of element nodes (typically -1, +1 for 2-node elements)
            xi_nodes = self._get_node_natural_coordinates(n_nodes_elem)
            
            # Get Gauss point data for this element using the mapped index
            elem_strains = self.gaussian_results.strain[gauss_idx]
            elem_stresses = self.gaussian_results.stress[gauss_idx]
            elem_energies = self.gaussian_results.internal_energy_density[gauss_idx]
            
            # Get Gauss point natural coordinates
            xi_gauss = [gp.xi for gp in elem_obj.gauss_data]
            
            # Project each field quantity to nodes
            for node_idx, node_id in enumerate(node_ids):
                xi_node = xi_nodes[node_idx]
                
                # Extrapolate from Gauss points to this node
                strain_node = self._extrapolate_to_node(
                    xi_node, xi_gauss, elem_strains
                )
                stress_node = self._extrapolate_to_node(
                    xi_node, xi_gauss, elem_stresses
                )
                energy_node = self._extrapolate_scalar_to_node(
                    xi_node, xi_gauss, elem_energies
                )
                
                # Accumulate contributions (will average later for shared nodes)
                nodal_strain[node_id] += strain_node
                nodal_stress[node_id] += stress_node
                nodal_energy[node_id] += energy_node
                nodal_weight[node_id] += 1.0
        
        # Normalize by weight to average shared nodes
        nonzero_mask = nodal_weight > 0
        nodal_strain[nonzero_mask] /= nodal_weight[nonzero_mask, np.newaxis]
        nodal_stress[nonzero_mask] /= nodal_weight[nonzero_mask, np.newaxis]
        nodal_energy[nonzero_mask] /= nodal_weight[nonzero_mask]
        
        self.logger.info(f"  ✅ Projected results to {self.n_nodes} nodes")
        
        return NodalResults(
            strain=nodal_strain,
            stress=nodal_stress,
            strain_energy_density=nodal_energy
        )
    
    def _get_element_node_ids(self, elem_id: int) -> List[int]:
        """Get node IDs for an element."""
        # Find element index in dictionary
        elem_idx = np.where(self.element_dictionary["ids"] == elem_id)[0][0]
        node_ids = self.element_dictionary["connectivity"][elem_idx]
        return node_ids.tolist()
    
    def _get_node_natural_coordinates(self, n_nodes: int) -> np.ndarray:
        """
        Get natural coordinates of element nodes.
        
        For 2-node elements: [-1, +1]
        For 3-node elements: [-1, 0, +1]
        etc.
        """
        if n_nodes == 2:
            return np.array([-1.0, 1.0])
        elif n_nodes == 3:
            return np.array([-1.0, 0.0, 1.0])
        else:
            # General case: equally spaced in [-1, 1]
            return np.linspace(-1.0, 1.0, n_nodes)
    
    def _extrapolate_to_node(
        self,
        xi_node: float,
        xi_gauss: List[float],
        gauss_values: List[np.ndarray]
    ) -> np.ndarray:
        """
        Extrapolate vector field from Gauss points to a node.
        
        Uses shape function extrapolation if available, otherwise
        uses polynomial interpolation based on number of Gauss points.
        
        Parameters
        ----------
        xi_node : float
            Natural coordinate of target node
        xi_gauss : List[float]
            Natural coordinates of Gauss points
        gauss_values : List[np.ndarray]
            Field values at Gauss points (each is a vector)
        
        Returns
        -------
        np.ndarray
            Extrapolated field value at node
        """
        n_gauss = len(xi_gauss)
        n_components = gauss_values[0].shape[0]
        
        if n_gauss == 1:
            # Single Gauss point: constant extrapolation
            return gauss_values[0].copy()
        elif n_gauss == 2:
            # Two Gauss points: linear interpolation/extrapolation
            # Using Lagrange interpolation: N1(ξ) = (ξ - ξ2)/(ξ1 - ξ2), N2(ξ) = (ξ - ξ1)/(ξ2 - ξ1)
            xi1, xi2 = xi_gauss[0], xi_gauss[1]
            N1 = (xi_node - xi2) / (xi1 - xi2) if abs(xi1 - xi2) > 1e-10 else 0.0
            N2 = (xi_node - xi1) / (xi2 - xi1) if abs(xi2 - xi1) > 1e-10 else 0.0
            
            result = N1 * gauss_values[0] + N2 * gauss_values[1]
            return result
        else:
            # Multiple Gauss points: use Lagrange polynomial interpolation
            result = np.zeros(n_components)
            for i, xi_g in enumerate(xi_gauss):
                # Lagrange basis function L_i(ξ)
                Li = 1.0
                for j, xi_gj in enumerate(xi_gauss):
                    if i != j:
                        if abs(xi_g - xi_gj) > 1e-10:
                            Li *= (xi_node - xi_gj) / (xi_g - xi_gj)
                        else:
                            Li = 0.0
                            break
                result += Li * gauss_values[i]
            return result
    
    def _extrapolate_scalar_to_node(
        self,
        xi_node: float,
        xi_gauss: List[float],
        gauss_values: List[float]
    ) -> float:
        """
        Extrapolate scalar field from Gauss points to a node.
        
        Parameters
        ----------
        xi_node : float
            Natural coordinate of target node
        xi_gauss : List[float]
            Natural coordinates of Gauss points
        gauss_values : List[float]
            Scalar values at Gauss points
        
        Returns
        -------
        float
            Extrapolated scalar value at node
        """
        n_gauss = len(xi_gauss)
        
        if n_gauss == 1:
            return float(gauss_values[0])
        elif n_gauss == 2:
            # Linear interpolation
            xi1, xi2 = xi_gauss[0], xi_gauss[1]
            N1 = (xi_node - xi2) / (xi1 - xi2) if abs(xi1 - xi2) > 1e-10 else 0.0
            N2 = (xi_node - xi1) / (xi2 - xi1) if abs(xi2 - xi1) > 1e-10 else 0.0
            return float(N1 * gauss_values[0] + N2 * gauss_values[1])
        else:
            # Lagrange polynomial interpolation
            result = 0.0
            for i, xi_g in enumerate(xi_gauss):
                Li = 1.0
                for j, xi_gj in enumerate(xi_gauss):
                    if i != j:
                        if abs(xi_g - xi_gj) > 1e-10:
                            Li *= (xi_node - xi_gj) / (xi_g - xi_gj)
                        else:
                            Li = 0.0
                            break
                result += Li * gauss_values[i]
            return float(result)
