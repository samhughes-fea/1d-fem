# processing_OOP\static\results\compute_tertiary\section_force.py

import numpy as np
from typing import List, Optional
import logging

logger = logging.getLogger(__name__)


class ComputeSectionForce:
    """
    Computes internal section force resultants [N, Vy, Vz, T, My, Mz] at Gauss points.

    Section forces represent the **physical stress-resultant vectors** at specific
    points along the element, describing the internal state due to deformation.
    They are essential for **stress recovery**, **failure checks**, and **engineering
    interpretation**.

    Unlike nodal forces (K_e @ U_e), which express nodal reactions, section forces
    reflect local **internal equilibrium** across the cross-section:

        ε = B @ U_e                → strain at Gauss point
        σ = D @ ε                 → stress from material law
        [N, Vy, Vz, T, My, Mz]    → integrated from σ over cross-section

    For 1D beam theory, these are directly obtained from the stress vector.

    Parameters
    ----------
    stress_gauss : List[List[np.ndarray]]
        Stress tensors at Gauss points for all elements
        Shape: List[element] -> List[gauss_point] -> np.ndarray(6,)
    
    Attributes
    ----------
    stress_gauss : List[List[np.ndarray]]
        Input stress data
    """

    def __init__(self, stress_gauss: List[List[np.ndarray]]):
        self.stress_gauss = stress_gauss

    def compute(self) -> List[List[np.ndarray]]:
        """
        Compute section force resultants at each Gauss point for all elements.

        Returns
        -------
        List[List[np.ndarray]]
            Section forces [N, Vy, Vz, T, My, Mz] per Gauss point per element
            Shape: List[element] -> List[gauss_point] -> np.ndarray(shape=(6,))
        """
        section_forces = []

        for elem_idx, elem_stresses in enumerate(self.stress_gauss):
            elem_section_forces = []
            
            for gp_idx, stress in enumerate(elem_stresses):
                # For beam elements, stress vector is already in section force form
                # stress = [σ_axial, τ_shear_y, τ_shear_z, τ_torsion, σ_bend_y, σ_bend_z]
                # maps directly to section forces:
                # section_force = [N, Vy, Vz, T, My, Mz]
                
                section_force = self._stress_to_section_force(stress)
                elem_section_forces.append(section_force)

            section_forces.append(elem_section_forces)

        logger.info(f"✅ Computed section forces for {len(section_forces)} elements")
        return section_forces

    def _stress_to_section_force(self, stress: np.ndarray) -> np.ndarray:
        """
        Convert stress tensor to section force resultants.

        For beam elements, the stress vector components typically represent
        integrated quantities:
        - stress[0] → N  (axial force)
        - stress[1] → Vy (shear force y)
        - stress[2] → Vz (shear force z)
        - stress[3] → T  (torque)
        - stress[4] → My (bending moment y)
        - stress[5] → Mz (bending moment z)

        Parameters
        ----------
        stress : np.ndarray
            Stress vector at Gauss point (shape: (6,))

        Returns
        -------
        np.ndarray
            Section force vector [N, Vy, Vz, T, My, Mz] (shape: (6,))
        """
        # For typical beam formulations, this is a direct mapping
        # For more complex elements, integration over cross-section would be needed
        return stress.copy()


class ComputeNodalSectionForce:
    """
    Projects Gaussian section forces to nodal values using shape function
    extrapolation or averaging schemes.

    This is useful for visualization and post-processing at nodes rather
    than at integration points.
    """

    def __init__(
        self,
        section_forces_gauss: List[List[np.ndarray]],
        connectivity: np.ndarray,
        n_nodes: int,
        method: str = "average"
    ):
        """
        Parameters
        ----------
        section_forces_gauss : List[List[np.ndarray]]
            Section forces at Gauss points
        connectivity : np.ndarray
            Element connectivity array (n_elements x 2)
        n_nodes : int
            Total number of nodes in mesh
        method : str
            Projection method ('average' or 'extrapolate')
        """
        self.section_forces_gauss = section_forces_gauss
        self.connectivity = connectivity
        self.n_nodes = n_nodes
        self.method = method

    def compute(self) -> np.ndarray:
        """
        Project Gauss point section forces to nodes.

        Returns
        -------
        np.ndarray
            Nodal section forces (n_nodes x 6)
        """
        nodal_forces = np.zeros((self.n_nodes, 6))
        nodal_counts = np.zeros(self.n_nodes)

        if self.method == "average":
            return self._average_projection(nodal_forces, nodal_counts)
        elif self.method == "extrapolate":
            return self._extrapolate_projection()
        else:
            raise ValueError(f"Unknown projection method: {self.method}")

    def _average_projection(
        self,
        nodal_forces: np.ndarray,
        nodal_counts: np.ndarray
    ) -> np.ndarray:
        """Simple averaging of adjacent element contributions."""
        for elem_idx, elem_section_forces in enumerate(self.section_forces_gauss):
            # Get element nodes
            node_ids = self.connectivity[elem_idx]
            
            # Average all Gauss point values for this element
            elem_avg = np.mean(elem_section_forces, axis=0)
            
            # Distribute to both nodes
            for node_id in node_ids:
                nodal_forces[node_id] += elem_avg
                nodal_counts[node_id] += 1

        # Normalize by contribution count
        mask = nodal_counts > 0
        nodal_forces[mask] /= nodal_counts[mask, np.newaxis]

        return nodal_forces

    def _extrapolate_projection(self) -> np.ndarray:
        """
        Extrapolate from Gauss points to nodes using shape functions.
        This is more accurate but requires knowing the shape functions.
        """
        # Placeholder for extrapolation method
        # Would require access to element shape functions
        raise NotImplementedError("Extrapolation method not yet implemented")

