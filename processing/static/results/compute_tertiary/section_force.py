# processing\static\results\compute_tertiary\section_force.py

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

    Note on Euler-Bernoulli
    -----------------------
    Euler-Bernoulli elements have no shear strain (γ_xy = γ_xz = 0), so
    σ = D @ ε yields V_y = V_z = 0. The section force outputs Vy, Vz will
    therefore be zero for EB. Shear force in EB is given by equilibrium
    (V = dM/dx), not by this constitutive output. Shear-deformable elements
    (Timoshenko, Levinson) produce non-zero Vy, Vz from D @ ε.

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
                stress = np.asarray(stress)
                # Section force [N, Vy, Vz, T, My, Mz] is only defined for 6-component beam stress.
                # Bar (2) and truss (3) components are not reordered; output zeros for those elements.
                if stress.size != 6:
                    section_force = np.zeros(6, dtype=stress.dtype)
                else:
                    # Formulation (B, D) outputs stress conjugate to ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]:
                    # stress = [N, M_y, M_z, V_y, V_z, T]. Reorder to [N, Vy, Vz, T, My, Mz].
                    section_force = self._stress_to_section_force(stress)
                elem_section_forces.append(section_force)

            section_forces.append(elem_section_forces)

        logger.info(f"✅ Computed section forces for {len(section_forces)} elements")
        return section_forces

    def _stress_to_section_force(self, stress: np.ndarray) -> np.ndarray:
        """
        Convert formulation stress resultants to section force convention.

        The formulation (B, D) uses strain ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x], so
        σ = D @ ε yields stress in order [N, M_y, M_z, V_y, V_z, T]:
        - stress[0] = N   (axial)
        - stress[1] = M_y (bending about y → x-z plane)
        - stress[2] = M_z (bending about z → x-y plane; activated by u_y)
        - stress[3] = V_y (shear in y; activated by u_y)
        - stress[4] = V_z (shear in z)
        - stress[5] = T  (torsion)

        Section force / CSV / plot convention is [N, Vy, Vz, T, My, Mz].

        Parameters
        ----------
        stress : np.ndarray
            Formulation stress at Gauss point (shape: (6,)) in order [N, M_y, M_z, V_y, V_z, T]

        Returns
        -------
        np.ndarray
            Section force vector [N, Vy, Vz, T, My, Mz] (shape: (6,))
        """
        return np.array([
            stress[0],   # N
            stress[3],   # Vy
            stress[4],   # Vz
            stress[5],   # T
            stress[1],   # My
            stress[2],   # Mz
        ], dtype=stress.dtype)


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

