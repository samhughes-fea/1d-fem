# processing\static\results\compute_secondary\nodal_force.py

"""
OPTIONAL / UNWIRED: Nodal internal forces (F_int = K_e @ U_e) are not in the current pipeline.

Primary already provides K_e, U_e, R_global, and R_residual. Nodal internal forces
F_int = K_e @ U_e could be added to the PRIMARY pipeline (they use only primary data)
if needed for equilibrium checks or reporting. This module is kept for reference
but is not called by any orchestrator. To use: wire into primary results or a
dedicated step after disassembly.
"""

import numpy as np
from typing import List, Optional
import logging

logger = logging.getLogger(__name__)


class ComputeNodalForce:
    """
    OPTIONAL / UNWIRED. Computes the internal (nodal) force vector F_int for finite elements.

    Not called by any orchestrator. If needed, wire into primary (uses K_e, U_e only).

    This vector is the **element nodal internal force vector** that expresses
    the internal forces the element "pushes back with" at its nodes in response
    to deformation. It is computed via:

        F_int = K_e @ U_e

    where:
        - K_e is the element stiffness matrix
        - U_e is the element nodal displacement vector

    These internal nodal forces are **virtual reaction forces** that, when
    assembled globally, enforce global equilibrium:

        ∑_e (K_e @ U_e) = F_ext

    While the vector has a nodal structure, its values represent
    **stress-resultants** derived from the strain-displacement matrix and
    constitutive behavior. They capture how the material and geometry resist
    external loads in terms of nodal degrees of freedom.

    This is computed at the **nodal resolution** (as opposed to Gaussian
    resolution for section forces).

    Parameters
    ----------
    K_e_list : List[np.ndarray]
        Element stiffness matrices for all elements
    U_e_list : List[np.ndarray]
        Element displacement vectors for all elements
    """

    def __init__(
        self,
        K_e_list: List[np.ndarray],
        U_e_list: List[np.ndarray]
    ):
        self.K_e_list = K_e_list
        self.U_e_list = U_e_list

    def compute(self) -> List[np.ndarray]:
        """
        Compute internal nodal force vectors for all elements.

        Returns
        -------
        List[np.ndarray]
            Internal force vectors F_int = K_e @ U_e for each element
            Shape: List[element] -> np.ndarray(n_dof_per_element,)
        """
        if len(self.K_e_list) != len(self.U_e_list):
            raise ValueError(
                f"Mismatch: {len(self.K_e_list)} stiffness matrices "
                f"but {len(self.U_e_list)} displacement vectors"
            )

        internal_forces = []

        for elem_idx, (K_e, U_e) in enumerate(zip(self.K_e_list, self.U_e_list)):
            # Convert sparse matrix to dense if needed
            if hasattr(K_e, 'toarray'):
                K_e = K_e.toarray()

            # Compute F_int = K_e @ U_e
            F_int = K_e @ U_e

            internal_forces.append(F_int)

        logger.info(f"✅ Computed nodal internal forces for {len(internal_forces)} elements")
        return internal_forces


class ComputeGlobalNodalForce:
    """
    Assembles element nodal forces into a global nodal force vector.

    This performs the assembly operation:
        F_global_internal[dof] = ∑_e F_int_e[local_dof]

    for all elements connected to each DOF.
    """

    def __init__(
        self,
        nodal_forces_element: List[np.ndarray],
        assembly_map: List,
        n_global_dofs: int
    ):
        """
        Parameters
        ----------
        nodal_forces_element : List[np.ndarray]
            Element nodal force vectors
        assembly_map : List
            DOF assembly mapping (element -> global)
        n_global_dofs : int
            Total number of global DOFs
        """
        self.nodal_forces_element = nodal_forces_element
        self.assembly_map = assembly_map
        self.n_global_dofs = n_global_dofs

    def assemble(self) -> np.ndarray:
        """
        Assemble element nodal forces to global vector.

        Returns
        -------
        np.ndarray
            Global nodal internal force vector (n_global_dofs,)
        """
        F_global_internal = np.zeros(self.n_global_dofs)

        for elem_idx, (F_e, map_entry) in enumerate(
            zip(self.nodal_forces_element, self.assembly_map)
        ):
            # Get global DOF indices for this element
            global_dofs = map_entry.global_dofs
            
            # Add element contribution to global vector
            F_global_internal[global_dofs] += F_e

        logger.info(f"✅ Assembled global nodal force vector (shape: {F_global_internal.shape})")
        return F_global_internal


class ComputeNodalReactionBalance:
    """
    Computes the nodal reaction balance for equilibrium verification.

    At each node, equilibrium requires:
        F_external + F_internal + F_reaction = 0

    This class computes and checks this balance to verify solution accuracy.
    """

    def __init__(
        self,
        F_external: np.ndarray,
        F_internal: np.ndarray,
        F_reaction: np.ndarray
    ):
        """
        Parameters
        ----------
        F_external : np.ndarray
            Global external force vector
        F_internal : np.ndarray
            Global internal force vector (assembled from elements)
        F_reaction : np.ndarray
            Global reaction force vector (at supports)
        """
        self.F_external = F_external
        self.F_internal = F_internal
        self.F_reaction = F_reaction

    def compute_residual(self) -> np.ndarray:
        """
        Compute nodal force residual.

        The residual represents the out-of-balance forces:
            Residual = F_external - F_internal - F_reaction

        For a converged solution, this should be near machine precision.

        Returns
        -------
        np.ndarray
            Nodal force residual vector
        """
        residual = self.F_external - self.F_internal - self.F_reaction
        
        max_residual = np.max(np.abs(residual))
        logger.info(f"  Max nodal force residual: {max_residual:.6e}")
        
        return residual

    def check_equilibrium(self, tolerance: float = 1e-6) -> bool:
        """
        Check if nodal forces satisfy equilibrium within tolerance.

        Parameters
        ----------
        tolerance : float
            Acceptable residual magnitude

        Returns
        -------
        bool
            True if equilibrium is satisfied
        """
        residual = self.compute_residual()
        max_residual = np.max(np.abs(residual))
        
        is_balanced = max_residual < tolerance
        
        if is_balanced:
            logger.info(f"✅ Nodal equilibrium satisfied (residual < {tolerance:.1e})")
        else:
            logger.warning(
                f"⚠️  Nodal equilibrium NOT satisfied! "
                f"Max residual: {max_residual:.6e} > {tolerance:.1e}"
            )
        
        return is_balanced

