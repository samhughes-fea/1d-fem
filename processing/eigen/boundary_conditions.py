# processing/eigen/boundary_conditions.py
"""
Apply boundary conditions to global K and M for §2 eigen problems (penalty method).

Copied and extended from static modification logic; no imports from ``processing.static``.
"""

import logging
import numpy as np
from scipy.sparse import csr_matrix
from typing import Optional, Sequence, Tuple, Union

_logger = logging.getLogger(__name__)

PENALTY = 1e36


def apply_boundary_conditions(
    K_global: csr_matrix,
    M_global: csr_matrix,
    fixed_dofs: Optional[Union[Sequence[int], np.ndarray]] = None,
    prescribed_displacements: Optional[dict] = None,
) -> Tuple[csr_matrix, csr_matrix, np.ndarray]:
    """
    Apply BCs to both K and M (penalty method): fixed DOFs get zero row/col,
    diagonal K = penalty, diagonal M = 1.0.

    Parameters
    ----------
    K_global, M_global : csr_matrix
        Global stiffness and mass matrices.
    fixed_dofs : array-like, optional
        DOF indices to fix. Default: first 6 (rigid body).
    prescribed_displacements : dict, optional
        Keys 'global_dof' and 'value'. Zero values are treated as fixed.

    Returns
    -------
    K_mod, M_mod : csr_matrix
    bc_dofs : np.ndarray
        All constrained DOF indices (fixed + prescribed zero).
    """
    n = K_global.shape[0]
    if M_global.shape != (n, n) or K_global.shape != (n, n):
        raise ValueError("K and M must be square and same size")

    if prescribed_displacements is not None:
        prescribed_dofs = np.asarray(prescribed_displacements["global_dof"], dtype=np.int32)
        prescribed_values = np.asarray(prescribed_displacements["value"], dtype=np.float64)
        zero_mask = np.abs(prescribed_values) < 1e-12
        prescribed_zero = prescribed_dofs[zero_mask]
        if fixed_dofs is not None:
            bc_dofs = np.unique(np.concatenate([np.asarray(fixed_dofs, dtype=np.int32), prescribed_zero]))
        else:
            bc_dofs = np.unique(prescribed_zero)
    else:
        bc_dofs = np.asarray(fixed_dofs if fixed_dofs is not None else range(6), dtype=np.int32)

    bc_dofs = np.unique(bc_dofs)
    if bc_dofs.size == 0:
        _logger.debug("No boundary conditions applied")
        return K_global.copy(), M_global.copy(), bc_dofs

    K_work = K_global.astype(np.float64).tolil()
    M_work = M_global.astype(np.float64).tolil()

    for dof in bc_dofs:
        K_work[dof, :] = 0.0
        K_work[:, dof] = 0.0
        K_work[dof, dof] = PENALTY
        M_work[dof, :] = 0.0
        M_work[:, dof] = 0.0
        M_work[dof, dof] = 1.0

    K_mod = K_work.tocsr()
    M_mod = M_work.tocsr()
    _logger.info("Eigen BCs applied: %s DOFs constrained", bc_dofs.size)
    return K_mod, M_mod, bc_dofs
