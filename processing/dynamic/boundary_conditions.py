# processing/dynamic/boundary_conditions.py
"""
Apply boundary conditions to K, M, and optionally C for dynamic analysis.
Copied and extended from modal boundary_conditions; no imports from processing.static.
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
    C_global: Optional[csr_matrix] = None,
    fixed_dofs: Optional[Union[Sequence[int], np.ndarray]] = None,
    prescribed_displacements: Optional[dict] = None,
) -> Tuple[csr_matrix, csr_matrix, Optional[csr_matrix], np.ndarray]:
    """
    Apply BCs to K, M, and optionally C (penalty method).

    Returns
    -------
    K_mod, M_mod : csr_matrix
    C_mod : csr_matrix or None
    bc_dofs : np.ndarray
    """
    n = K_global.shape[0]
    if M_global.shape != (n, n) or K_global.shape != (n, n):
        raise ValueError("K and M must be square and same size")
    if C_global is not None and C_global.shape != (n, n):
        raise ValueError("C must be same size as K and M")

    if prescribed_displacements is not None:
        prescribed_dofs = np.asarray(prescribed_displacements["global_dof"], dtype=np.int32)
        prescribed_values = np.asarray(prescribed_displacements["value"], dtype=np.float64)
        zero_mask = np.abs(prescribed_values) < 1e-12
        prescribed_zero = prescribed_dofs[zero_mask]
        bc_dofs = np.unique(
            np.concatenate([np.asarray(fixed_dofs or [], dtype=np.int32), prescribed_zero])
        )
    else:
        bc_dofs = np.asarray(fixed_dofs if fixed_dofs is not None else range(6), dtype=np.int32)

    bc_dofs = np.unique(bc_dofs)
    if bc_dofs.size == 0:
        _logger.debug("No boundary conditions applied")
        return K_global.copy(), M_global.copy(), C_global.copy() if C_global is not None else None, bc_dofs

    K_work = K_global.astype(np.float64).tolil()
    M_work = M_global.astype(np.float64).tolil()
    C_work = C_global.astype(np.float64).tolil() if C_global is not None else None

    for dof in bc_dofs:
        K_work[dof, :] = 0.0
        K_work[:, dof] = 0.0
        K_work[dof, dof] = PENALTY
        M_work[dof, :] = 0.0
        M_work[:, dof] = 0.0
        M_work[dof, dof] = 1.0
        if C_work is not None:
            C_work[dof, :] = 0.0
            C_work[:, dof] = 0.0
            C_work[dof, dof] = 1.0

    K_mod = K_work.tocsr()
    M_mod = M_work.tocsr()
    C_mod = C_work.tocsr() if C_work is not None else None
    _logger.info("Dynamic BCs applied: %s DOFs constrained", bc_dofs.size)
    return K_mod, M_mod, C_mod, bc_dofs
