# processing/modal/buckling.py
"""Linear buckling: global :math:`\\mathbf{K}_\\sigma` from prestress and eigenpairs for :math:`\\mathbf{K}_L + \\lambda \\mathbf{K}_\\sigma`."""

from __future__ import annotations

import logging
from typing import Any, List, Optional, Sequence, Tuple

import numpy as np
from scipy import linalg
from scipy.sparse import coo_matrix, csr_matrix

from processing.modal.boundary_conditions import apply_boundary_conditions
from processing.modal.assembly import _compute_local_global_dof_map, _scatter_element_matrix

logger = logging.getLogger(__name__)


def assemble_global_geometric_stiffness(
    elements: List[Any],
    U_global: np.ndarray,
    total_dof: int,
) -> csr_matrix:
    """
    Assemble :math:`\\mathbf{K}_g` (global geometric / stress stiffness) from element
    :meth:`linear_geometric_stiffness_matrix` and a reference displacement field.
    """
    if not elements:
        raise ValueError("No elements for geometric stiffness assembly")
    U = np.asarray(U_global, dtype=np.float64).ravel()
    if U.size != total_dof:
        raise ValueError(f"U_global size {U.size} != total_dof {total_dof}")

    local_global_dof_map = _compute_local_global_dof_map(elements, total_dof)
    rows_all: List[np.ndarray] = []
    cols_all: List[np.ndarray] = []
    data_all: List[np.ndarray] = []

    for elem, dof_map in zip(elements, local_global_dof_map):
        U_e = U[dof_map]
        if not hasattr(elem, "linear_geometric_stiffness_matrix"):
            raise TypeError(
                f"Element {elem.element_id} ({type(elem).__name__}) has no linear_geometric_stiffness_matrix — "
                "use linear Euler–Bernoulli or Timoshenko beams for buckling."
            )
        Kg = elem.linear_geometric_stiffness_matrix(U_e)
        r, c, d = _scatter_element_matrix(Kg, dof_map)
        rows_all.append(r)
        cols_all.append(c)
        data_all.append(d)

    Kg_global = coo_matrix(
        (
            np.concatenate(data_all),
            (np.concatenate(rows_all), np.concatenate(cols_all)),
        ),
        shape=(total_dof, total_dof),
        dtype=np.float64,
    ).tocsr()
    logger.info("Geometric stiffness assembled: nnz=%s", Kg_global.nnz)
    return Kg_global


def apply_buckling_boundary_conditions(
    K_global: csr_matrix,
    Kg_global: csr_matrix,
    prescribed_displacements: Optional[dict] = None,
    fixed_dofs: Optional[Sequence[int]] = None,
) -> Tuple[csr_matrix, csr_matrix, np.ndarray]:
    """
    Apply the **same** constrained DOFs as modal vibration: penalty on ``K``,
    zero rows/columns on ``Kg`` at constrained DOFs (no artificial geometric stiffness at supports).
    """
    n = K_global.shape[0]
    dummy_M = csr_matrix((n, n), dtype=np.float64)
    K_mod, _, bc_dofs = apply_boundary_conditions(
        K_global, dummy_M, fixed_dofs=fixed_dofs, prescribed_displacements=prescribed_displacements
    )
    if bc_dofs.size == 0:
        return K_global.copy(), Kg_global.copy(), bc_dofs

    Kg_work = Kg_global.astype(np.float64).tolil()
    for dof in bc_dofs:
        Kg_work[dof, :] = 0.0
        Kg_work[:, dof] = 0.0
        Kg_work[dof, dof] = 0.0
    Kg_mod = Kg_work.tocsr()
    return K_mod, Kg_mod, bc_dofs


def solve_linear_buckling_eigenpairs(
    K_mod: csr_matrix,
    Kg_mod: csr_matrix,
    num_modes: int,
    constrained_dofs: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve the generalized eigenproblem :math:`\\mathbf{K}\\mathbf{x} = \\mu (-\\mathbf{K}_g)\\mathbf{x}`
    associated with :math:`(\\mathbf{K} + \\lambda \\mathbf{K}_g)\\mathbf{x}=\\mathbf{0}`.

    Parameters
    ----------
    constrained_dofs
        Global indices where modal BCs were applied (penalty on ``K`` and zeroed rows/cols on ``Kg``).
        The eigenproblem is assembled on the **free** subspace so ``-K_g`` is not singular at fixed DOFs.

    Returns
    -------
    lambdas : np.ndarray
        Buckling load factors :math:`\\lambda_i = \\mu_i` (smallest positive controls first instability).
    vectors : np.ndarray
        Mode shapes (columns), same ordering as ``lambdas`` (length ``n`` global DOFs).
    """
    num_modes = max(int(num_modes), 1)
    n = K_mod.shape[0]
    if constrained_dofs is not None and np.asarray(constrained_dofs).size > 0:
        cd = np.unique(np.asarray(constrained_dofs, dtype=np.int32).ravel())
        mask = np.ones(n, dtype=bool)
        mask[cd] = False
        free = np.nonzero(mask)[0]
    else:
        free = np.arange(n, dtype=np.int32)

    if free.size == 0:
        raise RuntimeError("Buckling eigenproblem has no free DOFs after boundary conditions.")

    Kd = K_mod[np.ix_(free, free)].toarray()
    Gd = Kg_mod[np.ix_(free, free)].toarray()
    rhs = -Gd
    # Symmetrize for numerical stability
    Kd = 0.5 * (Kd + Kd.T)
    rhs = 0.5 * (rhs + rhs.T)

    # ``eigh`` requires a definite metric; :math:`-\\mathbf{K}_g` may be indefinite under compression.
    evals_c, evecs_free = linalg.eig(Kd, rhs)
    evals = np.asarray(evals_c, dtype=np.complex128)
    imag_tol = 1e-9 * np.maximum(1.0, np.abs(evals.real))
    stable = np.abs(evals.imag) <= imag_tol
    evals = np.real(evals[stable])
    evecs_free = np.real(evecs_free[:, stable])
    positive = evals > 0
    if not np.any(positive):
        raise RuntimeError(
            "No positive buckling eigenvalues — check prestress (compression) and boundary conditions."
        )
    idx_pos = np.where(positive)[0]
    order = idx_pos[np.argsort(evals[idx_pos])][:num_modes]
    lambdas = evals[order]
    vecs_free = evecs_free[:, order]
    vectors = np.zeros((n, vecs_free.shape[1]), dtype=np.float64)
    vectors[free, :] = vecs_free
    return lambdas, vectors
