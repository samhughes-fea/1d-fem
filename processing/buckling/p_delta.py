# processing/buckling/p_delta.py
"""Second-order (P-delta) static solve: ``(K + K_sigma(U)) U = F`` by fixed-point iteration.

``K_sigma`` is the geometric stiffness of the current section forces
(``assemble_global_geometric_stiffness``), so each iteration is a linear solve about the
previous state. The iteration contracts at the rate of the load ratio ``P / P_cr`` and
fails to converge at or past the first buckling load. Rotations stay small: the
displacements are not updated into the geometry (that is the large-displacement
analysis).

Public API
----------
PDeltaResult, solve_p_delta.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional, Sequence

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from processing.buckling.linear_buckling import assemble_global_geometric_stiffness


@dataclass
class PDeltaResult:
    """Displacements and the iteration record of a P-delta solve."""

    U: np.ndarray
    converged: bool
    iterations: int
    history: List[float] = field(default_factory=list)  # relative change of U per iteration


def _direct_solver(K, F: np.ndarray, free: np.ndarray) -> np.ndarray:
    U = np.zeros(F.shape[0])
    if sparse.issparse(K):
        Kf = sparse.csc_matrix(K)[free][:, free]
        U[free] = spla.spsolve(Kf, F[free])
    else:
        U[free] = np.linalg.solve(np.asarray(K)[np.ix_(free, free)], F[free])
    return U


def solve_p_delta(
    K,
    F: np.ndarray,
    elements: Sequence[Any],
    fixed_dofs: Sequence[int],
    *,
    tolerance: float = 1e-8,
    max_iterations: int = 50,
    U0: Optional[np.ndarray] = None,
    solver: Optional[Callable[[Any, np.ndarray], np.ndarray]] = None,
) -> PDeltaResult:
    """Solve ``(K + K_sigma(U)) U = F`` with homogeneous constraints on ``fixed_dofs``.

    Parameters
    ----------
    K
        Global linear stiffness (dense or sparse), without boundary conditions.
    F
        Global load vector.
    elements
        Elements with ``linear_geometric_stiffness_matrix``, in assembly order.
    fixed_dofs
        Global DOFs held at zero.
    tolerance
        Convergence on ``||U_k - U_{k-1}|| / ||U_k||``.
    U0
        Starting state (default: the linear solution).
    solver
        ``solver(K_total, F) -> U`` to use instead of a direct solve on the free DOFs
        (it must apply the constraints itself).

    Returns
    -------
    PDeltaResult
        ``converged`` is False when the iteration stalls or grows, typically at or past
        the first buckling load.
    """
    F = np.asarray(F, dtype=np.float64).ravel()
    n = F.shape[0]
    free = np.setdiff1d(np.arange(n), np.asarray(fixed_dofs, dtype=int))
    solve = solver or (lambda K_t, F_t: _direct_solver(K_t, F_t, free))
    U = np.asarray(U0, dtype=np.float64).ravel() if U0 is not None else solve(K, F)
    history: List[float] = []
    for it in range(1, int(max_iterations) + 1):
        Kg = assemble_global_geometric_stiffness(list(elements), U, n)
        K_t = (sparse.csr_matrix(K) + Kg) if sparse.issparse(K) else np.asarray(K) + Kg.toarray()
        U_new = solve(K_t, F)
        scale = max(np.linalg.norm(U_new), np.finfo(float).tiny)
        change = float(np.linalg.norm(U_new - U) / scale)
        history.append(change)
        U = U_new
        if not np.isfinite(change):
            return PDeltaResult(U, False, it, history)
        if change <= tolerance:
            return PDeltaResult(U, True, it, history)
        if it >= 3 and change >= history[-2] >= history[-3]:
            return PDeltaResult(U, False, it, history)  # not contracting: at or past buckling
    return PDeltaResult(U, False, int(max_iterations), history)
