# processing/dynamic/time_integration.py
"""
Time integration for M u'' + C u' + K u = F(t).
Newmark-beta (constant average acceleration) by default; no imports from processing.static.
"""

import logging
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from typing import Callable, Optional, Tuple

_logger = logging.getLogger(__name__)


def newmark_integrate(
    K: csr_matrix,
    M: csr_matrix,
    C: Optional[csr_matrix],
    u0: np.ndarray,
    v0: np.ndarray,
    t_grid: np.ndarray,
    F_func: Callable[[float], np.ndarray],
    beta: float = 0.25,
    gamma: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Newmark-beta time integration: M u'' + C u' + K u = F(t).

    Parameters
    ----------
    K, M : csr_matrix
        Stiffness and mass (modified by BCs).
    C : csr_matrix or None
        Damping (optional); if None, undamped.
    u0, v0 : np.ndarray
        Initial displacement and velocity (length n).
    t_grid : np.ndarray
        Time points (length n_steps + 1).
    F_func : callable
        F(t) returns force vector (length n).
    beta, gamma : float
        Newmark parameters; beta=0.25, gamma=0.5 gives constant average acceleration.

    Returns
    -------
    U : np.ndarray
        (n_steps + 1, n) displacements.
    V : np.ndarray
        (n_steps + 1, n) velocities.
    A : np.ndarray
        (n_steps + 1, n) accelerations.
    """
    n = M.shape[0]
    if K.shape != (n, n) or M.shape != (n, n):
        raise ValueError("K and M must be square and same size as u0")
    if u0.shape != (n,) or v0.shape != (n,):
        raise ValueError("u0 and v0 must have length n")
    if C is not None and C.shape != (n, n):
        raise ValueError("C must be (n, n) or None")

    n_steps = len(t_grid) - 1
    dt = t_grid[1] - t_grid[0] if n_steps > 0 else 0.0
    if np.any(np.abs(np.diff(t_grid) - dt) > 1e-12 * (dt + 1)):
        _logger.warning("Time grid may not be uniform; using variable dt per step")

    U = np.zeros((n_steps + 1, n), dtype=np.float64)
    V = np.zeros((n_steps + 1, n), dtype=np.float64)
    A = np.zeros((n_steps + 1, n), dtype=np.float64)

    U[0] = u0
    V[0] = v0
    if C is not None:
        A[0] = spsolve(M, F_func(t_grid[0]) - K @ u0 - C @ v0)
    else:
        A[0] = spsolve(M, F_func(t_grid[0]) - K @ u0)

    C_mat = C if C is not None else csr_matrix((n, n), dtype=np.float64)

    for i in range(n_steps):
        dt_i = t_grid[i + 1] - t_grid[i]
        t_half = t_grid[i] + 0.5 * dt_i
        F_next = F_func(t_grid[i + 1])

        a1 = 1.0 / (beta * dt_i * dt_i)
        a2 = 1.0 / (beta * dt_i)
        a3 = 1.0 / (2 * beta) - 1.0
        a4 = gamma / (beta * dt_i) - 1.0
        a5 = dt_i * (gamma / beta - 1.0)

        K_eff = K + a1 * M + a4 * C_mat
        R_eff = (
            F_next
            + M @ (a1 * U[i] + a2 * V[i] + a3 * A[i])
            + C_mat @ (a4 * U[i] + a5 * V[i] + dt_i * (gamma / (2 * beta) - 1.0) * A[i])
        )
        U[i + 1] = spsolve(K_eff, R_eff)
        A[i + 1] = a1 * (U[i + 1] - U[i]) - a2 * V[i] - a3 * A[i]
        V[i + 1] = V[i] + dt_i * ((1 - gamma) * A[i] + gamma * A[i + 1])

    return U, V, A
