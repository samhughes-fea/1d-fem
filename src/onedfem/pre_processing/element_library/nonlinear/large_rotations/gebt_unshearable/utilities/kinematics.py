# Unshearable GEBT: stretch and nodal curvature from the Battini chord map.

from __future__ import annotations

from typing import Tuple

import numpy as np

from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.utilities.local_map import (
    battini_p_batch,
)

_H2 = 1e-4


def extras_none(U: np.ndarray, p: np.ndarray) -> np.ndarray:
    return np.zeros((U.shape[0], 0), dtype=U.dtype)


def local_p(U: np.ndarray, X1: np.ndarray, X2: np.ndarray, E0: np.ndarray, l0: float, eid: int) -> np.ndarray:
    U = np.asarray(U, dtype=np.float64).reshape(1, -1)
    return battini_p_batch(U, X1, X2, E0, l0, extras_none, eid)[0]


def local_p_batch(U: np.ndarray, X1: np.ndarray, X2: np.ndarray, E0: np.ndarray, l0: float, eid: int) -> np.ndarray:
    return battini_p_batch(np.asarray(U), X1, X2, E0, l0, extras_none, eid)


def jacobian_p(U_e: np.ndarray, X1, X2, E0, l0, eid) -> np.ndarray:
    n = U_e.size
    h = 1e-30
    Z = np.tile(U_e.astype(np.complex128), (n, 1))
    Z[np.arange(n), np.arange(n)] += 1j * h
    return (local_p_batch(Z, X1, X2, E0, l0, eid).imag / h).T


def hessian_contraction(U_e: np.ndarray, f: np.ndarray, X1, X2, E0, l0, eid) -> np.ndarray:
    n = U_e.size
    hs = np.full(12, _H2)
    hs[[0, 1, 2, 6, 7, 8]] *= l0
    j, k = np.triu_indices(12)
    m = j.size
    D = np.zeros((4, m, n))
    for c, (sj, sk) in enumerate(((1, 1), (1, -1), (-1, 1), (-1, -1))):
        D[c, np.arange(m), j] += sj * hs[j]
        D[c, np.arange(m), k] += sk * hs[k]
    vals = (local_p_batch(U_e[None, :] + D.reshape(-1, n), X1, X2, E0, l0, eid) @ f).reshape(4, m)
    h2 = (vals[0] - vals[1] - vals[2] + vals[3]) / (4 * hs[j] * hs[k])
    Kg = np.zeros((n, n))
    Kg[j, k] = h2
    Kg[k, j] = h2
    return Kg


def local_force(p: np.ndarray, K_core: np.ndarray, f_E0: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    f = K_core @ p - f_E0
    return f, K_core
