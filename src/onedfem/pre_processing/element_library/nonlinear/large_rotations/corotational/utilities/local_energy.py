# Cubic/quartic local energy on corotational deformational DOFs p.

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple

import numpy as np


def cubic_quartic_force_and_tangent(
    p: np.ndarray,
    K_core: np.ndarray,
    C: np.ndarray,
    axial_pts: Optional[np.ndarray],
    G_gauss: Sequence[np.ndarray],
    gauss_w: Iterable[float],
    l0: float,
    f_E0: Optional[np.ndarray] = None,
    K_sigma0: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Local force ``f`` and Hessian ``H`` of ``½ pᵀ K p + ½ pᵀ K_σ(p) p + Π₄``.

    Optional constitutive precurvature: ``f -= f_E0`` (``Bᵀ D E_0``) and
    ``K_σ`` of the ``E_0`` prestress (``K_sigma0``), assembled with ``E_0`` held
    out of the linear map ``C``.
    """
    Ks = np.tensordot(p, C, axes=(0, 0))
    M = C @ p
    f = K_core @ p + Ks @ p + 0.5 * np.einsum("jkl,k,l->j", C, p, p)
    H = K_core + Ks + M + M.T
    if K_sigma0 is not None:
        f = f + K_sigma0 @ p
        H = H + K_sigma0
    if f_E0 is not None:
        f = f - np.asarray(f_E0, dtype=np.float64).ravel()[: f.size]
    if axial_pts is None:
        return f, H
    y, z, k = axial_pts[:, 0], axial_pts[:, 1], axial_pts[:, 2]
    detJ = 0.5 * l0
    for G, wg in zip(G_gauss, gauss_w):
        v1, w1, f1, _ = G @ p
        a = v1 - z * f1
        b = w1 + y * f1
        e2 = 0.5 * (a * a + b * b)
        Ja = np.stack([a, b, -z * a + y * b], axis=1)
        grad = (k * e2) @ Ja
        Hg = Ja.T @ (k[:, None] * Ja)
        ke = k * e2
        Hg += np.array(
            [
                [ke.sum(), 0.0, -(ke * z).sum()],
                [0.0, ke.sum(), (ke * y).sum()],
                [-(ke * z).sum(), (ke * y).sum(), (ke * (y * y + z * z)).sum()],
            ]
        )
        G3 = G[:3]
        f = f + wg * detJ * G3.T @ grad
        H = H + wg * detJ * G3.T @ Hg @ G3
    return f, H
