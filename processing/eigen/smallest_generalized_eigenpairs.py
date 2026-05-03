# processing/eigen/smallest_generalized_eigenpairs.py
"""
Smallest eigenpairs for the symmetric generalized pencil **K** x = λ **M** x.

Shared by §2 eigen-style workflows (including harmonic modal superposition) and
``simulation_runner.spectral.VibrationBucklingBackend`` so sparse retries stay aligned.
"""

from __future__ import annotations

import logging
from typing import Tuple

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

logger = logging.getLogger(__name__)


def smallest_generalized_eigenpairs(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    num_modes: int,
    *,
    dense_threshold: int = 512,
    context: str = "generalized eigenproblem",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Lowest ``num_modes`` eigenvalues λ and eigenvectors (columns) for **K** x = λ **M** x.

    Parameters
    ----------
    K_mod, M_mod
        Symmetric sparse matrices (CSR), same shape.
    num_modes
        Number of modes (clamped to ``[1, n]``).
    dense_threshold
        Use dense ``eigh`` when ``n <= dense_threshold``.
    context
        Prefix for log messages.

    Returns
    -------
    lam : ndarray
        Shape ``(num_modes,)``, sorted ascending.
    Phi : ndarray
        Shape ``(n, num_modes)``.
    """
    n = int(K_mod.shape[0])
    if n < 2:
        raise ValueError(f"{context}: need at least 2 global DOFs")
    num_modes = max(1, min(int(num_modes), n))
    if n <= dense_threshold:
        Kd = K_mod.toarray()
        Md = M_mod.toarray()
        Kd = 0.5 * (Kd + Kd.T)
        Md = 0.5 * (Md + Md.T)
        lam, Phi = eigh(Kd, Md)
        idx = np.argsort(lam)[:num_modes]
        return lam[idx].astype(np.float64), np.asarray(Phi[:, idx], dtype=np.float64)

    k = min(num_modes, n - 1)
    ncv = max(min(n - 1, max(2 * k + 1, 20)), k + 3)
    maxit = max(1000, n * 10)

    def _eigsh_sm() -> tuple[np.ndarray, np.ndarray]:
        return eigsh(
            K_mod,
            k=k,
            M=M_mod,
            which="SM",
            maxiter=maxit,
            ncv=ncv,
        )

    def _eigsh_shift_invert(sigma: float) -> tuple[np.ndarray, np.ndarray]:
        return eigsh(
            K_mod,
            k=k,
            M=M_mod,
            sigma=float(sigma),
            which="LM",
            maxiter=max(maxit, 3000),
            ncv=max(ncv, k * 4),
        )

    lam: np.ndarray | None = None
    Phi: np.ndarray | None = None
    try:
        lam, Phi = _eigsh_sm()
    except Exception as exc:
        logger.warning(
            "%s: eigsh(which='SM') failed (%s); retrying shift-invert near zero",
            context,
            exc,
        )
        diag = np.abs(K_mod.diagonal())
        sigma0 = float(np.median(diag[diag > 0.0])) * 1e-8 if np.any(diag > 0.0) else 1e-8
        try:
            lam_si, Phi_si = _eigsh_shift_invert(sigma0)
            take = min(num_modes, int(lam_si.shape[0]))
            idx = np.argsort(np.real(lam_si))[:take]
            lam = lam_si[idx]
            Phi = Phi_si[:, idx]
        except Exception as exc2:
            logger.warning(
                "%s: eigsh shift-invert failed (%s); dense fallback when feasible",
                context,
                exc2,
            )
            if n > 4096:
                raise RuntimeError(
                    f"{context}: sparse eigen solve failed and system is too large for dense "
                    "fallback; try a smaller mesh, fewer modes, or improve conditioning "
                    "(boundary conditions / scaling)"
                ) from exc2
            Kd = K_mod.toarray()
            Md = M_mod.toarray()
            Kd = 0.5 * (Kd + Kd.T)
            Md = 0.5 * (Md + Md.T)
            lam_full, Phi_full = eigh(Kd, Md)
            idx = np.argsort(lam_full)[:num_modes]
            lam = lam_full[idx]
            Phi = Phi_full[:, idx]

    assert lam is not None and Phi is not None
    return np.asarray(lam, dtype=np.float64), np.asarray(Phi, dtype=np.float64)
