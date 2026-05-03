# processing/harmonic/modal_superposition.py
"""Undamped modal basis + damped modal harmonic response (§4 optional path)."""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.sparse import csr_matrix

from processing.eigen.smallest_generalized_eigenpairs import smallest_generalized_eigenpairs


def undamped_natural_modes(
    K_mod: csr_matrix,
    M_mod: csr_matrix,
    num_modes: int,
    *,
    dense_threshold: int = 512,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Lowest undamped natural frequencies (rad/s) and mode shapes (columns), M-normalized.

    Uses dense ``eigh`` for ``n <= dense_threshold``, otherwise sparse ``eigsh``
    with shift-invert retry — same pencil solver as ``VibrationBucklingBackend``
    (see ``processing.eigen.smallest_generalized_eigenpairs``).
    """
    lam, Phi = smallest_generalized_eigenpairs(
        K_mod,
        M_mod,
        num_modes,
        dense_threshold=dense_threshold,
        context="harmonic modal basis",
    )
    omega = np.sqrt(np.maximum(np.real(lam), 0.0))
    return omega.astype(np.float64), np.asarray(Phi, dtype=np.float64)


def harmonic_displacement_modal_superposition(
    omega_natural_rad: np.ndarray,
    mode_shapes: np.ndarray,
    F: np.ndarray,
    frequencies_hz: np.ndarray,
    zeta: float,
    *,
    zeta_per_mode: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Steady-state response using classical viscous modal damping.

    Each mode r uses denominator (omega_r^2 - omega^2 + 2 i zeta_r omega_r omega).
    *mode_shapes* columns are M-orthonormal (from generalized eigenproblem).
    *F* may be complex (harmonic load phasing).

    When *zeta_per_mode* is set (length ``num_modes``), it overrides scalar *zeta* for each
    oscillator (e.g. ζ(f_r) sampled from a damping table at natural frequencies in Hz).
    """
    Fv = np.asarray(F, dtype=np.complex128).ravel()
    Phi = np.asarray(mode_shapes, dtype=np.float64)
    n, nm = Phi.shape
    if Fv.shape[0] != n:
        raise ValueError("F length must match mode row dimension")
    if omega_natural_rad.shape[0] != nm:
        raise ValueError("omega_natural_rad must have length num_modes")
    if zeta_per_mode is not None:
        zp = np.asarray(zeta_per_mode, dtype=np.float64).ravel()
        if zp.shape[0] != nm:
            raise ValueError("zeta_per_mode must have length num_modes")
        if np.any(zp < 0.0):
            raise ValueError("zeta_per_mode values must be non-negative")

    f_h = np.asarray(frequencies_hz, dtype=np.float64).ravel()
    nfreq = f_h.size
    Fm = Phi.T @ Fv
    U = np.zeros((n, nfreq), dtype=np.complex128)
    z0 = float(zeta)
    for k in range(nfreq):
        w = 2.0 * np.pi * float(f_h[k])
        w2 = w * w
        q = np.zeros(nm, dtype=np.complex128)
        for r in range(nm):
            wr = float(omega_natural_rad[r])
            zr = float(zeta_per_mode[r]) if zeta_per_mode is not None else z0
            den = wr * wr - w2 + 2.0j * zr * wr * w
            if abs(den) < 1e-30:
                q[r] = 0.0
            else:
                q[r] = Fm[r] / den
        U[:, k] = Phi @ q
    return U


def harmonic_truncation_metrics_vs_direct(
    U_modal: np.ndarray,
    U_direct: np.ndarray,
) -> tuple[float, float]:
    """
    Relative error norms comparing modal and direct displacement matrices (same layout).

    Returns
    -------
    mean_relative_l2 : float
        Mean over frequency columns of ``‖u_m - u_d‖₂ / max(‖u_d‖₂, ε)``.
    max_relative_l2 : float
        Maximum of the same ratio over columns.
    """
    Um = np.asarray(U_modal, dtype=np.complex128)
    Ud = np.asarray(U_direct, dtype=np.complex128)
    if Um.shape != Ud.shape:
        raise ValueError("U_modal and U_direct must have the same shape")
    eps = 1e-30
    ratios = []
    for j in range(Um.shape[1]):
        a = Um[:, j]
        b = Ud[:, j]
        nb = np.linalg.norm(b)
        ratios.append(float(np.linalg.norm(a - b) / max(nb, eps)))
    arr = np.asarray(ratios, dtype=np.float64)
    return float(np.mean(arr)), float(np.max(arr))
