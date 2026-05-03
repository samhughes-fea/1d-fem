# pre_processing/element_library/point_load_utils.py
"""Helpers for optional harmonic per-load phase (column 10 in point_load rows)."""

from __future__ import annotations

import numpy as np


def point_load_phase_rad(load_row: np.ndarray) -> float:
    """Return phase [rad] from optional 10th column; default 0."""
    lr = np.asarray(load_row, dtype=np.float64).ravel()
    if lr.size > 9:
        return float(lr[9])
    return 0.0


def add_phased_increment(Fe: np.ndarray, increment: np.ndarray, phase_rad: float) -> np.ndarray:
    """``Fe += increment * exp(i * phase_rad)``, promoting *Fe* to complex when needed."""
    inc = np.asarray(increment, dtype=np.float64).ravel()
    if abs(float(phase_rad)) < 1e-20 and Fe.dtype != np.complex128:
        return Fe + inc
    out = Fe.astype(np.complex128) if Fe.dtype != np.complex128 else Fe
    return out + inc.astype(np.complex128) * np.exp(1j * float(phase_rad))


def accumulate_phased_slice(
    Fe: np.ndarray,
    indices: np.ndarray | list[int],
    increment: np.ndarray,
    phase_rad: float,
) -> np.ndarray:
    """
    ``Fe[indices] += increment`` with optional complex phase on the increment.

    Promotes *Fe* to ``complex128`` when *phase_rad* is non-zero.
    """
    inc = np.asarray(increment, dtype=np.float64).ravel()
    idx = np.asarray(indices, dtype=np.intp).ravel()
    if abs(phase_rad) < 1e-20:
        Fe[idx] = Fe[idx] + inc
        return Fe
    if Fe.dtype != np.complex128:
        Fe = Fe.astype(np.complex128)
    Fe[idx] = Fe[idx] + inc.astype(np.complex128) * np.exp(1j * float(phase_rad))
    return Fe
