# processing/harmonic/damping_table.py
"""Frequency-dependent modal damping ratio ζ(f) from tabulated data."""

from __future__ import annotations

import os

import numpy as np


def load_zeta_vs_frequency_hz(path: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Load a two-column text table: frequency_hz, zeta.

    Lines starting with ``#`` and blank lines are ignored. Values must be
    monotonically increasing in frequency (validated).
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"damping zeta table not found: {path}")
    rows: list[tuple[float, float]] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            rows.append((float(parts[0]), float(parts[1])))
    if len(rows) < 2:
        raise ValueError(f"damping zeta table needs at least two samples, got {len(rows)} in {path}")
    f_hz = np.array([r[0] for r in rows], dtype=np.float64)
    zeta = np.array([r[1] for r in rows], dtype=np.float64)
    if np.any(zeta < 0.0):
        raise ValueError("damping table zeta values must be non-negative")
    if np.any(np.diff(f_hz) <= 0.0):
        raise ValueError("damping table frequencies_hz must be strictly increasing")
    return f_hz, zeta


def interpolate_zeta_hz(sample_hz: np.ndarray, f_table_hz: np.ndarray, zeta_table: np.ndarray) -> np.ndarray:
    """
    Piecewise linear interpolation of ζ(f); clamps outside the table range to endpoints.
    """
    x = np.asarray(sample_hz, dtype=np.float64).ravel()
    out = np.interp(x, f_table_hz, zeta_table, left=float(zeta_table[0]), right=float(zeta_table[-1]))
    return out.astype(np.float64)
