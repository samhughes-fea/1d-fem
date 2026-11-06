# pre_processing/element_library/nonlinear/large_rotations/corotational/utilities/frame.py
"""Chord-aligned orthonormal frames for corotational beams."""

from __future__ import annotations

import numpy as np


def orthogonal_basis_from_chord(d: np.ndarray) -> np.ndarray:
    """Rotation matrix R (3, 3) with columns e1, e2, e3; e1 aligns with chord ``d``."""
    L = float(np.linalg.norm(d))
    if L <= 1e-18:
        raise ValueError("Corotational chord length is degenerate.")
    e1 = d / L
    a = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    if abs(float(np.dot(e1, a))) > 0.95:
        a = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    e2 = np.cross(e1, a)
    e2 = e2 / np.linalg.norm(e2)
    e3 = np.cross(e1, e2)
    return np.column_stack((e1, e2, e3))


def block_twelve(R: np.ndarray) -> np.ndarray:
    """Map global nodal vectors to corotated axes: four R.T blocks (trans/rot at each node)."""
    Rt = R.T
    return np.kron(np.eye(4, dtype=np.float64), Rt)


def reference_frame(d: np.ndarray) -> np.ndarray:
    """Undeformed element frame: ``e1`` along the chord, ``e2`` towards global ``y`` (or ``z``)."""
    e1 = d / np.linalg.norm(d)
    ref = np.array([0.0, 1.0, 0.0]) if abs(e1[1]) < 0.9 else np.array([0.0, 0.0, 1.0])
    e2 = ref - (ref @ e1) * e1
    e2 /= np.linalg.norm(e2)
    return np.column_stack([e1, e2, np.cross(e1, e2)])
