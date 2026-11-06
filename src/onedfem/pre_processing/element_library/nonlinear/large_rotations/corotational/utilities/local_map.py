# Battini–Pacoste local deformational map p(U) for the 12 standard DOFs.

from __future__ import annotations

from typing import Callable

import numpy as np

from .rotation import exp_batch, log_batch


def vector_norm(a: np.ndarray) -> np.ndarray:
    """sqrt(a . a) along the last axis (complex-step safe)."""
    return np.sqrt(np.sum(a * a, axis=-1))


def battini_p_batch(
    U: np.ndarray,
    X1: np.ndarray,
    X2: np.ndarray,
    E0: np.ndarray,
    l0: float,
    extras_fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    element_id: int,
) -> np.ndarray:
    """Local deformational DOFs ``p`` (m, n) from global states ``U`` (m, n).

    Columns 0–11 follow the Battini–Pacoste chord map; columns 12+ come from ``extras_fn(U, p)``.
    """
    d = (X2 + U[:, 6:9]) - (X1 + U[:, 0:3])
    l = vector_norm(d)
    if np.any(np.abs(l) <= 1e-12 * l0):
        raise ValueError(f"Element {element_id}: corotational chord collapsed")
    R1 = exp_batch(U[:, 3:6]) @ E0
    R2 = exp_batch(U[:, 9:12]) @ E0
    e1 = d / l[:, None]
    q = 0.5 * (R1[:, :, 1] + R2[:, :, 1])
    e3 = np.cross(e1, q)
    e3 = e3 / vector_norm(e3)[:, None]
    Rr = np.stack([e1, np.cross(e3, e1), e3], axis=2)
    RrT = np.transpose(Rr, (0, 2, 1))
    p = np.zeros(U.shape, dtype=U.dtype)
    p[:, 3:6] = log_batch(RrT @ R1)
    p[:, 6] = l - l0
    p[:, 9:12] = log_batch(RrT @ R2)
    p[:, 12:] = extras_fn(U, p)
    return p
