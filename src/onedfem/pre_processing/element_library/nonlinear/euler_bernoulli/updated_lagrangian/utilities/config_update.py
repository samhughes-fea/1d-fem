# Updated-Lagrangian configuration update — last converged geometry becomes the new reference.

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ConfigurationUpdateOperator:
    """Push nodal translations from ``U_e`` into the current chord endpoints.

    The runner owns *when* this runs (after a converged load increment). The
    element owns *how* the chord, ``L``, and operators are rebuilt.
    """

    def apply_translations(self, node_coords: np.ndarray, U_e: np.ndarray) -> np.ndarray:
        coords = np.asarray(node_coords, dtype=np.float64).copy()
        u = np.asarray(U_e, dtype=np.float64).ravel()
        coords[0] += u[0:3]
        coords[1] += u[6:9]
        return coords

    def chord_length(self, node_coords: np.ndarray) -> float:
        c = np.asarray(node_coords, dtype=np.float64)
        L = float(np.linalg.norm(c[1] - c[0]))
        if L <= 1e-18:
            raise ValueError("Updated-Lagrangian chord length is degenerate.")
        return L
