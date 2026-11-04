# pre_processing\element_library\bar\utilities\local_frame.py

"""
Axial direction cosines and 4×12 L matrix for 3D bar (axial + torsion only).
"""

import numpy as np


def direction_cosines(node_coords: np.ndarray) -> np.ndarray:
    """Unit vector (cx, cy, cz) along element."""
    d = node_coords[1] - node_coords[0]
    L = np.linalg.norm(d)
    if L < 1e-12:
        raise ValueError("Element length too small")
    return d / L


def build_L_matrix_4x12(axial: np.ndarray) -> np.ndarray:
    """
    Build 4×12 L: local DOF (u1_axial, u2_axial, θx1, θx2) from global 12 DOF.
    u_local = L @ u_global.
    """
    cx, cy, cz = axial
    L = np.zeros((4, 12), dtype=np.float64)
    L[0, 0:3] = (cx, cy, cz)
    L[1, 6:9] = (cx, cy, cz)
    L[2, 3] = 1.0
    L[3, 9] = 1.0
    return L
