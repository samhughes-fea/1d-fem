# pre_processing\element_library\truss\utilities\local_frame.py

"""
Local coordinate frame and mapping matrix L for 3D truss (axial + transverse + torsion).
"""

import numpy as np


def direction_cosines_and_transverse(node_coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute axial direction cosines and one transverse direction (local y) for the element.

    Parameters
    ----------
    node_coords : np.ndarray, shape (2, 3)
        Coordinates of node 1 and node 2 (rows).

    Returns
    -------
    axial : np.ndarray, shape (3,)
        Unit vector along element (x2 - x1) / L. Direction cosines (cx, cy, cz).
    transverse : np.ndarray, shape (3,)
        Unit vector perpendicular to axial (local y). Used for transverse shear DOF.
    """
    d = node_coords[1] - node_coords[0]
    L = np.linalg.norm(d)
    if L < 1e-12:
        raise ValueError("Element length too small for local frame")
    axial = d / L

    # Reference vector to build transverse: use global Z unless axial is parallel to Z
    ref = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    if np.abs(np.dot(axial, ref)) > 0.99:
        ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    transverse = np.cross(axial, ref)
    tnorm = np.linalg.norm(transverse)
    if tnorm < 1e-12:
        ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        transverse = np.cross(axial, ref)
        tnorm = np.linalg.norm(transverse)
    transverse = transverse / tnorm
    return axial, transverse


def build_L_matrix_6x12(axial: np.ndarray, transverse: np.ndarray) -> np.ndarray:
    """
    Build the 6×12 matrix L that maps global 12 DOF to local 6 DOF.

    Local DOF order: (u1_axial, u2_axial), (v1_trans, v2_trans), (θx1, θx2).
    Global DOF order: node1 [ux, uy, uz, θx, θy, θz], node2 [ux, uy, uz, θx, θy, θz].

    So u_local = L @ u_global (6,) = L (6×12) @ u_global (12,).

    Returns
    -------
    L : np.ndarray, shape (6, 12)
    """
    cx, cy, cz = axial
    ty_x, ty_y, ty_z = transverse
    L = np.zeros((6, 12), dtype=np.float64)
    # Row 0: u1_axial = cx*ux1 + cy*uy1 + cz*uz1
    L[0, 0:3] = (cx, cy, cz)
    # Row 1: u2_axial = cx*ux2 + cy*uy2 + cz*uz2
    L[1, 6:9] = (cx, cy, cz)
    # Row 2: v1_trans = transverse component at node 1
    L[2, 0:3] = (ty_x, ty_y, ty_z)
    # Row 3: v2_trans at node 2
    L[3, 6:9] = (ty_x, ty_y, ty_z)
    # Row 4: θx1 -> global DOF 3
    L[4, 3] = 1.0
    # Row 5: θx2 -> global DOF 9
    L[5, 9] = 1.0
    return L
