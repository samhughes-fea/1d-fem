"""
Directional modal metrics without importing processing.static.

Uses a unit pattern on translational +Z DOFs. Any node-major layout with the six standard
DOFs first works (6, 7 with warping, more with cross-section modes; see ``dof_layout``).
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix


def global_uz_unit_pattern(total_dof: int, dof_per_node: int) -> np.ndarray:
    """Return a unit-norm vector activating uz on every node (zero if unsupported layout)."""
    if dof_per_node < 6 or total_dof % dof_per_node != 0:
        return np.zeros(total_dof, dtype=np.float64)
    r = np.zeros(total_dof, dtype=np.float64)
    r[2::dof_per_node] = 1.0
    nrm = float(np.linalg.norm(r))
    if nrm < 1e-30:
        return r
    return r / nrm


def modal_effective_mass_fraction_z(
    M_mod: csr_matrix,
    mode_shapes: np.ndarray,
    *,
    dof_per_node: int,
) -> np.ndarray | None:
    """
    Per mode *j*, scalar ((φ_jᵀ M r)²) / ((rᵀ M r)(φ_jᵀ M φ_j))) with *r* = global uz unit pattern.

    Returns None if the mesh layout is not node-major with at least six DOFs per node.
    """
    phi = np.asarray(mode_shapes, dtype=np.float64)
    if phi.ndim != 2:
        return None
    n, nm = phi.shape
    r = global_uz_unit_pattern(n, dof_per_node)
    if float(np.linalg.norm(r)) < 1e-30:
        return None
    Mr = M_mod @ r
    rMr = float(np.dot(r, Mr))
    if abs(rMr) < 1e-30:
        return None
    out = np.empty(nm, dtype=np.float64)
    for j in range(nm):
        v = phi[:, j]
        mjj = float(np.dot(v, M_mod @ v))
        if abs(mjj) < 1e-30:
            out[j] = 0.0
            continue
        lj = float(np.dot(v, Mr))
        out[j] = (lj * lj) / (rMr * mjj)
    return out
