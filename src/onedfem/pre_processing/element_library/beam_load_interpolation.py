"""Work-conjugate load interpolation for two-node, six-DOF-per-node beams.

A load ``P = [F_x, F_y, F_z, M_x, M_y, M_z]`` at ``xi`` (global, right-handed, like the
nodal DOFs) does virtual work ``dU^T N_f(xi)^T P``, where ``N_f`` (6, 12) returns the
element's own displacement and rotation fields at ``xi`` from ``U_e``. So the nodal
loads are ``N_f^T P`` for a point load and ``int N_f^T q dx`` for a distributed one.

- Euler-Bernoulli (``bending="hermite"``): ``v`` and ``w`` are Hermite cubics, the
  rotations are their slopes, ``theta_z = v'`` and ``theta_y = -w'``. A moment at a node
  goes straight to that node's rotation; a uniform line load gives the consistent
  ``q L / 2`` forces and ``q L^2 / 12`` end moments.
- Timoshenko (``bending="linear"``): every field is linear in ``xi`` and the rotations
  are independent fields.

Axial displacement and twist are linear in both. Element shape-function tensors are
not used here: they carry the element's own sign conventions (``theta_y`` negated) and
the Euler-Bernoulli rotation rows hold Hermite displacement shapes, not rotations.
"""

from __future__ import annotations

import numpy as np

__all__ = ["beam_field_interpolation", "point_load_vector", "distributed_load_vector"]


def _hermite(xi: float, L: float):
    """Hermite shapes and their x-slopes for ``[d1, slope1, d2, slope2]``."""
    H = np.array([
        0.25 * (1 - xi) ** 2 * (2 + xi),
        L / 8 * (1 - xi) ** 2 * (1 + xi),
        0.25 * (1 + xi) ** 2 * (2 - xi),
        -L / 8 * (1 + xi) ** 2 * (1 - xi),
    ])
    dH = np.array([
        0.75 * (xi**2 - 1),
        L / 8 * (3 * xi**2 - 2 * xi - 1),
        0.75 * (1 - xi**2),
        L / 8 * (3 * xi**2 + 2 * xi - 1),
    ]) * (2.0 / L)
    return H, dH


def beam_field_interpolation(xi: float, L: float, bending: str = "hermite") -> np.ndarray:
    """``N_f`` (6, 12): ``[u, v, w, theta_x, theta_y, theta_z](xi) = N_f U_e``, right-handed."""
    xi = float(xi)
    L1, L2 = 0.5 * (1 - xi), 0.5 * (1 + xi)
    N = np.zeros((6, 12))
    for r in (0, 3):
        N[r, [r, r + 6]] = (L1, L2)
    if bending == "linear":
        for r in (1, 2, 4, 5):
            N[r, [r, r + 6]] = (L1, L2)
        return N
    if bending != "hermite":
        raise ValueError(f"bending must be 'hermite' or 'linear', got {bending!r}")
    H, dH = _hermite(xi, L)
    v_cols = [1, 5, 7, 11]           # v1, theta_z1, v2, theta_z2   (v' = theta_z)
    w_cols = [2, 4, 8, 10]           # w1, theta_y1, w2, theta_y2   (w' = -theta_y)
    sign_w = np.array([1.0, -1.0, 1.0, -1.0])
    N[1, v_cols] = H
    N[5, v_cols] = dH
    N[2, w_cols] = H * sign_w
    N[4, w_cols] = -dH * sign_w
    return N


def point_load_vector(xi: float, L: float, load6, bending: str = "hermite") -> np.ndarray:
    """Nodal loads (12,) of a point load ``load6`` at ``xi``."""
    return beam_field_interpolation(xi, L, bending).T @ np.asarray(load6)


def distributed_load_vector(xi_g, w_g, q_g, L: float, bending: str = "hermite") -> np.ndarray:
    """Nodal loads (12,) of line loads ``q_g`` (n_g, 6) at Gauss points ``xi_g``, weights ``w_g``."""
    F = np.zeros(12, dtype=np.result_type(float, np.asarray(q_g)))
    for xi, w, q in zip(xi_g, w_g, q_g):
        F += beam_field_interpolation(xi, L, bending).T @ q * (w * L / 2)
    return F
