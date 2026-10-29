# pre_processing/element_library/linear/beam/second_order.py
"""Second-order (geometric) stiffness of the linear beam elements from their section forces.

The section is rigid in its own plane, so a point at ``(y, z)`` moves by
``u_y = v - z phi`` and ``u_z = w + y phi`` (right-handed nodal ``v``, ``w`` and
``phi = theta_x`` of the reference axis). The quadratic part of the Green strain is

    eps_xx^(2)   = 1/2 [(v' - z phi')^2 + (w' + y phi')^2]
    gamma_xs^(2) = (v' - z phi')(-t_z phi) + (w' + y phi')(t_y phi)

along a wall with tangent ``t`` (and the same across it with the normal ``n``). Its work
against the current section stresses is ``1/2 g^T S g`` per length with
``g = [v', w', phi', phi]`` and

    S = [[ N,    0,   -S_z, -V_z],
         [ 0,    N,    S_y,  V_y],
         [-S_z,  S_y,  W,    Q_r],
         [-V_z,  V_y,  Q_r,  0  ]]

where ``[N, S_y, S_z, W, V_y, V_z, Q_r]`` are the section terms of
``SectionAnalysis.geometric_terms`` (``S_z = M_y``, ``S_y = -M_z``,
``W = int sigma (y^2 + z^2)``). Every term is about the reference axis, so the
shear-centre and tension-centre offsets come from the coupled section stiffness and
from the stresses themselves, not from shifted formulas. The Wagner term ``W`` includes
the polar term, the Wagner coefficients of the bending moments and, for Benscoter
sections, the bimoment's warping stress. ``K_sigma = int G^T S G dx`` with ``G`` the
element's interpolation of ``g``.

Public API
----------
GEOMETRIC_TERMS, geometric_section_matrix, GeometricStiffnessMixin.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

GEOMETRIC_TERMS: Tuple[str, ...] = ("N", "S_y", "S_z", "W", "V_y", "V_z", "Q_r")

#: Gauss points for K_sigma: exact for Hermite slopes against linear section forces
_GAUSS_ORDER = 4


def geometric_section_matrix(c: np.ndarray) -> np.ndarray:
    """``S`` (4, 4) on ``[v', w', phi', phi]`` from the section terms ``c`` (7,)."""
    N, Sy, Sz, W, Vy, Vz, Qr = (float(v) for v in c)
    return np.array([
        [N, 0.0, -Sz, -Vz],
        [0.0, N, Sy, Vy],
        [-Sz, Sy, W, Qr],
        [-Vz, Vy, Qr, 0.0],
    ])


def _linear_rows(L: float, xi: float, n_dof: int) -> np.ndarray:
    """``G`` (4, n_dof) for the two-node Timoshenko family (linear fields, node-major DOFs).

    The slopes are the section rotation plus the element's constant shear strain,
    ``v' = theta_z + gamma_xy`` and ``w' = -theta_y + gamma_xz`` with ``gamma`` taken at the
    centre (as the reduced shear rule sees it). That is the chord slope plus the linear
    part of the rotation, the slope of the quadratic ("linked") deflection the element
    implies, and it converges much faster in lateral buckling than the chord slope alone.
    """
    G = np.zeros((4, n_dof))
    G[0, [1, 7]] = [-1.0 / L, 1.0 / L]
    G[0, [5, 11]] = [-0.5 * xi, 0.5 * xi]
    G[1, [2, 8]] = [-1.0 / L, 1.0 / L]
    G[1, [4, 10]] = [0.5 * xi, -0.5 * xi]
    G[2, [3, 9]] = [-1.0 / L, 1.0 / L]
    G[3, [3, 9]] = [0.5 * (1.0 - xi), 0.5 * (1.0 + xi)]
    return G


class GeometricStiffnessMixin:
    """``linear_geometric_stiffness_matrix`` for an element that defines ``_geometric_state``.

    ``_geometric_state(xi, U_e)`` returns ``(G, c)``: the rows (4, n_dof) of
    ``[v', w', phi', phi]`` and the section terms (7,) at natural coordinate ``xi``.
    """

    def _geometric_state(self, xi: float, U_e: np.ndarray):  # pragma: no cover - interface
        raise NotImplementedError

    def linear_geometric_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """Geometric stiffness ``K_sigma`` (n_dof, n_dof) about the state ``U_e``.

        The section forces come from the linear strains of ``U_e`` (less any
        precurvature), so ``K_sigma`` is linear in ``U_e``. See
        ``docs/element_library/second_order.md``.
        """
        U_e = np.asarray(U_e, dtype=np.float64).ravel()
        if U_e.size != self._n_dof:
            raise ValueError(f"expected {self._n_dof} element DOFs, got {U_e.size}")
        detJ = self.jacobian_determinant
        K = np.zeros((self._n_dof, self._n_dof))
        for xi, wg in zip(*np.polynomial.legendre.leggauss(_GAUSS_ORDER)):
            G, c = self._geometric_state(float(xi), U_e)
            K += G.T @ geometric_section_matrix(c) @ G * (wg * detJ)
        return 0.5 * (K + K.T)

    def section_geometric_terms(self, U_e: np.ndarray, xi: float = 0.0) -> np.ndarray:
        """Section terms ``[N, S_y, S_z, W, V_y, V_z, Q_r]`` at ``xi`` for the state ``U_e``."""
        return self._geometric_state(float(xi), np.asarray(U_e, dtype=np.float64).ravel())[1]

    def _geometric_terms_matrix(self) -> np.ndarray:
        H = self.section_analysis.geometric_terms
        if H is None:
            raise ValueError(f"Element {self.element_id}: the section has no second-order terms")
        return np.asarray(H, dtype=np.float64)
