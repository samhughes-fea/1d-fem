# pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/benscoter/linear_gbt_3D.py
"""
2-node 3D GBT beam: the Benscoter beam plus ``m`` cross-section distortional modes.

**DOFs:** per node the six standard DOFs, then ``[psi, a_1, a_1', ..., a_m, a_m']``
(``1 + 2m`` extra DOFs). Element order: the 12 standard DOFs (node 1, node 2), then the
extra DOFs of node 1, then those of node 2. Mode amplitudes ``a_k(x)`` are cubic Hermite
in ``x`` (the Vlasov modes need ``a_k''``).

**Strains:** ``[ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x, d, ψ′, a''(m), a'(m), a(m)]``; ``D`` is the
element-signed ``SectionAnalysis.gbt_matrices(m)`` stiffness. With all ``a_k = 0`` the
element is the Benscoter element.

**Quadrature:** a full rule of at least 4 points (``a_k²`` is degree 6), with the
``{γ_xy, γ_xz, d}`` block at the shear-block order as in the Benscoter element.

**Mass:** ``∫ N_qᵀ M_s N_q dx`` with the coupled section mass ``M_s`` for
``[u, v, w, θ_x, θ_y, θ_z, ψ, a, a′]`` (rigid fields and ``ψ`` linear, ``a`` Hermite).

**Supports:** a diaphragm (rib) at a node is ``a_k = 0`` there (DOFs ``Q2, Q4, …``);
``a_k′ = 0`` also stops the mode warping (``Q3, Q5, …``).
"""

from __future__ import annotations

import numpy as np

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.benscoter.linear_benscoter_3D import (
    LinearBenscoterBeamElement3D,
)

_N_STD = 12
_REDUCED_ROWS = [3, 4, 6]


def gbt_modes_of_row(element_dictionary: dict, idx: int) -> int:
    modes = element_dictionary.get("gbt_modes")
    return 0 if modes is None else int(np.asarray(modes).ravel()[idx])


def _hermite_x(xi: float, L: float):
    """Hermite shapes for ``(a_1, a_1', a_2, a_2')`` on ``x = (1 + xi) L / 2``: values, d/dx, d2/dx2."""
    x = 0.5 * (1.0 + xi)
    H = np.array([1 - 3 * x**2 + 2 * x**3, L * (x - 2 * x**2 + x**3), 3 * x**2 - 2 * x**3, L * (-x**2 + x**3)])
    dH = np.array([-6 * x + 6 * x**2, L * (1 - 4 * x + 3 * x**2), 6 * x - 6 * x**2, L * (-2 * x + 3 * x**2)]) / L
    d2H = np.array([-6 + 12 * x, L * (-4 + 6 * x), 6 - 12 * x, L * (-2 + 6 * x)]) / L**2
    return H, dH, d2H


class LinearGBTBeamElement3D(LinearBenscoterBeamElement3D):
    """Benscoter beam with conventional GBT cross-section distortion modes.

    The number of modes comes from ``element_dictionary["gbt_modes"]`` (set from the
    ``[modes]`` column of ``element_section.txt``).
    """

    element_type_name = "GBT-3D"

    @classmethod
    def _nodal_dof_count(cls, element_dictionary: dict, idx: int, warp_mesh: bool) -> int:
        return 7 + 2 * gbt_modes_of_row(element_dictionary, idx)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.n_modes = (self.dof_per_node - 7) // 2
        if self.n_modes < 1:
            raise ValueError(
                f"Element {self.element_id}: LinearGBTBeamElement3D needs at least one cross-section mode "
                "([modes] column of element_section.txt)."
            )
        K, M = self.section_analysis.gbt_matrices(self.n_modes)
        self._D_gbt, self._M_gbt = K, M

    # ------------------------------------------------------------------ #
    def _D_matrix_7x7(self) -> np.ndarray:
        return self._D_gbt

    def _full_order(self) -> int:
        return max(self._timoshenko_orders.max_full_order, 4)

    def _mode_cols(self, k: int):
        """Element columns of ``(a_k, a_k')`` at node 1 and node 2."""
        e = self._n_extra
        c1 = _N_STD + 1 + 2 * k
        return c1, c1 + 1, c1 + e, c1 + e + 1

    def _B_matrix_7x14(self, dN_dξ, d2N_dξ2, N, xi=None):
        B8 = super()._B_matrix_7x14(dN_dξ, d2N_dξ2, N, xi=xi)
        n_gauss = B8.shape[0]
        xi = np.zeros(n_gauss) if xi is None else np.atleast_1d(xi)
        m = self.n_modes
        B = np.zeros((n_gauss, 8 + 3 * m, self._n_dof))
        B[:, :8, :] = B8
        for g in range(n_gauss):
            H, dH, d2H = _hermite_x(float(xi[g]), self.L)
            for k in range(m):
                cols = list(self._mode_cols(k))
                B[g, 8 + k, cols] = d2H
                B[g, 8 + m + k, cols] = dH
                B[g, 8 + 2 * m + k, cols] = H
        return B

    def linear_geometric_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        raise NotImplementedError("LinearGBTBeamElement3D is linear only (second order arrives later).")

    # ------------------------------------------------------------------ #
    def _field_interpolation(self, xi: float) -> np.ndarray:
        """``N_q`` (7 + 2m, n_dof) for ``[u, v, w, θ_x, θ_y, θ_z, ψ, a(m), a'(m)]`` at ``xi``."""
        m = self.n_modes
        Nq = np.zeros((7 + 2 * m, self._n_dof))
        L1, L2 = 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)
        Nq[:6, :6] = L1 * np.eye(6)
        Nq[:6, 6:12] = L2 * np.eye(6)
        p1, p2 = self._psi_cols
        Nq[6, p1], Nq[6, p2] = L1, L2
        H, dH, _ = _hermite_x(xi, self.L)
        for k in range(m):
            cols = list(self._mode_cols(k))
            Nq[7 + k, cols] = H
            Nq[7 + m + k, cols] = dH
        return Nq

    def element_mass_matrix(self):
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        M_e = np.zeros((self._n_dof, self._n_dof))
        detJ = self.jacobian_determinant
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(4)):
            Nq = self._field_interpolation(xi_g)
            M_e += Nq.T @ self._M_gbt @ Nq * w_g * detJ
        return MassObject(element_id=self.element_id, element_type=self.element_type_name, M_e=0.5 * (M_e + M_e.T))
