# pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/benscoter/linear_gbt_3D.py
"""
2-node 3D GBT beam: plane-section beam kinematics, warping modes and cross-section distortion modes.

**DOFs:** per node the six standard DOFs, then ``[b_1, ..., b_p, a_1, a_1', ..., a_m, a_m']``
(``p + 2m`` extra DOFs). ``b_j`` are the amplitudes of the warping modes ``Omega_j(s)``
(linear in ``x``): the first ``n_sv`` span the Saint-Venant warping of the section, the
rest are the shear-lag modes asked for. ``a_k`` are the amplitudes of the conventional
distortional modes (cubic Hermite in ``x``, since the Vlasov modes need ``a_k''``).
Element order: the 12 standard DOFs (node 1, node 2), then the extra DOFs of node 1,
then those of node 2.

**Strains:** ``[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x, b'(p), b(p), a''(m), a'(m), a(m)]``;
``D`` is the element-signed ``SectionAnalysis.gbt_matrices(m, p_extra)`` stiffness. The
beam strains act on plane sections; all warping is carried by ``b``. With free warping
(no condition on ``b``) and tip loads the element reproduces the Saint-Venant
(Timoshenko) beam; ``b = 0`` at a support keeps the section plane there and gives shear
lag and the warping-restraint end effect.

**Quadrature:** a full rule of at least 4 points (``a_k^2`` is degree 6); the block of
the linear shear strains ``{gamma_xy, gamma_xz, b}`` is integrated at the shear-block order
(default 1), as for the Timoshenko element, so the element does not lock.

**Mass:** ``int N_q^T M_s N_q dx`` with the coupled section mass ``M_s`` for
``[u, v, w, theta_x, theta_y, theta_z, b(p), a(m), a'(m)]``.

**Supports:** warping restraint is ``b_j = 0`` (DOFs ``Q1 .. Qp``); a diaphragm (rib)
is ``a_k = 0`` (``Q(p+1), Q(p+3), ...``), and ``a_k' = 0`` also stops the mode warping.
"""

from __future__ import annotations

import numpy as np

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)

_N_STD = 12


def gbt_modes_of_row(element_dictionary: dict, idx: int, key: str = "gbt_modes") -> int:
    modes = element_dictionary.get(key)
    return 0 if modes is None else int(np.asarray(modes).ravel()[idx])


def _hermite_x(xi: float, L: float):
    """Hermite shapes for ``(a_1, a_1', a_2, a_2')`` on ``x = (1 + xi) L / 2``: values, d/dx, d2/dx2."""
    x = 0.5 * (1.0 + xi)
    H = np.array([1 - 3 * x**2 + 2 * x**3, L * (x - 2 * x**2 + x**3), 3 * x**2 - 2 * x**3, L * (-x**2 + x**3)])
    dH = np.array([-6 * x + 6 * x**2, L * (1 - 4 * x + 3 * x**2), 6 * x - 6 * x**2, L * (-2 * x + 3 * x**2)]) / L
    d2H = np.array([-6 + 12 * x, L * (-4 + 6 * x), 6 - 12 * x, L * (-2 + 6 * x)]) / L**2
    return H, dH, d2H


class LinearGBTBeamElement3D(LinearTimoshenkoBeamElement3D):
    """Generalised Beam Theory element for thin-walled anisotropic sections.

    The number of distortional and extra shear-lag modes comes from
    ``element_dictionary["gbt_modes"]`` and ``["gbt_shear_modes"]`` (the ``[modes]`` and
    ``[shear_modes]`` columns of ``element_section.txt``); the Saint-Venant warping
    modes are always carried. ``element_dictionary["extra_dofs"]`` must give
    ``n_sv + shear_modes + 2 modes`` for the row (``attach_section_analyses`` sets it).
    """

    supports_element_frame = True    # mode amplitudes are scalars: they map to themselves

    element_type_name = "GBT-3D"
    _thin_walled_warping_dof = True

    @classmethod
    def _nodal_dof_count(cls, element_dictionary: dict, idx: int, warp_mesh: bool) -> int:
        n = gbt_modes_of_row(element_dictionary, idx, "extra_dofs")
        if n < 1:
            raise ValueError(
                "LinearGBTBeamElement3D needs element_dictionary['extra_dofs'] (set by "
                "attach_section_analyses from element_section.txt)."
            )
        return 6 + n

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        ed = kwargs["element_dictionary"]
        eid = self.element_id
        idx = int(np.where(np.asarray(ed["ids"]) == eid)[0][0])
        if not self.has_thin_walled_section:
            raise ValueError(
                f"Element {eid}: LinearGBTBeamElement3D needs a thin-walled section "
                "(laminate.txt, section_geometry.txt and element_section.txt)."
            )
        sa = self.section_analysis
        if sa.gbt_stiffness is None:
            raise ValueError(f"Element {eid}: the section has no GBT matrices.")
        self.n_modes = gbt_modes_of_row(ed, idx)
        self.n_shear = gbt_modes_of_row(ed, idx, "gbt_shear_modes")
        self.n_warping = sa.gbt_sv_modes + self.n_shear
        if self.dof_per_node - 6 != self.n_warping + 2 * self.n_modes:
            raise ValueError(
                f"Element {eid}: {self.dof_per_node - 6} extra DOFs per node, the section needs "
                f"{self.n_warping} warping + 2 x {self.n_modes} mode DOFs."
            )
        K, M = sa.gbt_matrices(self.n_modes, self.n_shear)
        self._D_gbt, self._M_gbt = K, M
        self._load_rows = sa.gbt_indices(self.n_modes, self.n_shear)[1]
        p = self.n_warping
        self._reduced_rows = np.r_[3, 4, 6 + p + np.arange(p)]

    # ------------------------------------------------------------------ #
    # DOF columns
    # ------------------------------------------------------------------ #
    @property
    def _n_extra(self) -> int:
        return self.dof_per_node - 6

    def _warping_cols(self, j: int):
        """Element columns of ``b_j`` at node 1 and node 2."""
        c1 = _N_STD + j
        return c1, c1 + self._n_extra

    def _mode_cols(self, k: int):
        """Element columns of ``(a_k, a_k')`` at node 1 and node 2."""
        e = self._n_extra
        c1 = _N_STD + self.n_warping + 2 * k
        return c1, c1 + 1, c1 + e, c1 + e + 1

    # ------------------------------------------------------------------ #
    # Operators
    # ------------------------------------------------------------------ #
    def _D_matrix_7x7(self) -> np.ndarray:
        return self._D_gbt

    def D_matrix(self) -> np.ndarray:
        return self._D_gbt

    def _full_order(self) -> int:
        return max(self._timoshenko_orders.max_full_order, 4)

    def _B_at(self, xi_g: float):
        """B (n_strain, n_dof) at ``xi_g``, with the shape functions for the Gauss cache."""
        N, dN, d2N = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B6 = self.strain_displacement_operator.physical_coordinate_form(dN, d2N, N)[0]
        m, p = self.n_modes, self.n_warping
        B = np.zeros((6 + 2 * p + 3 * m, self._n_dof))
        B[:6, :_N_STD] = B6
        L1, L2 = 0.5 * (1.0 - xi_g), 0.5 * (1.0 + xi_g)
        for j in range(p):
            c = list(self._warping_cols(j))
            B[6 + j, c] = [-1.0 / self.L, 1.0 / self.L]
            B[6 + p + j, c] = [L1, L2]
        H, dH, d2H = _hermite_x(float(xi_g), self.L)
        oa = 6 + 2 * p
        for k in range(m):
            c = list(self._mode_cols(k))
            B[oa + k, c] = d2H
            B[oa + m + k, c] = dH
            B[oa + 2 * m + k, c] = H
        return B, N, dN

    # ------------------------------------------------------------------ #
    # Stiffness and loads
    # ------------------------------------------------------------------ #
    def element_stiffness_matrix(self):
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        D = self._D_gbt
        r = self._reduced_rows
        D_red = np.zeros_like(D)
        D_red[np.ix_(r, r)] = D[np.ix_(r, r)]
        detJ = self.jacobian_determinant
        Ke = np.zeros((self._n_dof, self._n_dof))
        gauss_cache = []
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(self._full_order())):
            B, N, dN = self._B_at(xi_g)
            Ke += B.T @ (D - D_red) @ B * w_g * detJ
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g), weight=float(w_g), B_matrix=B.copy(), D_matrix=D.copy(),
                jacobian=float(detJ), shape_functions=N.copy(), shape_derivatives=dN.copy(),
            ))
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(self._timoshenko_orders.shear_block)):
            B, _, _ = self._B_at(xi_g)
            Ke += B.T @ D_red @ B * w_g * detJ
        Ke = 0.5 * (Ke + Ke.T)
        if self.logger_operator:
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} {self.element_type_name} K_e {Ke.shape} ===")
            self.logger_operator.log_matrix("stiffness", Ke, {"name": "Element Stiffness Matrix"})
            self.logger_operator.flush("stiffness")
        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=Ke,
            gauss_data=gauss_cache,
            integration_scheme="Gauss-Legendre (Selective)",
        )

    def _precurvature_equivalent_load(self) -> np.ndarray:
        if not np.any(self._E_0_voigt):
            return np.zeros(self._n_dof)
        E0 = np.zeros(self._D_gbt.shape[0])
        E0[:6] = self._E_0_voigt
        f_e = np.zeros(self._n_dof)
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(self._full_order())):
            B, _, _ = self._B_at(xi_g)
            f_e += B.T @ self._D_gbt @ E0 * w_g * self.jacobian_determinant
        return f_e

    def linear_geometric_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        raise NotImplementedError(
            "LinearGBTBeamElement3D has no geometric stiffness: its distortional modes need a wall-level "
            "K_sigma (local and distortional buckling). Use LinearBenscoterBeamElement3D for global buckling."
        )

    # ------------------------------------------------------------------ #
    # Fields, wall loads and mass
    # ------------------------------------------------------------------ #
    def _field_interpolation(self, xi: float) -> np.ndarray:
        """``N_q`` (6 + p + 2m, n_dof) for ``[u, v, w, theta_x, theta_y, theta_z, b(p), a(m), a'(m)]``."""
        m, p = self.n_modes, self.n_warping
        Nq = np.zeros((6 + p + 2 * m, self._n_dof))
        L1, L2 = 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)
        Nq[:6, :6] = L1 * np.eye(6)
        Nq[:6, 6:12] = L2 * np.eye(6)
        for j in range(p):
            Nq[6 + j, list(self._warping_cols(j))] = [L1, L2]
        H, dH, _ = _hermite_x(xi, self.L)
        for k in range(m):
            c = list(self._mode_cols(k))
            Nq[6 + p + k, c] = H
            Nq[6 + p + m + k, c] = dH
        return Nq

    def _wall_load_vector(self) -> np.ndarray:
        """``F_e = int N_q^T f dx`` for the uniform generalised wall load ``f``."""
        Fe = np.zeros(self._n_dof)
        f = self._wall_load_resultant()
        if f is None:
            return Fe
        f = np.asarray(f)[self._load_rows]
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(4)):
            Fe += self._field_interpolation(xi_g).T @ f * w_g * self.jacobian_determinant
        return Fe

    def element_mass_matrix(self):
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        M_e = np.zeros((self._n_dof, self._n_dof))
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(4)):
            Nq = self._field_interpolation(xi_g)
            M_e += Nq.T @ self._M_gbt @ Nq * w_g * self.jacobian_determinant
        return MassObject(element_id=self.element_id, element_type=self.element_type_name, M_e=0.5 * (M_e + M_e.T))
