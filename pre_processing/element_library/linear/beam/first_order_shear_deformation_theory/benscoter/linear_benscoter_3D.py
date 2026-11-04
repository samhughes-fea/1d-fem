# pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/benscoter/linear_benscoter_3D.py
"""
2-node 3D Benscoter beam: Timoshenko kinematics plus an independent warping amplitude ``psi``.

**Tensors:** U_e (14,) with the 12 standard DOFs node-major (u_x, u_y, u_z, θ_x, θ_y, θ_z) and
``psi`` at DOFs 12 (node 1) and 13 (node 2); K_e (14, 14). Per Gauss point B (8, 14), D (8, 8),
ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x, d, ψ′] with ``d = ψ − θ_x′`` and ``ψ′ = dψ/dx``.

**Constitutive:** ``D`` is the element-signed ``SectionAnalysis.benscoter_stiffness``. Its leading
6×6 block is the Saint-Venant section stiffness; ``D[6, 6]`` is the warping shear stiffness and
``D[7, 7]`` the primary warping stiffness ``E Γ``.

**Quadrature:** a full rule at ``max_full_order`` gives every term, then the block of the linear
strains ``{γ_xy, γ_xz, d}`` is replaced by its value at ``shear_block`` order (default 1). With
linear ``ψ`` and ``θ_x`` the ``d`` row has the same locking mechanism as the shear rows.

**Mass:** rigid-section mass of the SectionAnalysis on the standard DOFs plus
``∫ m_ω N_ψ N_ψ dx`` on ``ψ`` (``m_ω = ∫ ρ_A ω² ds``).

Requires a thin-walled section (``element_section.txt``) and a warping mesh (``[warping] = 1``).
Root warping restraint is ``psi = 0`` (DOF ``W``/``CHI``); a free end needs no condition
(zero bimoment is natural).
"""

from __future__ import annotations

import numpy as np

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)

_N_STD = 12
_N_DOF = 14
_N_STRAIN = 8
#: strain rows integrated with the reduced (shear) rule: gamma_xy, gamma_xz, d
_REDUCED_ROWS = [3, 4, 6]


class LinearBenscoterBeamElement3D(LinearTimoshenkoBeamElement3D):
    """Timoshenko beam with Benscoter restrained warping for thin-walled anisotropic sections.

    Notes
    -----
    The warping displacement of the section is the Saint-Venant warping of the six
    beam strains plus ``ω(s) d`` with ``d = ψ − θ_x′``, so the total torsional warping
    amplitude is ``ψ``. Wall strains gain ``ε_xx += ω ψ′`` and ``γ_xs += ω_s d``. With free
    warping ``ψ`` follows ``θ_x′`` and the element reproduces the Timoshenko element; a
    restrained end (``ψ = 0``) gives the boundary layer of length ``1/k`` with
    ``k² = S GJ / ((GJ + S) E Γ)``.
    """

    element_type_name = "Benscoter-3D"
    _thin_walled_warping_dof = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        eid = self.element_id
        if not self.has_thin_walled_section:
            raise ValueError(
                f"Element {eid}: LinearBenscoterBeamElement3D needs a thin-walled section "
                "(laminate.txt, section_geometry.txt and element_section.txt)."
            )
        if self._n_dof < _N_DOF or not self._warp_stiff:
            raise ValueError(
                f"Element {eid}: LinearBenscoterBeamElement3D needs the warping DOF; set [warping] = 1."
            )
        if self.section_analysis.benscoter_stiffness is None:
            raise ValueError(
                f"Element {eid}: the section has no Benscoter stiffness (a 2D section mesh, or a thin-walled "
                "section that cannot shear, e.g. a flat strip)."
            )

    # ------------------------------------------------------------------ #
    # Operators
    # ------------------------------------------------------------------ #
    @property
    def _n_extra(self) -> int:
        """Extra DOFs per node used by this element (1 for Benscoter: ``psi``)."""
        return self.dof_per_node - 6

    @property
    def _psi_cols(self) -> tuple:
        """Element columns of ``psi`` at node 1 and node 2."""
        return (_N_STD, _N_STD + self._n_extra)
    def _D_matrix_7x7(self) -> np.ndarray:  # noqa: D401  (parent hook name)
        """Element-signed Benscoter stiffness (8, 8)."""
        return np.array(self.section_analysis.benscoter_stiffness, dtype=np.float64)

    def D_matrix(self) -> np.ndarray:
        return self._D_matrix_7x7()

    def _B_matrix_7x14(self, dN_dξ: np.ndarray, d2N_dξ2: np.ndarray, N: np.ndarray, xi=None) -> np.ndarray:
        """B (n_gauss, 8, 14). ``xi`` (n_gauss,) is needed for the ``d`` row; defaults to the centre."""
        n_gauss = dN_dξ.shape[0]
        B6 = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)
        xi = np.zeros(n_gauss) if xi is None else np.atleast_1d(np.asarray(xi, dtype=np.float64))
        p1, p2 = self._psi_cols
        B = np.zeros((n_gauss, _N_STRAIN, self._n_dof), dtype=np.float64)
        B[:, :6, :_N_STD] = B6
        for g in range(n_gauss):
            B[g, 6, :] = -B[g, 5, :]
            B[g, 6, p1] += 0.5 * (1.0 - xi[g])
            B[g, 6, p2] += 0.5 * (1.0 + xi[g])
            B[g, 7, p1] = -1.0 / self.L
            B[g, 7, p2] = 1.0 / self.L
        return B

    def _full_order(self) -> int:
        return self._timoshenko_orders.max_full_order

    def _B_at(self, xi_g: float):
        N, dN, d2N = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        return self._B_matrix_7x14(dN, d2N, N, xi=np.array([xi_g]))[0], N, dN

    # ------------------------------------------------------------------ #
    # Stiffness
    # ------------------------------------------------------------------ #
    def element_stiffness_matrix(self):
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        D = self._D_matrix_7x7()
        D_red = np.zeros_like(D)
        D_red[np.ix_(_REDUCED_ROWS, _REDUCED_ROWS)] = D[np.ix_(_REDUCED_ROWS, _REDUCED_ROWS)]
        # rows 0-7 are the Benscoter strains; any further rows are integrated with the full rule
        detJ = self.jacobian_determinant
        orders = self._timoshenko_orders

        Ke = np.zeros((self._n_dof, self._n_dof), dtype=np.float64)
        gauss_cache = []
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(self._full_order())):
            B, N, dN = self._B_at(xi_g)
            Ke += B.T @ (D - D_red) @ B * w_g * detJ
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g), weight=float(w_g), B_matrix=B.copy(), D_matrix=D.copy(),
                jacobian=float(detJ), shape_functions=N.copy(), shape_derivatives=dN.copy(),
            ))
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(orders.shear_block)):
            B, _, _ = self._B_at(xi_g)
            Ke += B.T @ D_red @ B * w_g * detJ
        Ke = 0.5 * (Ke + Ke.T)

        if self.logger_operator:
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} {self.element_type_name} K_e {Ke.shape} ===")
            self.logger_operator.log_matrix("stiffness", D, {"name": f"Section stiffness D {D.shape}"})
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
            return np.zeros(self._n_dof, dtype=np.float64)
        D = self._D_matrix_7x7()
        E0 = np.zeros(D.shape[0])
        E0[:6] = self._E_0_voigt
        detJ = self.jacobian_determinant
        f_e = np.zeros(self._n_dof, dtype=np.float64)
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(self._full_order())):
            B, _, _ = self._B_at(xi_g)
            f_e += B.T @ D @ E0 * w_g * detJ
        return f_e

    def _wall_load_vector(self) -> np.ndarray:
        Fe = super()._wall_load_vector()
        f = self._wall_load_resultant()
        if f is not None:
            p1, p2 = self._psi_cols
            Fe[p1] += 0.5 * self.L * f[6]
            Fe[p2] += 0.5 * self.L * f[6]
        return Fe

    def _geometric_state(self, xi: float, U_e: np.ndarray):
        """Timoshenko rows of ``[v', w', phi', phi]`` and the section terms with the bimoment.

        The walls carry the Saint-Venant stresses of ``K psi6`` plus the warping stresses
        of ``[d, psi']``; the ``psi'`` part is the bimoment's axial stress, which adds to
        the Wagner term.
        """
        from pre_processing.element_library.linear.beam.second_order import _linear_rows

        sa = self.section_analysis
        e = self._B_at(xi)[0] @ U_e
        e[_REDUCED_ROWS] = self._B_at(0.0)[0][_REDUCED_ROWS] @ U_e  # as the reduced rule sees them
        e[:6] -= self._E_0_voigt
        c = self._geometric_terms_matrix() @ (sa.stiffness @ e[:6]) + sa.benscoter_geometric_terms @ e[6:8]
        return _linear_rows(self.L, xi, self._n_dof), c

    # ------------------------------------------------------------------ #
    # Mass
    # ------------------------------------------------------------------ #
    def element_mass_matrix(self):
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        M_e = np.zeros((_N_DOF, _N_DOF), dtype=np.float64)
        if self.section_analysis.mass:
            M_e[:_N_STD, :_N_STD] = self._section_mass_matrix()
        m_w = float(self.section_analysis.warping_inertia)
        # consistent linear-Lagrange mass on psi
        M_e[12:, 12:] = m_w * self.L / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
        return MassObject(element_id=self.element_id, element_type=self.element_type_name, M_e=M_e)
