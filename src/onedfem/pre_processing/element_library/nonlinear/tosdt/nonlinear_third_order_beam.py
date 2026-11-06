# Shared TL pipeline for Reddy / Levinson (16 DOF Hermite + Timoshenko Green–Lagrange on 12).

from __future__ import annotations

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.tosdt.third_order_beam import (
    N_DOF,
    ThirdOrderShearBeamElement3D,
    _N_GAUSS,
)
from onedfem.pre_processing.element_library.nonlinear.tosdt.utilities import (
    GeometricStiffnessOperator,
    GreenLagrangeStrainOperator,
    StressResultantOperator,
)


class NonlinearThirdOrderShearBeamElement3D(ThirdOrderShearBeamElement3D):
    """Total Lagrangian third-order beam: linear tosdt ``B8``/``D8`` plus Timoshenko ``E_nl`` on 12 DOFs."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.green_lagrange_strain_operator = GreenLagrangeStrainOperator(
            element_length=self.L, include_shear=True
        )
        self.stress_resultant_operator = StressResultantOperator()
        self.geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._K_0_nl = None

    def _dxi_dx(self) -> float:
        return 2.0 / self.L

    def _get_K_0(self) -> np.ndarray:
        if self._K_0_nl is None:
            self._K_0_nl = np.asarray(super().element_stiffness_matrix().K_e, dtype=np.float64)
        return self._K_0_nl

    def _tl_at(self, U_e: np.ndarray, xi: float):
        B_lin = self._B_at(xi)
        N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi]))
        dxi_dx = self._dxi_dx()
        dN_dx = dN_dξ[0] * dxi_dx
        d2N_dx2 = d2N_dξ2[0] * (dxi_dx**2)
        U12 = U_e[:12]
        E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx, d2N_dx2, U12, N[0])
        B_nl12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
            dN_dx, d2N_dx2, U12, N[0]
        )
        E = B_lin @ U_e
        E[:6] += E_nl
        B_nl = np.zeros_like(B_lin)
        B_nl[:6, :12] = B_nl12
        return B_lin, B_nl, E, dN_dx, d2N_dx2, N[0]

    def _test_and_D(self, B: np.ndarray, D: np.ndarray):
        if self.variational:
            return B, D
        return B[:6], D[:6, :]

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        U_e = np.asarray(U_e, dtype=np.float64).reshape(N_DOF)
        D = self.D_matrix()
        detJ = self.jacobian_determinant
        F = np.zeros(N_DOF, dtype=np.float64)
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(_N_GAUSS)):
            B_lin, B_nl, E, *_ = self._tl_at(U_e, float(xi_g))
            B_tot = B_lin + B_nl
            T, Dt = self._test_and_D(B_tot, D)
            F += T.T @ (Dt @ E) * w_g * detJ
        return F

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        U_e = np.asarray(U_e, dtype=np.float64).reshape(N_DOF)
        D = self.D_matrix()
        detJ = self.jacobian_determinant
        K_0 = self._get_K_0()
        K_delta = np.zeros((N_DOF, N_DOF), dtype=np.float64)
        K_sigma = np.zeros((N_DOF, N_DOF), dtype=np.float64)
        op = self.green_lagrange_strain_operator
        for xi_g, w_g in zip(*np.polynomial.legendre.leggauss(_N_GAUSS)):
            B_lin, B_nl, E, dN_dx, d2N_dx2, N = self._tl_at(U_e, float(xi_g))
            B_tot = B_lin + B_nl
            if self.variational:
                K_delta += (B_tot.T @ D @ B_tot - B_lin.T @ D @ B_lin) * w_g * detJ
                S6 = (D @ E)[:6]
            else:
                Dt = D[:6, :]
                T_tot, T_lin = B_tot[:6], B_lin[:6]
                K_delta += (T_tot.T @ Dt @ B_tot - T_lin.T @ Dt @ B_lin) * w_g * detJ
                S6 = (Dt @ E)[:6]
            G = np.zeros((12, 12), dtype=np.float64)
            for j in range(12):
                e = np.zeros(12, dtype=np.float64)
                e[j] = 1.0
                G[:, j] = op.nonlinear_strain_displacement_gradient(dN_dx, d2N_dx2, e, N).T @ S6
            K_sigma[:12, :12] += G * w_g * detJ
        return K_0 + K_delta + K_sigma
