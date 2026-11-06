# pre_processing/element_library/nonlinear/timoshenko/total_lagrangian/nonlinear_timoshenko_3D.py
"""
2-node 3D Timoshenko beam, Total Lagrangian geometric nonlinearity.

Linear operators, selective ``K_0``, loads, mass and warping come from
``LinearTimoshenkoBeamElement3D``. This module adds Green–Lagrange ``K_delta``,
Hessian ``K_sigma``, and ``F_int``.

**See Also:** ``docs/element_library/total_lagrangian_beam_formulation.md``.
"""

from typing import List, Tuple

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.utilities.green_lagrange_strain import (
    GreenLagrangeStrainOperator,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.utilities.stress_resultant import (
    StressResultantOperator,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.utilities.tl_green_lagrange_voigt import (
    chord_frame_green_lagrange_voigt_timoshenko_12,
)

_NL_TS_W_STD = 12
_NL_TS_W_STRAIN = 7
_NL_TS_W_DOF = 14


def _B_nl_6x12_to_7x14(B_nl_6_12: np.ndarray) -> np.ndarray:
    """Pad TL nonlinear ``B`` from (6, 12) to (7, 14); warping row/cols remain zero."""
    B = np.zeros((_NL_TS_W_STRAIN, _NL_TS_W_DOF), dtype=np.float64)
    B[:6, :_NL_TS_W_STD] = B_nl_6_12
    return B


class NonlinearTimoshenkoBeamElement3D(LinearTimoshenkoBeamElement3D):
    """TL Timoshenko: selective ``K_0`` plus ``K_delta`` and Hessian ``K_sigma``.

    Newton ``K_sigma`` is ``_geometric_stiffness_at``; ``GeometricStiffnessOperator`` is
    the beam-column / prestress path.
    """

    element_type_name = "Timoshenko-3D-Nonlinear"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.green_lagrange_strain_operator = GreenLagrangeStrainOperator(
            element_length=self.L, include_shear=True
        )
        self.stress_resultant_operator = StressResultantOperator()
        self.geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._K_0: np.ndarray | None = None

    @property
    def dxi_dx(self) -> float:
        return 2.0 / self.L

    def _tl_voigt_strain_at_gauss(self, U_e: np.ndarray, xi_g: float) -> np.ndarray:
        """Six-component Voigt strain ``E`` at ``xi_g`` (GESDB may override)."""
        return chord_frame_green_lagrange_voigt_timoshenko_12(self, U_e, xi_g)

    def _gp_tl_kinematics(
        self, U_e: np.ndarray, xi_g: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        U_use = np.asarray(U_e, dtype=np.float64).ravel()
        if U_use.size >= _NL_TS_W_DOF:
            U_use = U_use[:_NL_TS_W_STD]
        N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        E = self._tl_voigt_strain_at_gauss(U_use, xi_g)
        dN_dx = dN_dξ[0] * self.dxi_dx
        return N, dN_dξ, d2N_dξ2, B, E, dN_dx

    def _geometric_stiffness_at(
        self, dN_dx: np.ndarray, d2N_dx2: np.ndarray, N: np.ndarray, S: np.ndarray
    ) -> np.ndarray:
        """``sum_k S_k d²E_k/dU²`` at one point (12, 12)."""
        G = np.zeros((12, 12), dtype=np.float64)
        op = self.green_lagrange_strain_operator
        for j in range(12):
            e = np.zeros(12, dtype=np.float64)
            e[j] = 1.0
            G[:, j] = op.nonlinear_strain_displacement_gradient(dN_dx, d2N_dx2, e, N).T @ S[:6]
        return G

    def _material_rules(self, D: np.ndarray) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """``[(xi, w, D_rest), (xi_s, w_s, D_shear)]`` matching linear selective ``K_0``."""
        D_shear = np.zeros_like(D)
        D_shear[3:5, 3:5] = D[3:5, 3:5]
        xi, w = self.integration_points
        xi_s, w_s = np.polynomial.legendre.leggauss(self._timoshenko_orders.shear_block)
        return [(xi, w, D - D_shear), (xi_s, w_s, D_shear)]

    def _get_K_0(self) -> np.ndarray:
        if self._K_0 is None:
            self._K_0 = np.asarray(super().element_stiffness_matrix().K_e, dtype=np.float64)
        return self._K_0

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """``K_T = K_0 + K_delta(U_e) + K_sigma(U_e)``."""
        if self._n_dof == _NL_TS_W_DOF:
            return self._tangent_stiffness_matrix_warping(U_e)
        K_0 = self._get_K_0()
        D = self.material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        U_e = np.asarray(U_e, dtype=np.float64).reshape(12)
        K_sigma = np.zeros((12, 12), dtype=np.float64)
        K_delta = np.zeros((12, 12), dtype=np.float64)
        for xi_r, w_r, D_part in self._material_rules(D):
            for xi_g, w_g in zip(xi_r, w_r):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B_lin = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
                dN_dx = dN_dξ[0] * self.dxi_dx
                d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
                B_nl = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                    dN_dx, d2N_dx2, U_e, N[0]
                )
                B_tot = B_lin + B_nl
                K_delta += (B_tot.T @ D_part @ B_tot - B_lin.T @ D_part @ B_lin) * w_g * detJ
                E = self._tl_voigt_strain_at_gauss(U_e, xi_g)
                S = D_part @ (E - self._E_0_voigt)
                K_sigma += self._geometric_stiffness_at(dN_dx, d2N_dx2, N[0], S) * w_g * detJ
        return K_0 + K_delta + K_sigma

    def _tangent_stiffness_matrix_warping(self, U_e: np.ndarray) -> np.ndarray:
        U_e = np.asarray(U_e, dtype=np.float64).reshape(_NL_TS_W_DOF)
        U12 = U_e[:_NL_TS_W_STD]
        K_0 = self._get_K_0()
        D7 = self._D_matrix_7x7()
        detJ = self.jacobian_determinant
        K_sigma = np.zeros((_NL_TS_W_DOF, _NL_TS_W_DOF), dtype=np.float64)
        K_delta = np.zeros((_NL_TS_W_DOF, _NL_TS_W_DOF), dtype=np.float64)
        for xi_r, w_r, D_part in self._material_rules(D7):
            for xi_g, w_g in zip(xi_r, w_r):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B_lin_7 = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
                dN_dx = dN_dξ[0] * self.dxi_dx
                d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
                B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                    dN_dx, d2N_dx2, U12, N[0]
                )
                B_nl_7 = _B_nl_6x12_to_7x14(B_nl_12)
                B_tot = B_lin_7 + B_nl_7
                K_delta += (B_tot.T @ D_part @ B_tot - B_lin_7.T @ D_part @ B_lin_7) * w_g * detJ
                E6 = self._tl_voigt_strain_at_gauss(U12, xi_g)
                S6 = D_part[:6, :6] @ (E6 - self._E_0_voigt)
                K_sigma[:_NL_TS_W_STD, :_NL_TS_W_STD] += (
                    self._geometric_stiffness_at(dN_dx, d2N_dx2, N[0], S6) * w_g * detJ
                )
        return K_0 + K_delta + K_sigma

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        if self._n_dof == _NL_TS_W_DOF:
            return self._internal_force_vector_warping(U_e)
        D = self.material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        F_int = np.zeros(12, dtype=np.float64)
        for xi_r, w_r, D_part in self._material_rules(D):
            for xi_g, w_g in zip(xi_r, w_r):
                N, dN_dξ, d2N_dξ2, B_lin, E, dN_dx = self._gp_tl_kinematics(U_e, xi_g)
                d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
                B_nl = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                    dN_dx, d2N_dx2, U_e, N[0]
                )
                S = D_part @ (E - self._E_0_voigt)
                F_int += ((B_lin + B_nl).T @ S) * w_g * detJ
        return F_int

    def _internal_force_vector_warping(self, U_e: np.ndarray) -> np.ndarray:
        U_e = np.asarray(U_e, dtype=np.float64).reshape(_NL_TS_W_DOF)
        U12 = U_e[:_NL_TS_W_STD]
        D7 = self._D_matrix_7x7()
        E0_7 = np.concatenate([self._E_0_voigt, np.zeros(1, dtype=np.float64)])
        detJ = self.jacobian_determinant
        F_int = np.zeros(_NL_TS_W_DOF, dtype=np.float64)
        for xi_r, w_r, D_part in self._material_rules(D7):
            for xi_g, w_g in zip(xi_r, w_r):
                N, dN_dξ, d2N_dξ2, B_lin, E6, dN_dx = self._gp_tl_kinematics(U12, xi_g)
                d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
                B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                    dN_dx, d2N_dx2, U12, N[0]
                )
                B_lin_7 = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
                B_nl_7 = _B_nl_6x12_to_7x14(B_nl_12)
                E_warp = float(np.dot(B_lin_7[6], U_e))
                E = np.concatenate([E6, np.array([E_warp], dtype=np.float64)])
                S = D_part @ (E - E0_7)
                F_int += ((B_lin_7 + B_nl_7).T @ S) * w_g * detJ
        return F_int

    def strain_at_gauss_points(self, U_e: np.ndarray) -> List[np.ndarray]:
        U = np.asarray(U_e, dtype=np.float64).ravel()
        xi, w = self.integration_points
        result: List[np.ndarray] = []
        if U.size == _NL_TS_W_DOF:
            for xi_g, w_g in zip(xi, w):
                N, dN_dξ, d2N_dξ2, _, E6, _ = self._gp_tl_kinematics(U[:_NL_TS_W_STD], xi_g)
                B_lin_7 = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
                E_warp = float(np.dot(B_lin_7[6], U))
                result.append(np.concatenate([E6, np.array([E_warp], dtype=np.float64)]))
            return result
        for xi_g, w_g in zip(xi, w):
            _, _, _, _, E, _ = self._gp_tl_kinematics(U_e, xi_g)
            result.append(np.asarray(E, dtype=np.float64))
        return result
