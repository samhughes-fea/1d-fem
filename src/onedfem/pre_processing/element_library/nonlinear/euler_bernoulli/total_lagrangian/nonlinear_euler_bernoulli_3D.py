# pre_processing/element_library/nonlinear/euler_bernoulli/total_lagrangian/nonlinear_euler_bernoulli_3D.py
"""
2-node 3D Euler-Bernoulli beam, Total Lagrangian (TL) geometric nonlinearity.

Linear operators, loads, mass and warping come from
``LinearEulerBernoulliBeamElement3D``. This module adds Green–Lagrange Gauss
assembly for ``F_int`` and ``K_T``.

**See Also:** ``docs/element_library/total_lagrangian_beam_formulation.md``.
"""

from typing import List

import numpy as np

from onedfem.pre_processing.element_library.linear.beam.zosdt.euler_bernoulli.linear_euler_bernoulli_3D import (
    LinearEulerBernoulliBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.total_lagrangian.utilities import (
    GeometricStiffnessOperator,
    GreenLagrangeStrainOperator,
    StressResultantOperator,
)

_NL_W_STD = 12
_NL_W_DOF = 14


def _B_nl_12_to_14(B_nl_6_12: np.ndarray) -> np.ndarray:
    B = np.zeros((7, _NL_W_DOF), dtype=np.float64)
    B[:6, :_NL_W_STD] = B_nl_6_12
    return B


class NonlinearEulerBernoulliBeamElement3D(LinearEulerBernoulliBeamElement3D):
    """TL EB: Green–Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums."""

    #: The tangent and internal force are not rotated by ``framed_element``: refuse rotated or
    #: curved meshes rather than inherit the linear parent's frame support.
    supports_element_frame = False
    bending_interpolation = "hermite"

    element_type_name = "Euler-Bernoulli-3D-Nonlinear"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.green_lagrange_strain_operator = GreenLagrangeStrainOperator(element_length=self.L)
        self.stress_resultant_operator = StressResultantOperator()
        self.geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)
        self._K_0: np.ndarray | None = None

    def _get_K_0(self) -> np.ndarray:
        """Linear material stiffness (cached); 12-DOF path only."""
        if self._n_dof != 12:
            raise NotImplementedError("_get_K_0 is defined for the 12-DOF EB path only.")
        if self._K_0 is None:
            self._K_0 = np.asarray(super().element_stiffness_matrix().K_e, dtype=np.float64)
        return self._K_0

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """``K_T = K_mat(U_e) + K_sigma(U_e)`` with ``B_tot = B_lin + B_nl``."""
        if self._n_dof == 12:
            U_e = np.asarray(U_e, dtype=np.float64).reshape(12)
            D = self.material_stiffness_operator.assembly_form()
            xi, w = self.integration_points
            detJ = self.jacobian_determinant
            dξ_dx = 2.0 / self.L
            d2ξ_dx2 = 4.0 / (self.L ** 2)
            n_g = len(xi)
            N_gp = np.zeros(n_g, dtype=np.float64)
            M_y_gp = np.zeros(n_g, dtype=np.float64)
            M_z_gp = np.zeros(n_g, dtype=np.float64)
            dN_dx_list = []
            for k, (xi_g, w_g) in enumerate(zip(xi, w)):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                dN_dx = dN_dξ.copy() * dξ_dx
                dN_dx[:, 3] = dN_dξ[:, 3] * dξ_dx
                d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
                E_lin = self.green_lagrange_strain_operator.strain_linear_part(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                E = E_lin + E_nl
                N_i, M_y_i, M_z_i = self.stress_resultant_operator.section_forces_from_strain(
                    E - self._E_0_voigt, D
                )
                N_gp[k] = N_i
                M_y_gp[k] = M_y_i
                M_z_gp[k] = M_z_i
                dN_dx_list.append(dN_dx[0])
            dN_dx_arr = np.stack(dN_dx_list, axis=0)
            K_sigma = self.geometric_stiffness_operator.assemble_K_sigma(
                N_gp, M_y_gp, M_z_gp, w, dN_dx_arr, detJ
            )
            K_mat = np.zeros((12, 12), dtype=np.float64)
            for xi_g, w_g in zip(xi, w):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                dN_dx = dN_dξ.copy() * dξ_dx
                dN_dx[:, 3] = dN_dξ[:, 3] * dξ_dx
                d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
                B_lin = self.green_lagrange_strain_operator.linearized_strain_displacement(
                    dN_dx[0], d2N_dx2[0]
                )
                B_nl = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                B = B_lin + B_nl
                K_mat += B.T @ D @ B * w_g * detJ
            return K_mat + K_sigma

        U_e = np.asarray(U_e, dtype=np.float64).reshape(_NL_W_DOF)
        U12 = U_e[:_NL_W_STD]
        assert self.warping_material_stiffness_operator is not None
        assert self.warping_strain_displacement_operator is not None
        D7 = self.warping_material_stiffness_operator.assembly_form()
        D6 = D7[:6, :6]
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        dξ_dx = 2.0 / self.L
        d2ξ_dx2 = 4.0 / (self.L ** 2)
        n_g = len(xi)
        N_gp = np.zeros(n_g, dtype=np.float64)
        M_y_gp = np.zeros(n_g, dtype=np.float64)
        M_z_gp = np.zeros(n_g, dtype=np.float64)
        dN_dx_list = []
        for k, (xi_g, w_g) in enumerate(zip(xi, w)):
            _N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            dN_dx = dN_dξ.copy() * dξ_dx
            dN_dx[:, 3] = dN_dξ[:, 3] * dξ_dx
            d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
            E_lin = self.green_lagrange_strain_operator.strain_linear_part(dN_dx[0], d2N_dx2[0], U12)
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx[0], d2N_dx2[0], U12)
            E6 = E_lin + E_nl
            N_i, M_y_i, M_z_i = self.stress_resultant_operator.section_forces_from_strain(
                E6 - self._E_0_voigt, D6
            )
            N_gp[k] = N_i
            M_y_gp[k] = M_y_i
            M_z_gp[k] = M_z_i
            dN_dx_list.append(dN_dx[0])
        dN_dx_arr = np.stack(dN_dx_list, axis=0)
        K_sigma_12 = self.geometric_stiffness_operator.assemble_K_sigma(
            N_gp, M_y_gp, M_z_gp, w, dN_dx_arr, detJ
        )
        K_sigma = np.zeros((_NL_W_DOF, _NL_W_DOF), dtype=np.float64)
        K_sigma[:_NL_W_STD, :_NL_W_STD] = K_sigma_12
        K_mat = np.zeros((_NL_W_DOF, _NL_W_DOF), dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            dN_dx = dN_dξ.copy() * dξ_dx
            dN_dx[:, 3] = dN_dξ[:, 3] * dξ_dx
            d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
            B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx[0], d2N_dx2[0], U12
            )
            B_lin_7 = self.warping_strain_displacement_operator.physical_coordinate_form(
                dN_dξ, d2N_dξ2
            )[0]
            B_tot = B_lin_7 + _B_nl_12_to_14(B_nl_12)
            K_mat += B_tot.T @ D7 @ B_tot * w_g * detJ
        return K_mat + K_sigma

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        """``F_int += B_tot.T @ S * w_g * detJ`` with ``S = D @ (E_lin + E_nl)``."""
        if self._n_dof == 12:
            U_e = np.asarray(U_e, dtype=np.float64).reshape(12)
            D = self.material_stiffness_operator.assembly_form()
            xi, w = self.integration_points
            detJ = self.jacobian_determinant
            dξ_dx = 2.0 / self.L
            d2ξ_dx2 = 4.0 / (self.L ** 2)
            F_int = np.zeros(12, dtype=np.float64)
            for xi_g, w_g in zip(xi, w):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                dN_dx = dN_dξ.copy() * dξ_dx
                d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
                B_lin = self.green_lagrange_strain_operator.linearized_strain_displacement(
                    dN_dx[0], d2N_dx2[0]
                )
                B_nl = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                E_lin = self.green_lagrange_strain_operator.strain_linear_part(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                S = D @ (E_lin + E_nl - self._E_0_voigt)
                F_int += ((B_lin + B_nl).T @ S) * w_g * detJ
            return F_int

        U_e = np.asarray(U_e, dtype=np.float64).reshape(_NL_W_DOF)
        U12 = U_e[:_NL_W_STD]
        assert self.warping_material_stiffness_operator is not None
        assert self.warping_strain_displacement_operator is not None
        D7 = self.warping_material_stiffness_operator.assembly_form()
        E0_7 = np.concatenate([self._E_0_voigt, np.zeros(1, dtype=np.float64)])
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        dξ_dx = 2.0 / self.L
        d2ξ_dx2 = 4.0 / (self.L ** 2)
        F_int = np.zeros(_NL_W_DOF, dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            dN_dx = dN_dξ.copy() * dξ_dx
            d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
            B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx[0], d2N_dx2[0], U12
            )
            B_lin_7 = self.warping_strain_displacement_operator.physical_coordinate_form(
                dN_dξ, d2N_dξ2
            )[0]
            B_tot = B_lin_7 + _B_nl_12_to_14(B_nl_12)
            E_lin = self.green_lagrange_strain_operator.strain_linear_part(dN_dx[0], d2N_dx2[0], U12)
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx[0], d2N_dx2[0], U12)
            E_warp = float(np.dot(B_lin_7[6], U_e))
            E = np.concatenate([E_lin + E_nl, np.array([E_warp], dtype=np.float64)])
            F_int += (B_tot.T @ (D7 @ (E - E0_7))) * w_g * detJ
        return F_int

    def strain_at_gauss_points(self, U_e: np.ndarray) -> List[np.ndarray]:
        """Green–Lagrange ``E_lin + E_nl`` at each Gauss station (first 12 DOFs)."""
        U12 = np.asarray(U_e, dtype=np.float64).ravel()[:12]
        xi, w = self.integration_points
        dξ_dx = 2.0 / self.L
        d2ξ_dx2 = 4.0 / (self.L ** 2)
        result = []
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            dN_dx = dN_dξ.copy() * dξ_dx
            d2N_dx2 = d2N_dξ2.copy() * d2ξ_dx2
            E_lin = self.green_lagrange_strain_operator.strain_linear_part(
                dN_dx[0], d2N_dx2[0], U12
            )
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(
                dN_dx[0], d2N_dx2[0], U12
            )
            result.append(np.asarray(E_lin + E_nl, dtype=np.float64))
        return result
