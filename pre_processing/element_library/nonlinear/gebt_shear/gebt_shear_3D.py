# pre_processing/element_library/nonlinear/gebt_shear/gebt_shear_3D.py
"""
2-node 3D GEBT shear beam (Phase 3a): shear-deformable geometrically exact beam theory.

K_T = K_0 + K_σ, F_int = ∫ Bᵀ S. At U_e=0, tangent stiffness equals linear Timoshenko K_e.
Uses same linear Timoshenko B, D and Total Lagrangian operators (Green–Lagrange strain,
geometric stiffness) so limit test passes; full current-config GEBT kinematics can be added later.
"""

import logging
from typing import List, Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.timoshenko.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.shape_function_registry import get_shape_function_operator
from pre_processing.element_library.linear.timoshenko.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.nonlinear.timoshenko.utilities import (
    GreenLagrangeStrainOperator,
    StressResultantOperator,
    GeometricStiffnessOperator,
)

logger = logging.getLogger(__name__)


class GEBTShearBeamElement3D(Element1DBase):
    """
    2-node 3D GEBT shear beam: finite rotations and shear; K_T(U_e), F_int(U_e).

    At U_e=0, tangent stiffness equals linear Timoshenko K_e (same B, D, quadrature).
    Uses Total Lagrangian strain (Green–Lagrange) and geometric stiffness for nonlinear
    response; compatible with nonlinear static runner and build_converged_formulation_cache.
    """

    element_type_name = "GEBTShear-3D"

    def __init__(
        self,
        *,
        element_id: int,
        element_dictionary: dict,
        grid_dictionary: dict,
        section_dictionary: dict,
        material_dictionary: dict,
        point_load_array: np.ndarray,
        distributed_load_array: np.ndarray,
        job_results_dir: str,
        quadrature_order: int = 3,
    ):
        super().__init__(
            element_id=element_id,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
            dof_per_node=6,
        )
        self.quadrature_order = max(int(quadrature_order), 2)
        self.node_coords = self.grid_array
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = float(grid_dictionary["coordinates"][:, 0].min())
        self.x_global_end = float(grid_dictionary["coordinates"][:, 0].max())
        self._validate_element_properties()
        self._assert_logging_ready()

        self.shape_function_operator = get_shape_function_operator(self.__class__.__name__, self.L)
        self.strain_displacement_operator = StrainDisplacementOperator(element_length=self.L)
        self.material_stiffness_operator = MaterialStiffnessOperator(
            youngs_modulus=self.E,
            shear_modulus=self.G,
            cross_section_area=self.A,
            moment_inertia_y=self.I_y,
            moment_inertia_z=self.I_z,
            torsion_constant=self.J_t,
            shear_correction_factor=self.kappa,
        )
        self.green_lagrange_strain_operator = GreenLagrangeStrainOperator(
            element_length=self.L,
            include_shear=True,
        )
        self.stress_resultant_operator = StressResultantOperator()
        self.geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)

        self._K_0: np.ndarray | None = None

    def _validate_element_properties(self) -> None:
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4 or self.section_array.size not in (5, 7):
            raise ValueError("Material/section arrays not properly initialised")
        conn = tuple(self.element_array[1:3])
        logger.debug(
            "Element %s geometry initialised\n"
            "• Connectivity: %s\n• Length: %.4e\n• Start/End X: %.4e / %.4e",
            self.element_id, conn, self.L, self.x_start, self.x_end,
        )

    @property
    def A(self) -> float:
        return float(self.section_array[0])

    @property
    def I_y(self) -> float:
        return float(self.section_array[2])

    @property
    def I_z(self) -> float:
        return float(self.section_array[3])

    @property
    def J_t(self) -> float:
        return float(self.section_array[4])

    @property
    def kappa(self) -> float:
        return float(self.section_array[5]) if self.section_array.size >= 7 else 5.0 / 6.0

    @property
    def E(self) -> float:
        return float(self.material_array[0])

    @property
    def G(self) -> float:
        return float(self.material_array[1])

    @property
    def jacobian_determinant(self) -> float:
        return self.L / 2.0

    @property
    def integration_points(self) -> Tuple[np.ndarray, np.ndarray]:
        return np.polynomial.legendre.leggauss(self.quadrature_order)

    def element_mass_matrix(self):
        """
        Reference-configuration consistent mass (same shape functions as linear Timoshenko) for modal/dynamic.
        """
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        rho = float(self.material_array[3])
        mu = np.zeros(12, dtype=np.float64)
        for i in (0, 1, 2, 6, 7, 8):
            mu[i] = rho * self.A
        for i in (3, 9):
            mu[i] = rho * self.J_t
        for i in (4, 10):
            mu[i] = rho * self.I_y
        for i in (5, 11):
            mu[i] = rho * self.I_z
        M_e = np.zeros((12, 12), dtype=np.float64)
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        for xi_g, w_g in zip(xi, w):
            N, _, _ = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            Ng = N[0]
            for i in range(12):
                for j in range(12):
                    mij = 0.5 * (mu[i] + mu[j])
                    M_e[i, j] += mij * float(np.dot(Ng[i, :], Ng[j, :])) * w_g * detJ
        return MassObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            M_e=M_e,
        )

    def _get_K_0(self) -> np.ndarray:
        """Material stiffness K_0 with same selective integration as linear Timoshenko (1-point shear, bending order)."""
        if self._K_0 is not None:
            return self._K_0
        D = self.material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        axial_order = max(int(self.element_array[3]), 1)
        bending_y_order = max(int(self.element_array[4]), 2)
        bending_z_order = max(int(self.element_array[5]), 2)
        shear_y_order = max(int(self.element_array[6]), 2) if self.element_array[6] > 0 else 2
        shear_z_order = max(int(self.element_array[7]), 2) if self.element_array[7] > 0 else 2
        torsion_order = max(int(self.element_array[8]), 1)
        bending_order = max(bending_y_order, bending_z_order)
        max_order = max(axial_order, bending_order, shear_y_order, shear_z_order, torsion_order)
        xi_full, w_full = np.polynomial.legendre.leggauss(max_order)
        Ke_full = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi_full, w_full):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            Ke_full += B.T @ D @ B * w_g * detJ
        xi_bending, w_bending = np.polynomial.legendre.leggauss(bending_order)
        Ke_bending_block = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi_bending, w_bending):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            B_bending = B[[1, 2], :]
            Ke_bending_block += B_bending.T @ np.diag([D[1, 1], D[2, 2]]) @ B_bending * w_g * detJ
        shear_order = 1
        xi_shear, w_shear = np.polynomial.legendre.leggauss(shear_order)
        Ke_shear_block = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi_shear, w_shear):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            B_shear = B[[3, 4], :]
            Ke_shear_block += B_shear.T @ np.diag([D[3, 3], D[4, 4]]) @ B_shear * w_g * detJ
        Ke_full_bending = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi_full, w_full):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            B_bending = B[[1, 2], :]
            Ke_full_bending += B_bending.T @ np.diag([D[1, 1], D[2, 2]]) @ B_bending * w_g * detJ
        Ke_full_shear = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi_full, w_full):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            B_shear = B[[3, 4], :]
            Ke_full_shear += B_shear.T @ np.diag([D[3, 3], D[4, 4]]) @ B_shear * w_g * detJ
        self._K_0 = Ke_full - Ke_full_bending - Ke_full_shear + Ke_bending_block + Ke_shear_block
        return self._K_0

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        K_0 = self._get_K_0()
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        dξ_dx = 2.0 / self.L
        N_sum, M_y_sum, M_z_sum = 0.0, 0.0, 0.0
        n_g = len(xi)
        dN_dx_list = []
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            dN_dx = dN_dξ.copy() * dξ_dx
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            E_lin = (B @ U_e).ravel()
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx[0], U_e)
            E = E_lin + E_nl
            N_i, M_y_i, M_z_i = self.stress_resultant_operator.section_forces_from_strain(E, D)
            N_sum += N_i * w_g
            M_y_sum += M_y_i * w_g
            M_z_sum += M_z_i * w_g
            dN_dx_list.append(dN_dx[0])
        N_avg = N_sum * detJ / (2.0 * detJ) if n_g > 0 else N_sum
        M_y_avg = M_y_sum * detJ / (2.0 * detJ) if n_g > 0 else M_y_sum
        M_z_avg = M_z_sum * detJ / (2.0 * detJ) if n_g > 0 else M_z_sum
        dN_dx_arr = np.stack(dN_dx_list, axis=0)
        K_sigma = self.geometric_stiffness_operator.assemble_K_sigma(
            N_avg, M_y_avg, M_z_avg, xi, w, dN_dx_arr, detJ
        )
        return K_0 + K_sigma

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        dξ_dx = 2.0 / self.L
        F_int = np.zeros(12, dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            E_lin = (B @ U_e).ravel()
            dN_dx = dN_dξ[0] * dξ_dx
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx, U_e)
            E = E_lin + E_nl
            S = D @ E
            F_int += (B.T @ S) * w_g * detJ
        return F_int

    def strain_at_gauss_points(self, U_e: np.ndarray) -> List[np.ndarray]:
        xi, w = self.integration_points
        dξ_dx = 2.0 / self.L
        result = []
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            E_lin = (B @ U_e).ravel()
            dN_dx = dN_dξ[0] * dξ_dx
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx, U_e)
            E = E_lin + E_nl
            result.append(np.asarray(E, dtype=np.float64))
        return result

    def element_stiffness_matrix(self):
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData
        self._assert_logging_ready()
        K_e = self.tangent_stiffness_matrix(np.zeros(12, dtype=np.float64))
        xi, w = self.integration_points
        D = self.material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        gauss_cache = []
        for g, (xi_g, w_g) in enumerate(zip(xi, w)):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            gauss_cache.append(
                StiffnessGaussPointData(
                    xi=float(xi_g),
                    weight=float(w_g),
                    B_matrix=B.copy(),
                    D_matrix=D.copy(),
                    jacobian=float(detJ),
                    shape_functions=N.copy(),
                    shape_derivatives=dN_dξ.copy(),
                )
            )
        op = self.shape_function_operator
        evaluate_shape_functions = lambda xi_val: op.natural_coordinate_form(np.asarray(xi_val))
        N_coeffs, dN_coeffs, d2N_coeffs = self._build_shape_function_coefficients_b2()
        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=K_e,
            gauss_data=gauss_cache,
            integration_scheme="Gauss-Legendre",
            evaluate_shape_functions=evaluate_shape_functions,
            shape_function_N_coefficients=N_coeffs,
            shape_function_dN_dxi_coefficients=dN_coeffs,
            shape_function_d2N_dxi2_coefficients=d2N_coeffs,
        )

    def _build_shape_function_coefficients_b2(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        c = np.zeros((12, 6, 4), dtype=np.float64)
        dc = np.zeros((12, 6, 4), dtype=np.float64)
        d2c = np.zeros((12, 6, 4), dtype=np.float64)
        c[0, 0, 0], c[0, 0, 1] = 0.5, -0.5
        c[6, 0, 0], c[6, 0, 1] = 0.5, 0.5
        dc[0, 0, 0] = -0.5
        dc[6, 0, 0] = 0.5
        c[1, 1], c[7, 1] = c[0, 0].copy(), c[6, 0].copy()
        dc[1, 1], dc[7, 1] = dc[0, 0].copy(), dc[6, 0].copy()
        c[5, 5], c[11, 5] = c[0, 0].copy(), c[6, 0].copy()
        dc[5, 5], dc[11, 5] = dc[0, 0].copy(), dc[6, 0].copy()
        c[2, 2], c[8, 2] = c[1, 1].copy(), c[7, 1].copy()
        dc[2, 2], dc[8, 2] = dc[1, 1].copy(), dc[7, 1].copy()
        c[4, 4], c[10, 4] = -c[5, 5].copy(), -c[11, 5].copy()
        dc[4, 4], dc[10, 4] = -dc[5, 5].copy(), -dc[11, 5].copy()
        c[3, 3], c[9, 3] = c[0, 0].copy(), c[6, 0].copy()
        dc[3, 3], dc[9, 3] = dc[0, 0].copy(), dc[6, 0].copy()
        return c, dc, d2c

    def element_force_vector(self):
        from pre_processing.element_library.gauss_point_data import ForceObject, ForceGaussPointData
        from pre_processing.element_library.linear.timoshenko.utilities.interpolate_loads import LoadInterpolationOperator
        self._assert_logging_ready()
        Fe = np.zeros(12, dtype=np.float64)
        gauss_cache = []
        x_start = float(self.x_start)
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        if self.distributed_load_array.size > 0:
            interpolator = LoadInterpolationOperator(
                distributed_loads_array=self.distributed_load_array,
                boundary_mode="error",
                interpolation_order="cubic",
                n_gauss_points=self.quadrature_order,
            )
            x_gauss = (xi + 1) * (self.L / 2) + x_start
            q_gauss = interpolator.interpolate(x_gauss)
            N = np.stack([
                self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))[0][0]
                for xi_g in xi
            ])
            Fe_dist = np.einsum("gij,gj,g->i", N, q_gauss, w) * (self.L / 2)
            Fe += Fe_dist
            for g, (xi_g, w_g) in enumerate(zip(xi, w)):
                gauss_cache.append(
                    ForceGaussPointData(
                        xi=float(xi_g),
                        weight=float(w_g),
                        shape_functions=N[g].copy(),
                        jacobian=float(detJ),
                        distributed_load=q_gauss[g].copy() if g < len(q_gauss) else None,
                    )
                )
        if self.point_load_array.size > 0:
            for load in self.point_load_array:
                x_p = float(load[0])
                F_p = load[3:9].astype(np.float64)
                if x_start <= x_p <= x_start + self.L:
                    xi_p = 2 * (x_p - x_start) / self.L - 1
                    N_p = self.shape_function_operator.natural_coordinate_form(np.array([xi_p]))[0][0]
                    Fe_trans = N_p[[0, 1, 2, 6, 7, 8], :3] @ F_p[:3]
                    Fe_rot = N_p[[3, 4, 5, 9, 10, 11], 3:] @ F_p[3:]
                    Fe[[0, 1, 2, 6, 7, 8]] += Fe_trans
                    Fe[[3, 4, 5, 9, 10, 11]] += Fe_rot
        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=gauss_cache,
            point_loads=self.point_load_array.copy() if self.point_load_array.size else None,
        )
