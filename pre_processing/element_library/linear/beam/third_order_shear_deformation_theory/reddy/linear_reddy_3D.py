# pre_processing/element_library/linear/reddy/linear_reddy_3D.py
"""
2-node 3D Reddy beam (variationally consistent third-order shear).

**Tensors:** Same sizes as Levinson — ``U_e`` (12,), ``K_e`` (12,12), ``F_e`` (12,), ``B`` (6,12), ``D`` (6,6), ``eps``, ``S`` —
per ``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md``. ``D`` matches Levinson (``G*A`` shear, no ``kappa``).

**Weak forms (Gauss, xi in [-1, 1]):** ``K_e += B.T @ D @ B * w_g * detJ`` with selective bending/shear rules as Levinson;
``F_dist += w_g * N.T @ q * detJ``; ``F_point = N.T @ P``; ``M_e`` per ``FORMULATION_DOCSTRING_STANDARDS.md``. ``detJ = L/2``.

**Kinematics:** Same Voigt strain layout as Levinson; ``B`` includes section parameter ``alpha`` (``StrainDisplacementOperator``),
from ``section_array[6]`` or ``I_z/A``. Operators re-exported from Levinson where noted in ``utilities/``.

**Quadrature:** Selective integration in ``element_stiffness_matrix`` matches Levinson-style bending/shear split.

**Public API:** ``element_stiffness_matrix`` → ``ElementObject``; ``element_force_vector`` → ``ForceObject``;
``element_mass_matrix`` → ``MassObject``.
"""

import numpy as np
from typing import Tuple

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.shape_function_registry import get_shape_function_operator
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy.utilities.interpolate_loads import LoadInterpolationOperator

import logging

logger = logging.getLogger(__name__)


class LinearReddyBeamElement3D(Element1DBase):
    """
    **Identity:** 2 nodes, 6 DOF/node, 12 global DOFs (standard beam ordering).

    **Tensors / kinematics:** **Same Voigt row meaning as Levinson**; Reddy–Levinson distinction is in the
    **strain–displacement** definitions (α-term). **Selective integration** in ``element_stiffness_matrix`` matches
    the Levinson-style bending/shear split (see implementation).
    """

    element_type_name = "Reddy-3D"

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

        if quadrature_order is None or quadrature_order == 3:
            axial_order = self.element_array[3]
            bending_y_order = self.element_array[4]
            bending_z_order = self.element_array[5]
            shear_y_order = self.element_array[6]
            shear_z_order = self.element_array[7]
            torsion_order = self.element_array[8]
            if shear_y_order == 0:
                shear_y_order = 3
            if shear_z_order == 0:
                shear_z_order = 3
            max_order = max(
                axial_order, bending_y_order, bending_z_order,
                shear_y_order, shear_z_order, torsion_order,
            )
            self.quadrature_order = max(max_order, 3)
        else:
            self.quadrature_order = quadrature_order

        self.node_coords = self.grid_array
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = float(grid_dictionary["coordinates"][:, 0].min())
        self.x_global_end = float(grid_dictionary["coordinates"][:, 0].max())

        self._validate_element_properties()
        self._assert_logging_ready()

        alpha_coeff = float(self.section_array[6]) if self.section_array.size >= 7 else (
            self.I_z / self.A if self.A > 0 else 0.0
        )

        self.shape_function_operator = get_shape_function_operator(self.__class__.__name__, self.L)
        self.strain_displacement_operator = StrainDisplacementOperator(
            element_length=self.L,
            alpha_coefficient=alpha_coeff,
        )
        self.material_stiffness_operator = MaterialStiffnessOperator(
            youngs_modulus=self.E,
            shear_modulus=self.G,
            cross_section_area=self.A,
            moment_inertia_y=self.I_y,
            moment_inertia_z=self.I_z,
            torsion_constant=self.J_t,
        )

    def _validate_element_properties(self) -> None:
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4 or self.section_array.size not in (5, 7):
            raise ValueError("Material/section arrays not properly initialised")

    @property
    def A(self) -> float:
        return self.section_array[0]

    @property
    def I_y(self) -> float:
        return self.section_array[2]

    @property
    def I_z(self) -> float:
        return self.section_array[3]

    @property
    def J_t(self) -> float:
        return self.section_array[4]

    @property
    def E(self) -> float:
        return self.material_array[0]

    @property
    def G(self) -> float:
        return self.material_array[1]

    @property
    def jacobian_determinant(self) -> float:
        return self.L / 2

    @property
    def integration_points(self) -> Tuple[np.ndarray, np.ndarray]:
        return np.polynomial.legendre.leggauss(self.quadrature_order)

    def element_stiffness_matrix(self):
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        D = self.material_stiffness_operator.assembly_form()
        axial_order = max(self.element_array[3], 1)
        bending_y_order = max(self.element_array[4], 2)
        bending_z_order = max(self.element_array[5], 2)
        shear_y_order = max(self.element_array[6], 3) if self.element_array[6] > 0 else 3
        shear_z_order = max(self.element_array[7], 3) if self.element_array[7] > 0 else 3
        torsion_order = max(self.element_array[8], 1)
        max_order = max(
            axial_order, bending_y_order, bending_z_order,
            shear_y_order, shear_z_order, torsion_order,
        )
        gauss_cache = []
        xi_full, w_full = np.polynomial.legendre.leggauss(max_order)
        Ke_full = np.zeros((12, 12), dtype=np.float64)
        detJ = self.jacobian_determinant

        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N, dN_dxi, d2N_dxi2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
            Ke_full += B.T @ D @ B * w_g * detJ
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g),
                weight=float(w_g),
                B_matrix=B.copy(),
                D_matrix=D.copy(),
                jacobian=float(detJ),
                shape_functions=N.copy(),
                shape_derivatives=dN_dxi.copy(),
            ))

        bending_order = max(bending_y_order, bending_z_order)
        xi_bending, w_bending = np.polynomial.legendre.leggauss(bending_order)
        Ke_bending_block = np.zeros((12, 12), dtype=np.float64)
        for g, (xi_g, w_g) in enumerate(zip(xi_bending, w_bending)):
            N, dN_dxi, d2N_dxi2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
            B_bending = B[[1, 2], :]
            D_bending = D[np.ix_([1, 2], [1, 2])]
            Ke_bending_block += B_bending.T @ D_bending @ B_bending * w_g * detJ

        shear_order = max(shear_y_order, shear_z_order)
        xi_shear, w_shear = np.polynomial.legendre.leggauss(shear_order)
        Ke_shear_block = np.zeros((12, 12), dtype=np.float64)
        for g, (xi_g, w_g) in enumerate(zip(xi_shear, w_shear)):
            N, dN_dxi, d2N_dxi2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
            B_shear = B[[3, 4], :]
            D_shear_diag = np.array([D[3, 3], D[4, 4]])
            Ke_shear_block += B_shear.T @ np.diag(D_shear_diag) @ B_shear * w_g * detJ

        Ke_full_bending = np.zeros((12, 12), dtype=np.float64)
        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N, dN_dxi, d2N_dxi2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
            Ke_full_bending += B[[1, 2], :].T @ D[np.ix_([1, 2], [1, 2])] @ B[[1, 2], :] * w_g * detJ
        Ke_full_shear = np.zeros((12, 12), dtype=np.float64)
        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N, dN_dxi, d2N_dxi2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
            Ke_full_shear += B[[3, 4], :].T @ np.diag([D[3, 3], D[4, 4]]) @ B[[3, 4], :] * w_g * detJ

        Ke = Ke_full - Ke_full_bending - Ke_full_shear + Ke_bending_block + Ke_shear_block
        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=Ke,
            gauss_data=gauss_cache,
            integration_scheme="Gauss-Legendre (Selective)",
        )

    def element_force_vector(self):
        from pre_processing.element_library.gauss_point_data import ForceObject, ForceGaussPointData

        self._assert_logging_ready()
        Fe = np.zeros(12, dtype=np.float64)
        gauss_cache = []
        detJ = self.jacobian_determinant

        if self.distributed_load_array.size > 0:
            xi_gauss, weights = self.integration_points
            x_gauss = (xi_gauss + 1) * (self.L / 2) + self.x_start
            interpolator = LoadInterpolationOperator(
                distributed_loads_array=self.distributed_load_array,
                boundary_mode="error",
                interpolation_order="cubic",
                n_gauss_points=self.quadrature_order,
            )
            q_gauss = interpolator.interpolate(x_gauss)
            N = np.stack([
                self.shape_function_operator.natural_coordinate_form(np.array([xi]))[0][0]
                for xi in xi_gauss
            ])
            Fe = np.einsum("gij,gj,g->i", N, q_gauss, weights) * (self.L / 2)
            for g, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
                gauss_cache.append(ForceGaussPointData(
                    xi=float(xi_g),
                    weight=float(w_g),
                    shape_functions=N[g].copy(),
                    jacobian=float(detJ),
                    distributed_load=q_gauss[g].copy() if g < len(q_gauss) else None,
                ))

        if self.point_load_array.size > 0:
            for load in self.point_load_array:
                x_p = float(load[0])
                F_p = load[3:9].astype(np.float64)
                if not self._is_point_in_element(x_p):
                    continue
                xi_p = 2 * (x_p - self.x_start) / self.L - 1
                N_p = self.shape_function_operator.natural_coordinate_form(np.array([xi_p]))[0][0]
                Fe[[0, 1, 2, 6, 7, 8]] += np.einsum("ij,j->i", N_p[[0, 1, 2, 6, 7, 8], :3], F_p[:3])
                Fe[[3, 4, 5, 9, 10, 11]] += np.einsum("ij,j->i", N_p[[3, 4, 5, 9, 10, 11], 3:], F_p[3:])

        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=gauss_cache,
            point_loads=self.point_load_array.copy() if self.point_load_array.size > 0 else None,
        )

    def _is_point_in_element(self, x: float) -> bool:
        tol = 1e-12 * self.L
        if np.isclose(self.x_end, self.x_global_end):
            return self.x_start - tol <= x <= self.x_end + tol
        return self.x_start - tol <= x < self.x_end

    def element_mass_matrix(self):
        """
        Consistent mass using Reddy/Levinson kinematic shape functions (same DOF layout as Timoshenko beam):
        translations ρA, torsion ρJ_t, bending rotations ρI_y / ρI_z.
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
