# pre_processing/element_library/linear/curved_beam/linear_curved_timoshenko_3D.py
"""2-node 3D curved Timoshenko beam with constant initial curvature κ0. K_e (12, 12), F_e (12,). When κ0=0 reduces to straight Timoshenko."""

import numpy as np
from typing import Tuple

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.curved_beam.utilities.B_matrix import CurvedStrainDisplacementOperator
from pre_processing.element_library.linear.curved_beam.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.curved_beam.utilities.interpolate_loads import LoadInterpolationOperator
from pre_processing.element_library.shape_function_registry import get_shape_function_operator
import logging

logger = logging.getLogger(__name__)


class LinearCurvedTimoshenkoBeamElement3D(Element1DBase):
    """
    2-node 3D curved Timoshenko beam. Constant curvature κ0 (1/R) in x-y plane;
    strain coupling: ε_s = du_x/ds - κ0*u_y, γ_xy = du_y/ds + κ0*u_x - θ_z.
    Reads curvature from element_dictionary["curvature"]; κ0=0 gives straight Timoshenko.
    """

    element_type_name = "CurvedTimoshenko-3D"

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
        idx = int(np.where(element_dictionary["ids"] == element_id)[0][0])
        self.curvature = float(
            element_dictionary.get("curvature", np.zeros(len(element_dictionary["ids"]), dtype=np.float64))[idx]
        )

        if quadrature_order is None or quadrature_order == 3:
            axial_order = self.element_array[3]
            bending_y_order = self.element_array[4]
            bending_z_order = self.element_array[5]
            shear_y_order = max(self.element_array[6], 2) if self.element_array[6] > 0 else 2
            shear_z_order = max(self.element_array[7], 2) if self.element_array[7] > 0 else 2
            torsion_order = self.element_array[8]
            max_order = max(
                axial_order, bending_y_order, bending_z_order,
                shear_y_order, shear_z_order, torsion_order,
            )
            self.quadrature_order = max(max_order, 2)
        else:
            self.quadrature_order = quadrature_order

        self.node_coords = self.grid_array
        self.L = np.linalg.norm(self.node_coords[1] - self.node_coords[0])
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = grid_dictionary["coordinates"][:, 0].min()
        self.x_global_end = grid_dictionary["coordinates"][:, 0].max()

        self._validate_element_properties()
        self._assert_logging_ready()

        self.shape_function_operator = get_shape_function_operator("LinearTimoshenkoBeamElement3D", self.L)
        self.strain_displacement_operator = CurvedStrainDisplacementOperator(
            element_length=self.L, curvature=self.curvature
        )
        self.material_stiffness_operator = MaterialStiffnessOperator(
            youngs_modulus=self.E,
            shear_modulus=self.G,
            cross_section_area=self.A,
            moment_inertia_y=self.I_y,
            moment_inertia_z=self.I_z,
            torsion_constant=self.J_t,
            shear_correction_factor=self.kappa,
            y_sc=self.y_sc,
            z_sc=self.z_sc,
        )

    def _validate_element_properties(self) -> None:
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4 or self.section_array.size not in (5, 7, 9, 10):
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
    def kappa(self) -> float:
        if self.section_array.size >= 7:
            return float(self.section_array[5])
        return 5.0 / 6.0

    @property
    def y_sc(self) -> float:
        if self.section_array.size >= 9:
            return float(self.section_array[7])
        return 0.0

    @property
    def z_sc(self) -> float:
        if self.section_array.size >= 9:
            return float(self.section_array[8])
        return 0.0

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
        Ke = np.zeros((12, 12), dtype=np.float64)
        D = self.material_stiffness_operator.assembly_form()
        xi_full, w_full = self.integration_points
        gauss_cache = []

        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            detJ = self.jacobian_determinant
            Ke += B.T @ D @ B * w_g * detJ
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g),
                weight=float(w_g),
                B_matrix=B.copy(),
                D_matrix=D.copy(),
                jacobian=float(detJ),
                shape_functions=N.copy(),
                shape_derivatives=dN_dξ.copy(),
            ))

        if self.logger_operator:
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} Curved Timoshenko K_e (κ0={self.curvature}) ===")
            self.logger_operator.log_matrix("stiffness", Ke, {"name": "Element Stiffness Matrix"})
            self.logger_operator.flush("stiffness")

        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=Ke,
            gauss_data=gauss_cache,
            integration_scheme="Gauss-Legendre",
        )

    def element_force_vector(self):
        from pre_processing.element_library.gauss_point_data import ForceObject, ForceGaussPointData

        self._assert_logging_ready()
        Fe = np.zeros(12, dtype=np.float64)
        gauss_cache = []
        detJ = self.jacobian_determinant

        if self.distributed_load_array.size > 0:
            interpolator = LoadInterpolationOperator(
                distributed_loads_array=self.distributed_load_array,
                boundary_mode="error",
                interpolation_order="cubic",
                n_gauss_points=self.quadrature_order,
            )
            xi_gauss, weights = self.integration_points
            x_gauss = (xi_gauss + 1) * (self.L / 2) + self.x_start
            q_gauss = interpolator.interpolate(x_gauss)
            if q_gauss.ndim == 1:
                q_gauss = q_gauss.reshape(1, -1)
            N_stack = np.stack([
                self.shape_function_operator.natural_coordinate_form(np.array([xi]))[0][0]
                for xi in xi_gauss
            ])
            Fe = np.einsum("gij,gj,g->i", N_stack, q_gauss, weights) * (self.L / 2)
            for gp, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
                gauss_cache.append(ForceGaussPointData(
                    xi=float(xi_g),
                    weight=float(w_g),
                    shape_functions=N_stack[gp].copy(),
                    jacobian=float(detJ),
                    distributed_load=q_gauss[gp].copy() if gp < len(q_gauss) else None,
                ))

        if self.point_load_array.size > 0:
            Fe += self._compute_point_load_contribution()

        if self.logger_operator:
            self.logger_operator.log_text("force", f"\n=== Element {self.element_id} Force Vector ===")
            self.logger_operator.log_matrix("force", Fe.reshape(1, -1), {"name": "Final Force Vector"})
            self.logger_operator.flush("force")

        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=gauss_cache,
            point_loads=self.point_load_array.copy() if self.point_load_array.size > 0 else None,
        )

    def _compute_point_load_contribution(self) -> np.ndarray:
        Fe = np.zeros(12, dtype=np.float64)
        for load in self.point_load_array:
            x_p = float(load[0])
            F_p = load[3:9].astype(np.float64)
            if x_p < self.x_start - 1e-9 or x_p > self.x_end + 1e-9:
                continue
            xi_p = 2.0 * (x_p - self.x_start) / self.L - 1.0
            N_p = self.shape_function_operator.natural_coordinate_form(np.array([xi_p]))[0][0]
            Fe[[0, 1, 2, 6, 7, 8]] += np.einsum("ij,j->i", N_p[[0, 1, 2, 6, 7, 8], :3], F_p[:3])
            Fe[[3, 4, 5, 9, 10, 11]] += np.einsum("ij,j->i", N_p[[3, 4, 5, 9, 10, 11], 3:], F_p[3:])
        return Fe
