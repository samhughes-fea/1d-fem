# pre_processing\element_library\truss\truss_3D.py

"""
2-node 3D Truss element: axial, transverse, and torsion.
"""

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.truss.utilities import (
    direction_cosines_and_transverse,
    build_L_matrix_6x12,
    ShapeFunctionOperator,
    StrainDisplacementOperator,
    MaterialStiffnessOperator,
    LoadInterpolationOperator,
)


class TrussElement3D(Element1DBase):
    """
    2-node 3D Truss element carrying axial, transverse, and torsion.

    Stiffness: EA/L (axial), GA/L (transverse), GJ_t/L (torsion).
    Local DOF: 6 (axial×2, transverse×2, torsion×2). K_e = L.T @ K_local @ L (12×12).
    """

    element_type_name = "Truss-3D"

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

        self.node_coords = self.grid_array
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        self._validate_element_properties()
        self._assert_logging_ready()

        axial, transverse = direction_cosines_and_transverse(self.node_coords)
        self._axial = axial
        self._transverse = transverse
        self._L_matrix = build_L_matrix_6x12(axial, transverse)
        self._shape_function_operator = ShapeFunctionOperator(element_length=self.L)
        self._strain_displacement_operator = StrainDisplacementOperator(
            element_length=self.L,
            axial=self._axial,
            transverse=self._transverse,
        )
        self._material_stiffness_operator = MaterialStiffnessOperator(
            youngs_modulus=self.E,
            shear_modulus=self.G,
            cross_section_area=self.A,
            torsion_constant=self.J_t,
            shear_correction_factor=self.kappa,
        )

    def _validate_element_properties(self) -> None:
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4 or self.section_array.size not in (5, 7):
            raise ValueError("Material/section arrays not properly initialised")

    @property
    def A(self) -> float:
        return float(self.section_array[0])

    @property
    def J_t(self) -> float:
        return float(self.section_array[4])

    @property
    def E(self) -> float:
        return float(self.material_array[0])

    @property
    def G(self) -> float:
        return float(self.material_array[1])

    @property
    def kappa(self) -> float:
        """Shear correction for transverse; from section if available."""
        if self.section_array.size >= 7:
            return float(self.section_array[5])
        return 5.0 / 6.0

    def _build_K_local(self) -> np.ndarray:
        """K_local 6×6: block_diag(axial 2×2, transverse 2×2, torsion 2×2)."""
        EA = self.E * self.A
        GA = self.kappa * self.G * self.A
        GJ = self.G * self.J_t
        L = self.L
        if L < 1e-12:
            return np.zeros((6, 6))
        k_axial = (EA / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
        k_trans = (GA / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
        k_torsion = (GJ / L) * np.array([[1, -1], [-1, 1]], dtype=np.float64)
        K_local = np.zeros((6, 6))
        K_local[0:2, 0:2] = k_axial
        K_local[2:4, 2:4] = k_trans
        K_local[4:6, 4:6] = k_torsion
        return K_local

    def element_stiffness_matrix(self):
        """Assemble K_e = L.T @ K_local @ L (12×12)."""
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        K_local = self._build_K_local()
        L_mat = self._L_matrix
        K_e = L_mat.T @ K_local @ L_mat

        if self.logger_operator:
            self.logger_operator.log_matrix("stiffness", K_e, {"name": "Truss Element Stiffness"})
            self.logger_operator.flush("stiffness")

        xi_g = np.array([0.0])
        N, dN_dξ, _ = self._shape_function_operator.natural_coordinate_form(xi_g)
        B = self._strain_displacement_operator.physical_coordinate_form(dN_dξ)[0]
        D = self._material_stiffness_operator.assembly_form()
        gauss_data = [
            StiffnessGaussPointData(
                xi=0.0,
                weight=2.0,
                B_matrix=B,
                D_matrix=D,
                jacobian=self.L / 2,
                shape_functions=N[0].copy(),
                shape_derivatives=dN_dξ[0].copy(),
            )
        ]
        op = self._shape_function_operator
        evaluate_shape_functions = lambda xi: op.natural_coordinate_form(np.asarray(xi))

        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=K_e,
            gauss_data=gauss_data,
            integration_scheme="exact",
            evaluate_shape_functions=evaluate_shape_functions,
            shape_function_N_coefficients=None,
            shape_function_dN_dxi_coefficients=None,
            shape_function_d2N_dxi2_coefficients=None,
        )

    def element_force_vector(self):
        """Element force vector F_e (12,). Point loads mapped to global DOF."""
        from pre_processing.element_library.gauss_point_data import ForceObject, ForceGaussPointData

        self._assert_logging_ready()
        F_e = np.zeros(12, dtype=np.float64)

        if self.point_load_array.size > 0:
            for row in self.point_load_array:
                x_pos = row[0]
                Fx, Fy, Fz = row[3], row[4], row[5]
                Mx = row[6]
                n1, n2 = self.element_array[1], self.element_array[2]
                coords = self.grid_array
                x1, x2 = coords[0, 0], coords[1, 0]
                L = self.L
                if L < 1e-12:
                    continue
                t = (x_pos - x1) / L if abs(x2 - x1) > 1e-12 else 0.5
                t = np.clip(t, 0.0, 1.0)
                N1, N2 = 1 - t, t
                F_e[0] += N1 * Fx
                F_e[1] += N1 * Fy
                F_e[2] += N1 * Fz
                F_e[3] += N1 * Mx
                F_e[6] += N2 * Fx
                F_e[7] += N2 * Fy
                F_e[8] += N2 * Fz
                F_e[9] += N2 * Mx

        point_loads_cache = self.point_load_array.copy() if self.point_load_array.size > 0 else None
        gauss_data = []
        if self.distributed_load_array.size > 0:
            interpolator = LoadInterpolationOperator(
                distributed_loads_array=self.distributed_load_array,
                boundary_mode="clamp",
                interpolation_order="linear",
                n_gauss_points=1,
            )
            x_mid = (self.node_coords[0, 0] + self.node_coords[1, 0]) / 2.0
            q = np.asarray(interpolator.interpolate(x_mid), dtype=np.float64)
            if q.shape == ():
                q = np.atleast_1d(q)
            N_gp, _, _ = self._shape_function_operator.natural_coordinate_form(np.array([0.0]))
            N = N_gp[0]
            F_e += (N @ q) * 2.0 * (self.L / 2)
            gauss_data.append(
                ForceGaussPointData(
                    xi=0.0,
                    weight=2.0,
                    shape_functions=N.copy(),
                    jacobian=self.L / 2,
                    distributed_load=q.copy(),
                )
            )

        if self.logger_operator:
            self.logger_operator.log_matrix("force", F_e.reshape(1, -1), {"name": "Truss Force Vector"})
            self.logger_operator.flush("force")

        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=F_e,
            gauss_data=gauss_data,
            point_loads=point_loads_cache,
        )
