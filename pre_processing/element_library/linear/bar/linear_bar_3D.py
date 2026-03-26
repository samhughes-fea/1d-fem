# pre_processing/element_library/linear/bar/linear_bar_3D.py
"""
2-node 3D bar: axial + torsion only (no bending/shear DOF coupling in ``B``).

**Tensors:** ``U_e`` (12,), ``K_e`` (12,12), ``F_e`` (12,) — six DOF/node for job compatibility.
Per Gauss point ``B`` (2, 12), ``D`` (2, 2) (``EA``, ``GJ_t``); ``eps`` (2,) = [axial strain, twist rate].
``N`` (12, 6) per GP (``linear/bar/utilities/shape_functions.py``). ``detJ = L/2``.

**Weak forms (Gauss, xi in [-1, 1]):** ``K_e += B.T @ D @ B * w_g * detJ``; ``F_dist += w_g * N.T @ q * detJ``;
``F_point = N.T @ P`` at load station; ``M_e`` consistent mass.

**Kinematics / frame:** Local operators use direction cosines along the chord; for 4-DOF axial–torsion global maps,
``linear/bar/utilities/local_frame.build_L_matrix_4x12`` builds ``L`` from the axial unit vector.

**Quadrature:** Order from argument or ``element_array`` (axial, torsion integration columns).

**Public API:** ``element_stiffness_matrix`` → ``ElementObject``; ``element_force_vector`` → ``ForceObject``;
``element_mass_matrix`` → ``MassObject`` where implemented.
"""

import numpy as np
from typing import Optional, Tuple

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.bar.utilities import (
    direction_cosines,
    StrainDisplacementOperator,
    MaterialStiffnessOperator,
    LoadInterpolationOperator,
)
from pre_processing.element_library.shape_function_registry import get_shape_function_operator


class LinearBarElement3D(Element1DBase):
    """
    2 nodes, 6 DOF/node; only axial and torsion strains are non-zero in ``B``.

    Notes
    -----
    Tensor shapes, ``detJ``, and ``build_L_matrix_4x12`` for 4-DOF global maps: module docstring.
    """

    element_type_name = "Bar-3D"

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
        quadrature_order: Optional[int] = None,
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
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self._validate_element_properties()
        self._assert_logging_ready()

        # Quadrature order: from integration_orders when not provided (axial, torsion)
        if quadrature_order is not None:
            self.quadrature_order = quadrature_order
        else:
            axial_order = int(self.element_array[3])
            torsion_order = int(self.element_array[8])
            self.quadrature_order = max(axial_order, torsion_order, 1)

        self._axial = direction_cosines(self.node_coords)
        self.shape_function_operator = get_shape_function_operator(self.__class__.__name__, self.L)
        self.strain_displacement_operator = StrainDisplacementOperator(
            element_length=self.L, axial=self._axial
        )
        self.material_stiffness_operator = MaterialStiffnessOperator(
            youngs_modulus=self.E,
            shear_modulus=self.G,
            cross_section_area=self.A,
            torsion_constant=self.J_t,
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
    def rho(self) -> float:
        """Mass density for consistent mass matrix."""
        return float(self.material_array[3])

    @property
    def integration_points(self) -> Tuple[np.ndarray, np.ndarray]:
        """Gauss-Legendre points and weights for stiffness/force quadrature."""
        return np.polynomial.legendre.leggauss(self.quadrature_order)

    @property
    def jacobian_determinant(self) -> float:
        return self.L / 2

    def _build_shape_function_coefficients_b2(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """B2 coefficient arrays (12, 6, 4) for N, dN/dξ, d²N/dξ² in monomial basis ξ^0..ξ^3. Bar: linear only."""
        c = np.zeros((12, 6, 4), dtype=np.float64)
        dc = np.zeros((12, 6, 4), dtype=np.float64)
        d2c = np.zeros((12, 6, 4), dtype=np.float64)
        # Axial: N1 = 0.5(1-ξ), N2 = 0.5(1+ξ) -> comp 0, DOF 0 and 6
        c[0, 0, 0], c[0, 0, 1] = 0.5, -0.5
        c[6, 0, 0], c[6, 0, 1] = 0.5, 0.5
        dc[0, 0, 0] = -0.5
        dc[6, 0, 0] = 0.5
        # Torsion: comp 3, DOF 3 and 9
        c[3, 3, 0], c[3, 3, 1] = 0.5, -0.5
        c[9, 3, 0], c[9, 3, 1] = 0.5, 0.5
        dc[3, 3, 0] = -0.5
        dc[9, 3, 0] = 0.5
        return c, dc, d2c

    def _is_point_in_element(self, x: float) -> bool:
        """True if x is within element bounds [x_start, x_end] with small tolerance."""
        tol = 1e-12 * max(self.L, 1.0)
        return (self.x_start - tol <= x <= self.x_end + tol)

    def _compute_point_load_contribution(self) -> np.ndarray:
        """Point load contribution using shape-function operator; all 6 load components."""
        Fe_point = np.zeros(12, dtype=np.float64)
        for load in self.point_load_array:
            x_p = float(load[0])
            F_p = load[3:9].astype(np.float64)
            if not self._is_point_in_element(x_p):
                continue
            xi_p = 2 * (x_p - self.x_start) / self.L - 1
            N_p = self.shape_function_operator.natural_coordinate_form(np.array([xi_p]))[0][0]
            Fe_trans = np.einsum("ij,j->i", N_p[[0, 1, 2, 6, 7, 8], :3], F_p[:3])
            Fe_rot = np.einsum("ij,j->i", N_p[[3, 4, 5, 9, 10, 11], 3:], F_p[3:])
            Fe_point[[0, 1, 2, 6, 7, 8]] += Fe_trans
            Fe_point[[3, 4, 5, 9, 10, 11]] += Fe_rot
        return Fe_point

    def _compute_distributed_load_contribution(self) -> Tuple[np.ndarray, list]:
        """Distributed load: ``F_dist += sum_g w_g * N.T @ q * detJ`` at Gauss ``xi_g``."""
        from pre_processing.element_library.gauss_point_data import ForceGaussPointData

        xi_gauss, weights = self.integration_points
        x_gauss = (xi_gauss + 1) * (self.L / 2) + self.x_start
        Fe_dist = np.zeros(12, dtype=np.float64)
        gauss_cache = []
        detJ = self.jacobian_determinant
        interpolator = LoadInterpolationOperator(
            distributed_loads_array=self.distributed_load_array,
            boundary_mode="clamp",
            interpolation_order="linear",
            n_gauss_points=self.quadrature_order,
        )
        q_gauss = interpolator.interpolate(x_gauss)
        if q_gauss.ndim == 1:
            q_gauss = q_gauss.reshape(1, -1)
        N_stack = np.stack([
            self.shape_function_operator.natural_coordinate_form(np.array([xi]))[0][0]
            for xi in xi_gauss
        ])
        Fe_dist = np.einsum("gij,gj,g->i", N_stack, q_gauss, weights) * detJ
        for g, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
            gauss_cache.append(ForceGaussPointData(
                xi=float(xi_g),
                weight=float(w_g),
                shape_functions=N_stack[g].copy(),
                jacobian=float(detJ),
                distributed_load=q_gauss[g].copy(),
            ))
        return Fe_dist, gauss_cache

    def element_stiffness_matrix(self):
        """Stiffness ``K_e += sum_g B.T @ D @ B * w_g * detJ`` (Gauss-Legendre on ``xi``)."""
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        Ke = np.zeros((12, 12), dtype=np.float64)
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        gauss_cache = []

        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, _ = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ)[0]
            detJ = self.jacobian_determinant
            Ke += B.T @ D @ B * w_g * detJ
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g),
                weight=float(w_g),
                B_matrix=B.copy(),
                D_matrix=D.copy(),
                jacobian=float(detJ),
                shape_functions=N[0].copy(),
                shape_derivatives=dN_dξ[0].copy(),
            ))

        if self.logger_operator:
            self.logger_operator.log_matrix("stiffness", Ke, {"name": "Bar Element Stiffness"})
            self.logger_operator.flush("stiffness")

        op = self.shape_function_operator
        evaluate_shape_functions = lambda xi: op.natural_coordinate_form(np.asarray(xi))
        N_coeffs, dN_coeffs, d2N_coeffs = self._build_shape_function_coefficients_b2()

        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=Ke,
            gauss_data=gauss_cache,
            integration_scheme="Gauss-Legendre",
            evaluate_shape_functions=evaluate_shape_functions,
            shape_function_N_coefficients=N_coeffs,
            shape_function_dN_dxi_coefficients=dN_coeffs,
            shape_function_d2N_dxi2_coefficients=d2N_coeffs,
        )

    def element_force_vector(self):
        """F_e = point loads (via shape-function operator) + distributed loads (multi-GP quadrature)."""
        from pre_processing.element_library.gauss_point_data import ForceObject

        self._assert_logging_ready()
        F_e = np.zeros(12, dtype=np.float64)
        gauss_cache = []

        if self.distributed_load_array.size > 0:
            Fe_dist, gauss_cache = self._compute_distributed_load_contribution()
            F_e += Fe_dist

        if self.point_load_array.size > 0:
            F_e += self._compute_point_load_contribution()
        point_loads_cache = self.point_load_array.copy() if self.point_load_array.size > 0 else None

        if self.logger_operator:
            self.logger_operator.log_matrix("force", F_e.reshape(1, -1), {"name": "Bar Force Vector"})
            self.logger_operator.flush("force")

        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=F_e,
            gauss_data=gauss_cache,
            point_loads=point_loads_cache,
        )

    def element_mass_matrix(self):
        """Consistent mass matrix M_e = rho * A * (L/2) * integral N @ N.T d xi (same DOF order as K_e)."""
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        M_e = np.zeros((12, 12), dtype=np.float64)
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        for xi_g, w_g in zip(xi, w):
            N, _, _ = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            N_g = N[0]
            M_e += self.rho * self.A * detJ * w_g * (N_g @ N_g.T)
        return MassObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            M_e=M_e,
        )
