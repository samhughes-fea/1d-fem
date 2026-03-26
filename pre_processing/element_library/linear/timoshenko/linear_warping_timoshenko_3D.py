# pre_processing/element_library/linear/timoshenko/linear_warping_timoshenko_3D.py
"""
2-node 3D Timoshenko beam with 7 DOF/node (six + warping chi) for Vlasov non-uniform torsion.

**Tensors:** ``U_e`` (14,) — standard Timoshenko layout for DOFs 0–11, chi at 12–13. Per Gauss point ``B`` (7, 14) —
first six rows from linear Timoshenko ``B`` (6, 12) on columns 0–11; row 6 is ``phi_x_prime`` (same as warping EB).
``D`` (7, 7) — 6x6 Timoshenko ``D`` (``kappa*G*A`` shear) plus ``D[6,6] = E*Gamma``. ``eps`` (7,), ``S = D @ eps``;
``detJ = L/2``.

**Weak forms (Gauss, xi in [-1, 1]):** ``K_e += B.T @ D @ B * w_g * detJ`` with ``D`` (7,7); ``F_dist``, ``F_point``, ``M_e`` as for warping EB.

**Kinematics / shapes:** Registry ``LinearTimoshenkoBeamElement3D`` → ``N`` (12, 6); warping extension in ``_N_14x6_at_xi``.

**Quadrature:** Same selective shear logic as linear Timoshenko where applicable (see ``element_stiffness_matrix``).

**Mass:** ``M_e`` (14, 14) consistent on ``N`` (14, 6) with ``rho*Gamma`` on warping DOFs.

**Public API:** ``element_stiffness_matrix`` → ``ElementObject``; ``element_force_vector`` → ``ForceObject``;
``element_mass_matrix`` → ``MassObject``.
"""

import numpy as np
from typing import Tuple

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.timoshenko.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.timoshenko.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.shape_function_registry import get_shape_function_operator
from pre_processing.element_library.linear.timoshenko.utilities.interpolate_loads import LoadInterpolationOperator
import logging

logger = logging.getLogger(__name__)

# Number of standard DOF and warping DOF
N_STANDARD_DOF = 12
N_WARPING_DOF = 2
N_DOF = N_STANDARD_DOF + N_WARPING_DOF  # 14
N_STRAIN = 7  # axial, κ_y, κ_z, γ_xy, γ_xz, φ_x, φ_x′


class LinearWarpingTimoshenkoBeamElement3D(Element1DBase):
    """
    2 nodes, 7 DOF/node; same chord-frame convention as linear Timoshenko.

    Notes
    -----
    Assembly: ``K_e += B.T @ D @ B * w_g * detJ`` with ``B`` (7,14), ``D`` (7,7).
    If ``quadrature_order`` is ``None`` or ``3``, derive from ``element_array`` (shear columns at least 2); else use the given integer.
    """

    element_type_name = "WarpingTimoshenko-3D"

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
            dof_per_node=7,
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
                shear_y_order, shear_z_order, torsion_order
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
        self.strain_displacement_operator = StrainDisplacementOperator(element_length=self.L)
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
        if self.material_array.size != 4:
            raise ValueError("Material array not properly initialised")
        if self.section_array.size not in (5, 7, 9, 10):
            raise ValueError("Section array must have 5, 7, 9 or 10 entries")

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
    def Gamma(self) -> float:
        """Warping constant Γ (general section integration); 0 if not in section_array."""
        if self.section_array.size >= 10:
            return float(self.section_array[9])
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

    def _B_warping_row(self, dξ_dx: float) -> np.ndarray:
        """Row 6 of B: φ_x′ = dθ_x/dx + dχ/dx. Linear shape for θ_x and warping: dN/dx = ±1/L."""
        row = np.zeros(N_DOF, dtype=np.float64)
        # θ_x at DOF 3 (node1), 9 (node2): linear N1=(1-ξ)/2, N2=(1+ξ)/2 => dN1/dx=-1/L, dN2/dx=1/L
        row[3] = -1.0 / self.L
        row[9] = 1.0 / self.L
        # Warping DOF 12, 13: same linear
        row[12] = -1.0 / self.L
        row[13] = 1.0 / self.L
        return row

    def _B_matrix_7x14(self, dN_dξ: np.ndarray, d2N_dξ2: np.ndarray, N: np.ndarray) -> np.ndarray:
        """Build B (7, 14) from Timoshenko B (6, 12) plus warping row."""
        n_gauss = dN_dξ.shape[0]
        B_6x12 = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)
        B = np.zeros((n_gauss, N_STRAIN, N_DOF), dtype=np.float64)
        B[:, :6, :12] = B_6x12
        for g in range(n_gauss):
            B[g, 6, :] = self._B_warping_row(self.strain_displacement_operator.dξ_dx)
        return B

    def _D_matrix_7x7(self) -> np.ndarray:
        """D (7,7): first 6×6 from Timoshenko, D[6,6] = E·Γ."""
        D6 = self.material_stiffness_operator.assembly_form()
        D = np.zeros((N_STRAIN, N_STRAIN), dtype=np.float64)
        D[:6, :6] = D6
        D[6, 6] = self.E * self.Gamma
        return D

    def element_stiffness_matrix(self):
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        Ke = np.zeros((N_DOF, N_DOF), dtype=np.float64)
        D = self._D_matrix_7x7()
        xi_full, w_full = np.polynomial.legendre.leggauss(self.quadrature_order)
        gauss_cache = []

        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N_12, dN_dξ_12, d2N_dξ2_12 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            N = N_12
            dN_dξ = dN_dξ_12
            d2N_dξ2 = d2N_dξ2_12
            B = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
            detJ = self.jacobian_determinant
            Ke += B.T @ D @ B * w_g * detJ
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g),
                weight=float(w_g),
                B_matrix=B.copy(),
                D_matrix=D.copy(),
                jacobian=float(detJ),
                shape_functions=None,
                shape_derivatives=None,
            ))

        if self.logger_operator:
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} Warping Timoshenko K_e (14,14) ===")
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
        from pre_processing.element_library.gauss_point_data import ForceObject

        self._assert_logging_ready()
        Fe = np.zeros(N_DOF, dtype=np.float64)

        if self.distributed_load_array.size > 0:
            Fe_dist, _ = self._compute_distributed_load_contribution()
            Fe[:N_STANDARD_DOF] += Fe_dist

        if self.point_load_array.size > 0:
            Fe[:N_STANDARD_DOF] += self._compute_point_load_contribution()

        if self.logger_operator:
            self.logger_operator.log_text("force", f"\n=== Element {self.element_id} Force Vector (14,) ===")
            self.logger_operator.log_matrix("force", Fe.reshape(1, -1), {"name": "Final Force Vector"})
            self.logger_operator.flush("force")

        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=[],
            point_loads=self.point_load_array.copy() if self.point_load_array.size > 0 else None,
        )

    def _N_14x6_at_xi(self, xi_g: float, N12: np.ndarray) -> np.ndarray:
        """Extend standard (12, 6) N to (14, 6); warping DOFs use linear Lagrange on θ_x component."""
        Nf = np.zeros((N_DOF, 6), dtype=np.float64)
        Nf[:N_STANDARD_DOF, :] = N12
        xi = float(xi_g)
        Nf[12, 3] = 0.5 * (1.0 - xi)
        Nf[13, 3] = 0.5 * (1.0 + xi)
        return Nf

    def element_mass_matrix(self):
        """
        Consistent mass: first 12 DOFs as straight Timoshenko; warping DOFs 12–13 with ρ·Γ
        and linear shape (θ_x component slot).
        """
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        rho = float(self.material_array[3])
        mu = np.zeros(N_DOF, dtype=np.float64)
        for i in (0, 1, 2, 6, 7, 8):
            mu[i] = rho * self.A
        for i in (3, 9):
            mu[i] = rho * self.J_t
        for i in (4, 10):
            mu[i] = rho * self.I_y
        for i in (5, 11):
            mu[i] = rho * self.I_z
        mu[12] = rho * self.Gamma
        mu[13] = rho * self.Gamma
        M_e = np.zeros((N_DOF, N_DOF), dtype=np.float64)
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        for xi_g, w_g in zip(xi, w):
            N12, _, _ = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            Ng = self._N_14x6_at_xi(xi_g, N12[0])
            for i in range(N_DOF):
                for j in range(N_DOF):
                    mij = 0.5 * (mu[i] + mu[j])
                    M_e[i, j] += mij * float(np.dot(Ng[i, :], Ng[j, :])) * w_g * detJ
        return MassObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            M_e=M_e,
        )

    def _compute_distributed_load_contribution(self):
        from pre_processing.element_library.gauss_point_data import ForceGaussPointData

        xi_gauss, weights = self.integration_points
        x_gauss = (xi_gauss + 1) * (self.L / 2) + self.x_start
        Fe_dist = np.zeros(N_STANDARD_DOF, dtype=np.float64)
        gauss_cache = []
        detJ = self.jacobian_determinant
        if self.distributed_load_array.size == 0:
            return Fe_dist, gauss_cache
        interpolator = LoadInterpolationOperator(
            distributed_loads_array=self.distributed_load_array,
            boundary_mode="error",
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
        Fe_dist = np.einsum("gij,gj,g->i", N_stack, q_gauss, weights) * (self.L / 2)
        for g, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
            gauss_cache.append(ForceGaussPointData(
                xi=float(xi_g),
                weight=float(w_g),
                shape_functions=N_stack[g].copy(),
                jacobian=float(detJ),
                distributed_load=q_gauss[g].copy() if g < len(q_gauss) else None,
            ))
        return Fe_dist, gauss_cache

    def _compute_point_load_contribution(self) -> np.ndarray:
        Fe = np.zeros(N_STANDARD_DOF, dtype=np.float64)
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
