# pre_processing/element_library/linear/euler_bernoulli/linear_warping_euler_bernoulli_3D.py
"""
2-node 3D Euler–Bernoulli beam with 7 DOF/node (standard six + warping intensity chi) for Vlasov non-uniform torsion.

**Tensors:** ``U_e`` (14,) node-major; DOFs 0–11 match linear EB, warping chi at indices 12–13 (one per node).
Per Gauss point ``B`` (7, 14) — first six rows from linear EB ``B`` on columns 0–11; row 6 is
``phi_x_prime = d(theta_x)/dx + d(chi)/dx`` (warping / bimoment strain). ``D`` (7, 7) — upper 6x6 is linear EB ``D``;
``D[6,6] = E*Gamma`` (bimoment stiffness). ``eps`` (7,), ``S = D @ eps``. ``detJ = L/2``.

**Weak forms (Gauss, xi in [-1, 1]):** ``K_e += B.T @ D @ B * w_g * detJ`` with ``B`` (7,14); distributed and point loads on first 12 standard DOFs; ``M_e`` uses extended ``N`` (14,6) per ``FORMULATION_DOCSTRING_STANDARDS.md``.

**Kinematics / shapes:** Registry ``LinearEulerBernoulliBeamElement3D`` supplies ``N`` (12, 6) per GP; warping DOFs use
linear Lagrange on the ``theta_x`` slot (``_N_14x6_at_xi``).

**Mass:** ``element_mass_matrix`` returns ``M_e`` (14, 14) — consistent mass on ``N`` (14, 6) with ``rho*Gamma`` on warping DOFs.

**Public API:** ``element_stiffness_matrix`` → ``ElementObject``; ``element_force_vector`` → ``ForceObject``;
``element_mass_matrix`` → ``MassObject``.
"""

import numpy as np
from typing import Optional, Tuple

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.shape_function_registry import get_shape_function_operator
from pre_processing.element_library.linear.euler_bernoulli.utilities.interpolate_loads import LoadInterpolationOperator
import logging

logger = logging.getLogger(__name__)

N_STANDARD_DOF = 12
N_WARPING_DOF = 2
N_DOF = N_STANDARD_DOF + N_WARPING_DOF
N_STRAIN = 7


class LinearWarpingEulerBernoulliBeamElement3D(Element1DBase):
    """
    2 nodes, 7 DOF/node; local ``x`` along chord.

    Notes
    -----
    **Contract/Diff vs 12-DOF EB:** ``U_e`` (14,) extends the standard (12,) packing with warping ``chi`` at indices 12–13;
    ``B`` (7,14) and ``D`` (7,7) embed linear EB on columns 0–11 and the upper 6x6 of ``D``; the seventh strain row
    (warping rate / bimoment strain) and ``D[6,6] = E*Gamma`` are the extension.

    Stiffness: ``K_e += B.T @ D @ B * w_g * detJ`` with ``B`` (7,14), ``D`` (7,7) as in the module docstring.
    ``Gamma`` from ``section_array[9]`` when length >= 10 (otherwise 0 — no bimoment stiffness).
    """

    element_type_name = "WarpingEulerBernoulli-3D"

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
            dof_per_node=7,
        )

        if quadrature_order is not None:
            self.quadrature_order = quadrature_order
        else:
            axial_order = int(self.element_array[3])
            bending_y_order = int(self.element_array[4])
            bending_z_order = int(self.element_array[5])
            torsion_order = int(self.element_array[8])
            load_order = int(self.element_array[9])
            self.quadrature_order = max(
                axial_order, bending_y_order, bending_z_order, torsion_order, load_order, 2
            )

        self.node_coords = self.grid_array
        self.L = np.linalg.norm(self.node_coords[1] - self.node_coords[0])
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = grid_dictionary["coordinates"][:, 0].min()
        self.x_global_end = grid_dictionary["coordinates"][:, 0].max()

        self._validate_element_properties()
        self._assert_logging_ready()

        self.shape_function_operator = get_shape_function_operator("LinearEulerBernoulliBeamElement3D", self.L)
        self.strain_displacement_operator = StrainDisplacementOperator(element_length=self.L)
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
    def Gamma(self) -> float:
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

    def _B_warping_row(self) -> np.ndarray:
        """Row 6 of B: φ_x′ = dθ_x/dx + dχ/dx. Linear: ±1/L for θ_x and warping DOFs."""
        row = np.zeros(N_DOF, dtype=np.float64)
        row[3] = -1.0 / self.L
        row[9] = 1.0 / self.L
        row[12] = -1.0 / self.L
        row[13] = 1.0 / self.L
        return row

    def _B_matrix_7x14(self, dN_dξ: np.ndarray, d2N_dξ2: np.ndarray) -> np.ndarray:
        """Build B (7, 14) from EB B (6, 12) plus warping row."""
        n_gauss = dN_dξ.shape[0]
        B_6x12 = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)
        B = np.zeros((n_gauss, N_STRAIN, N_DOF), dtype=np.float64)
        B[:, :6, :12] = B_6x12
        for g in range(n_gauss):
            B[g, 6, :] = self._B_warping_row()
        return B

    def _D_matrix_7x7(self) -> np.ndarray:
        """D (7,7): first 6×6 from EB (shear rows zero), D[6,6] = E·Γ."""
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
        xi_full, w_full = self.integration_points
        gauss_cache = []

        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N_12, dN_dξ_12, d2N_dξ2_12 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self._B_matrix_7x14(dN_dξ_12, d2N_dξ2_12)[0]
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
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} Warping EB K_e (14,14) ===")
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
        Fe = np.zeros(N_DOF, dtype=np.float64)

        gauss_cache = []
        if self.distributed_load_array.size > 0:
            Fe[:N_STANDARD_DOF], gauss_cache = self._compute_distributed_load_contribution()

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
            gauss_data=gauss_cache,
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
        Consistent mass: first 12 DOFs as straight Euler–Bernoulli; warping DOFs 12–13 with ρ·Γ
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
