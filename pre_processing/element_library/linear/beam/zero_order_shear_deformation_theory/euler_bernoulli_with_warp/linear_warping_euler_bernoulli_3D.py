# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/linear_warping_euler_bernoulli_3D.py
"""
2-node 3D Euler–Bernoulli beam with Vlasov warping (7 DOF/node: six standard + warping intensity χ per node).

**Tensors:** ``U_e`` (14,) node-major; ``K_e``, ``M_e`` (14, 14); ``F_e`` (14,). Per Gauss point ``B`` (7, 14),
``D`` (7, 7); ``ε`` (7,), ``S = D @ ε`` (7,). Voigt order for rows 0–5 follows
``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md``; row 6 is the warping-extension strain (see **Kinematics**).
``detJ = |J| = L/2``.

**Weak forms (Gauss, ξ in [-1, 1]):** ``K_e += B.T @ D @ B * w_g * detJ`` summed over Gauss points.
Distributed and point loads fill only the **first 12** standard DOFs using the registry EB shape tensor ``N`` (12, 6)
per point (same as linear EB). Consistent mass uses ``N`` extended to (14, 6); see ``utilities/shape_functions.py``.

**Kinematics:** Rows 0–5 of ``B`` on columns 0–11 match linear EB ``B`` (6, 12); column blocks 12–13 are the warping DOFs
for row 6. Strain row 6 is ``φ_x′ = ∂θ_x/∂x + ∂χ/∂x`` (bimoment-type warping rate), assembled as constant coefficients
on ``θ_x`` and ``χ`` DOFs in ``utilities/B_matrix.py``. Registry ``LinearEulerBernoulliBeamElement3D`` supplies
``N``, ``∂N/∂ξ``, ``∂²N/∂ξ²`` with batch slice (n_gp, 12, 6); χ rows for mass use
``utilities/shape_functions.extend_natural_shape_to_warping``.

**Constitutive:** ``D[:6, :6]`` is the linear EB material matrix; shear rows 3–4 remain zero. ``D[6, 6] = E·Γ`` with Γ
from ``section_array[9]`` when provided (else 0). See ``utilities/D_matrix.py``.

**Quadrature:** Gauss–Legendre order ``quadrature_order`` from the constructor or from ``element_array`` (axial,
bending_y, bending_z, torsion, load columns), same convention as linear EB.

**Public API:** ``element_stiffness_matrix`` → ``ElementObject``; ``element_force_vector`` → ``ForceObject``;
``element_mass_matrix`` → ``MassObject`` (consistent mass).
"""

import numpy as np
from typing import Optional, Tuple

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.shape_function_registry import get_shape_function_operator

from ..euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator
from ..euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator
from ..euler_bernoulli.utilities.interpolate_loads import LoadInterpolationOperator

from .utilities import (
    N_DOF,
    N_STANDARD_DOF,
    WarpingMaterialStiffnessOperator,
    WarpingStrainDisplacementOperator,
    extend_natural_shape_to_warping,
)
import logging

logger = logging.getLogger(__name__)


class LinearWarpingEulerBernoulliBeamElement3D(Element1DBase):
    """
    Two nodes; ``dof_per_node = 7``; local ``x`` along the chord (node 1 → node 2).

    **Identity / U_e:** Node-major displacement vector ``U_e ∈ ℝ^{14}``:

        [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, χ¹,
         u_x², u_y², u_z², θ_x², θ_y², θ_z², χ²]

    Indices 0–11 are the standard 12-DOF EB packing; indices 12–13 are warping intensities χ at nodes 1 and 2.

    **Tensors:** Same outer shapes as the module docstring: ``K_e``, ``M_e`` (14, 14); ``B`` (7, 14), ``D`` (7, 7);
    ``ε``, ``S`` (7,) with Voigt rows 0–5 per ``FORMULATION_DOCSTRING_STANDARDS`` and row 6 the warping strain.

    **Weak forms:** As the module **Weak forms** block; stiffness uses ``utilities/B_matrix.WarpingStrainDisplacementOperator``
    and ``utilities/D_matrix.WarpingMaterialStiffnessOperator``.

    **Constitutive / data:** ``Gamma`` (Γ) from ``section_array[9]`` when ``len(section_array) >= 10``, else 0.

    Notes
    -----
    **Contract vs 12-DOF EB:** The first 12 DOFs and first six strain rows/columns embed the linear EB baseline;
    the seventh strain row and DOFs 12–13 extend the model for warping, per
    ``FORMULATION_DOCSTRING_STANDARDS`` (extensions for ``(14,) U_e`` and ``(7, 14) B``).

    See Also
    --------
    LinearEulerBernoulliBeamElement3D
        12-DOF baseline in ``../euler_bernoulli/linear_euler_bernoulli_3D.py``.
    WarpingStrainDisplacementOperator, WarpingMaterialStiffnessOperator
        Warping ``B`` and ``D`` assembly in ``utilities/B_matrix.py`` and ``utilities/D_matrix.py``.
    extend_natural_shape_to_warping
        Extended ``N`` for consistent mass in ``utilities/shape_functions.py``.
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
        Voigt order for rows 0--5; extensions for ``(14,) U_e`` and ``(7, 14) B``.
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
        self.warping_strain_displacement_operator = WarpingStrainDisplacementOperator(
            element_length=self.L,
            base_strain_operator=self.strain_displacement_operator,
        )
        self.warping_material_stiffness_operator = WarpingMaterialStiffnessOperator(
            base_material_operator=self.material_stiffness_operator,
            youngs_modulus=self.E,
            warping_gamma=self.Gamma,
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

    def _precurvature_equivalent_load(self) -> np.ndarray:
        """``sum_g B.T @ D @ E_0 * w_g * detJ``; ``E_0`` is six beam rows plus zero warping strain."""
        if not np.any(self._E_0_voigt):
            return np.zeros(N_DOF, dtype=np.float64)
        E0 = np.concatenate([self._E_0_voigt, np.zeros(1, dtype=np.float64)])
        D = self.warping_material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        xi_full, w_full = self.integration_points
        f_e = np.zeros(N_DOF, dtype=np.float64)
        for xi_g, w_g in zip(xi_full, w_full):
            N_12, dN_dξ_12, d2N_dξ2_12 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.warping_strain_displacement_operator.physical_coordinate_form(dN_dξ_12, d2N_dξ2_12)[0]
            f_e += (B.T @ D @ E0) * w_g * detJ
        return f_e

    def element_stiffness_matrix(self):
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        Ke = np.zeros((N_DOF, N_DOF), dtype=np.float64)
        D = self.warping_material_stiffness_operator.assembly_form()
        xi_full, w_full = self.integration_points
        gauss_cache = []

        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N_12, dN_dξ_12, d2N_dξ2_12 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.warping_strain_displacement_operator.physical_coordinate_form(dN_dξ_12, d2N_dξ2_12)[0]
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

        Fe += self._precurvature_equivalent_load()

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
            Ng = extend_natural_shape_to_warping(N12[0], xi_g)
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
