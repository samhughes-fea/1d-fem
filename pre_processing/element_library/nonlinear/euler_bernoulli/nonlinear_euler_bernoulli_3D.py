# pre_processing/element_library/nonlinear/euler_bernoulli/nonlinear_euler_bernoulli_3D.py
"""
2-node 3D Euler-Bernoulli beam, Total Lagrangian (TL) geometric nonlinearity.

**12-DOF path:** ``U_e`` (12,); per Gauss ``B_lin`` (6, 12), ``B_nl`` (6, 12); ``D`` (6, 6); ``K_T`` (12, 12).

**14-DOF path (Vlasov warping):** When ``beam_warping.mesh_uses_warping_dof`` is true for this element
(``[warping]`` in ``element.txt`` / ``element_dictionary["warping"]``), ``U_e`` (14,); ``B_lin`` (7, 14),
``B_nl`` seventh row zero; ``D`` (7, 7) with ``D[6,6]=E·Γ_eff``; ``K_sigma`` embeds the 12×12 geometric
stiffness on the first 12 DOFs. Same policy as ``LinearEulerBernoulliBeamElement3D`` + warping.

**Weak forms:** See module-level notes in the removed standalone warping class (now merged here).

**See Also:** ``docs/element_library/total_lagrangian_beam_formulation.md``;
``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md``; ``docs/conventions/JOB_INPUT_BEAM_WARPING.md``.
"""

import logging
from typing import List, Optional, Tuple

import numpy as np

from pre_processing.element_library.beam_warping import (
    beam_warping_policy,
    enforce_strict_section_gamma,
    mesh_uses_warping_dof,
    section_gamma_from_section_array,
    warn_if_degenerate_warping_stiffness,
)
from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli_with_warp.utilities import (
    WarpingMaterialStiffnessOperator,
    WarpingStrainDisplacementOperator,
    extend_natural_shape_to_warping,
)
from pre_processing.element_library.shape_function_registry import get_shape_function_operator
from pre_processing.element_library.nonlinear.euler_bernoulli.utilities import (
    GreenLagrangeStrainOperator,
    StressResultantOperator,
    GeometricStiffnessOperator,
)

logger = logging.getLogger(__name__)

_NL_W_STD = 12
_NL_W_DOF = 14


def _B_nl_12_to_14(B_nl_6_12: np.ndarray) -> np.ndarray:
    B = np.zeros((7, _NL_W_DOF), dtype=np.float64)
    B[:6, :_NL_W_STD] = B_nl_6_12
    return B


class NonlinearEulerBernoulliBeamElement3D(Element1DBase):
    """
    TL EB element: Green-Lagrange strain, ``S = D @ E``, material and geometric tangents from Gauss sums.

    Optional **14 local DOFs** when the mesh allocates warping (χ) per node — same ``beam_warping`` policy as
    linear EB; use ``NonlinearEulerBernoulliBeamElement3D`` + ``[warping]`` in ``element.txt``.

    Notes
    -----
    Operators: ``green_lagrange_strain_operator``, ``stress_resultant_operator``, ``geometric_stiffness_operator``
    compose linear EB ``B_lin``, ``D``, and shapes. Shear rows in the engineering Voigt layout stay in the EB pattern;
    ``N``, ``M_y``, ``M_z`` for ``K_sigma`` from ``section_forces_from_strain(E, D)``. Full weak forms: module docstring.
    """

    element_type_name = "Euler-Bernoulli-3D-Nonlinear"

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
        """
        Initialize the nonlinear Euler-Bernoulli beam element.

        Warping DOFs follow ``mesh_uses_warping_dof(element_dictionary)`` (``[warping]`` column or legacy type
        name containing ``\"Warping\"`` when the column is absent).
        """
        idx = int(np.where(np.asarray(element_dictionary["ids"]) == element_id)[0][0])
        etype_str = str(element_dictionary["types"][idx])
        warp_mesh = mesh_uses_warping_dof(element_dictionary)
        dof_pn = 7 if warp_mesh else 6
        self._n_dof = dof_pn * 2

        super().__init__(
            element_id=element_id,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir,
            dof_per_node=dof_pn,
        )
        self._mesh_element_dictionary = element_dictionary
        self._mesh_grid_dictionary = grid_dictionary
        self._mesh_material_dictionary = material_dictionary
        self._mesh_section_dictionary = section_dictionary
        # Quadrature order: from element_array when not provided (same 7-column convention as linear EB)
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
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = float(grid_dictionary["coordinates"][:, 0].min())
        self.x_global_end = float(grid_dictionary["coordinates"][:, 0].max())
        self._warp_policy = beam_warping_policy(
            element_dictionary,
            idx,
            etype_str,
            section_gamma_from_section_array(self.section_array),
        )
        self._warp_mesh = self._warp_policy.mesh_allocates_chi_dof
        self._warp_stiff = self._warp_policy.warping_stiffness_on
        self._validate_element_properties()
        self._assert_logging_ready()

        # Linear operators (from linear/euler_bernoulli)
        self.shape_function_operator = get_shape_function_operator(self.__class__.__name__, self.L)
        self.strain_displacement_operator = StrainDisplacementOperator(element_length=self.L)
        self.material_stiffness_operator = MaterialStiffnessOperator(
            youngs_modulus=self.E,
            shear_modulus=self.G,
            cross_section_area=self.A,
            moment_inertia_y=self.I_y,
            moment_inertia_z=self.I_z,
            torsion_constant=self.J_t,
        )
        # Total Lagrangian operators
        self.green_lagrange_strain_operator = GreenLagrangeStrainOperator(element_length=self.L)
        self.stress_resultant_operator = StressResultantOperator()
        self.geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)

        # Cache K_0 (material stiffness, same as linear K_e); 12×12 path only
        self._K_0: np.ndarray | None = None

        self.warping_strain_displacement_operator: Optional[WarpingStrainDisplacementOperator] = None
        self.warping_material_stiffness_operator: Optional[WarpingMaterialStiffnessOperator] = None
        if self._n_dof == 14:
            enforce_strict_section_gamma(
                element_dictionary=element_dictionary,
                element_id=element_id,
                stiffness_on=self._warp_stiff,
                section_gamma=self._warp_policy.gamma_section,
            )
            warn_if_degenerate_warping_stiffness(
                stiffness_on=self._warp_stiff,
                section_gamma=self._warp_policy.gamma_section,
                element_id=element_id,
            )
            gamma_eff = self._warp_policy.gamma_effective
            self.warping_strain_displacement_operator = WarpingStrainDisplacementOperator(
                element_length=self.L,
                base_strain_operator=self.strain_displacement_operator,
            )
            self.warping_material_stiffness_operator = WarpingMaterialStiffnessOperator(
                base_material_operator=self.material_stiffness_operator,
                youngs_modulus=self.E,
                warping_gamma=gamma_eff,
            )

    def _validate_element_properties(self) -> None:
        """Validate critical element properties and log geometry."""
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4:
            raise ValueError("Material array not properly initialised")
        if self._warp_mesh:
            if self.section_array.size not in (5, 7, 9, 10):
                raise ValueError(
                    "Section array must have 5, 7, 9 or 10 entries when warping DOFs are used"
                )
        elif self.section_array.size not in (5, 7):
            raise ValueError("Material/section arrays not properly initialised")
        conn = tuple(self.element_array[1:3])
        logger.debug(
            "Element %s geometry initialised\n"
            "• Connectivity: %s\n• Length: %.4e\n• Start/End X: %.4e / %.4e",
            self.element_id, conn, self.L, self.x_start, self.x_end,
        )

    @property
    def Gamma(self) -> float:
        """Vlasov warping constant Γ [m⁶] from ``section_array`` index 9 when present."""
        if self.section_array.size >= 10:
            return float(self.section_array[9])
        return 0.0

    @property
    def A(self) -> float:
        """Cross-sectional area (m²)."""
        return float(self.section_array[0])

    @property
    def I_y(self) -> float:
        """Moment of inertia about y-axis (m⁴)."""
        return float(self.section_array[2])

    @property
    def I_z(self) -> float:
        """Moment of inertia about z-axis (m⁴)."""
        return float(self.section_array[3])

    @property
    def J_t(self) -> float:
        """Torsional constant (m⁴)."""
        return float(self.section_array[4])

    @property
    def E(self) -> float:
        """Young's modulus (Pa)."""
        return float(self.material_array[0])

    @property
    def G(self) -> float:
        """Shear modulus (Pa)."""
        return float(self.material_array[1])

    @property
    def jacobian_determinant(self) -> float:
        """Jacobian |J| = L/2 of the element coordinate mapping."""
        return self.L / 2.0

    @property
    def integration_points(self) -> Tuple[np.ndarray, np.ndarray]:
        """Gauss–Legendre quadrature points and weights."""
        return np.polynomial.legendre.leggauss(self.quadrature_order)

    def element_mass_matrix(self):
        """
        Reference-configuration consistent mass (same as linear Euler–Bernoulli) for modal/dynamic.
        """
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        if self._n_dof == 12:
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

        rho = float(self.material_array[3])
        mu = np.zeros(_NL_W_DOF, dtype=np.float64)
        for i in (0, 1, 2, 6, 7, 8):
            mu[i] = rho * self.A
        for i in (3, 9):
            mu[i] = rho * self.J_t
        for i in (4, 10):
            mu[i] = rho * self.I_y
        for i in (5, 11):
            mu[i] = rho * self.I_z
        g_m = self._warp_policy.gamma_effective
        mu[12] = rho * g_m
        mu[13] = rho * g_m
        M_e = np.zeros((_NL_W_DOF, _NL_W_DOF), dtype=np.float64)
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        for xi_g, w_g in zip(xi, w):
            N12, _, _ = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            Ng = extend_natural_shape_to_warping(N12[0], xi_g)
            for i in range(_NL_W_DOF):
                for j in range(_NL_W_DOF):
                    mij = 0.5 * (mu[i] + mu[j])
                    M_e[i, j] += mij * float(np.dot(Ng[i, :], Ng[j, :])) * w_g * detJ
        return MassObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            M_e=M_e,
        )

    def _get_K_0(self) -> np.ndarray:
        """Material stiffness ``K_0 = sum_g B_lin.T @ D @ B_lin * w_g * detJ`` (cached; linear ``B_lin``).

        Returns
        -------
        np.ndarray
            Shape (12, 12); same as linear EB element stiffness at U=0.
        """
        if self._n_dof != 12:
            raise NotImplementedError("_get_K_0 is defined for the 12-DOF EB path only.")
        if self._K_0 is not None:
            return self._K_0
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        K_0 = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)[0]
            K_0 += B.T @ D @ B * w_g * detJ
        self._K_0 = K_0
        return K_0

    def linear_geometric_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """
        Geometric stiffness for **linear** buckling theory about ``U_e``, matching
        :class:`LinearEulerBernoulliBeamElement3D` (for ``modal.buckling_prestress=nonlinear_static``).
        """
        from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
            LinearEulerBernoulliBeamElement3D,
        )

        ed = dict(self._mesh_element_dictionary)
        n_t = len(ed["types"])
        ed["types"] = np.array(["LinearEulerBernoulliBeamElement3D"] * n_t, dtype=object)
        lin = LinearEulerBernoulliBeamElement3D(
            element_id=self.element_id,
            element_dictionary=ed,
            grid_dictionary=self._mesh_grid_dictionary,
            material_dictionary=self._mesh_material_dictionary,
            section_dictionary=self._mesh_section_dictionary,
            point_load_array=self.point_load_array,
            distributed_load_array=self.distributed_load_array,
            job_results_dir=self.job_results_dir,
            quadrature_order=self.quadrature_order,
        )
        return lin.linear_geometric_stiffness_matrix(U_e)

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """
        Tangent stiffness ``K_T = K_mat(U_e) + K_sigma(U_e)``. ``K_mat += sum_g B_tot.T @ D @ B_tot * w_g * detJ``
        with ``B_tot = B_lin + B_nl`` (full nonlinear curvature); ``K_sigma`` is geometric stiffness.

        Parameters
        ----------
        U_e : np.ndarray
            Element displacement vector, shape (12,) or (14,) when warping DOFs are active.

        Returns
        -------
        np.ndarray
            Tangent stiffness matrix, shape (12, 12) or (14, 14).
        """
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
                dN_dx[:, 3] = dN_dξ[:, 3] * dξ_dx  # torsion
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
            B_lin_12 = self.green_lagrange_strain_operator.linearized_strain_displacement(
                dN_dx[0], d2N_dx2[0]
            )
            B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx[0], d2N_dx2[0], U12
            )
            B_lin_7 = self.warping_strain_displacement_operator.physical_coordinate_form(
                dN_dξ, d2N_dξ2
            )[0]
            B_nl_7 = _B_nl_12_to_14(B_nl_12)
            B_tot = B_lin_7 + B_nl_7
            K_mat += B_tot.T @ D7 @ B_tot * w_g * detJ
        return K_mat + K_sigma

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        """
        Internal (residual) force from Green–Lagrange strain and ``B_tot = B_lin + B_nl``.

        Accumulates ``F_int += B_tot.T @ S * w_g * detJ`` over Gauss points with ``S = D @ E``,
        ``E = E_lin + E_nl`` (Total Lagrangian EB). Matches the material tangent ``K_mat`` which
        uses the same ``B_tot``.
        """
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
                B_tot = B_lin + B_nl
                E_lin = self.green_lagrange_strain_operator.strain_linear_part(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(
                    dN_dx[0], d2N_dx2[0], U_e
                )
                E = E_lin + E_nl
                S = D @ (E - self._E_0_voigt)
                F_int += (B_tot.T @ S) * w_g * detJ
            return F_int

        U_e = np.asarray(U_e, dtype=np.float64).reshape(_NL_W_DOF)
        U12 = U_e[:_NL_W_STD]
        assert self.warping_material_stiffness_operator is not None
        assert self.warping_strain_displacement_operator is not None
        D7 = self.warping_material_stiffness_operator.assembly_form()
        D6 = D7[:6, :6]
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
            B_lin_12 = self.green_lagrange_strain_operator.linearized_strain_displacement(
                dN_dx[0], d2N_dx2[0]
            )
            B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx[0], d2N_dx2[0], U12
            )
            B_lin_7 = self.warping_strain_displacement_operator.physical_coordinate_form(
                dN_dξ, d2N_dξ2
            )[0]
            B_nl_7 = _B_nl_12_to_14(B_nl_12)
            B_tot = B_lin_7 + B_nl_7
            E_lin = self.green_lagrange_strain_operator.strain_linear_part(dN_dx[0], d2N_dx2[0], U12)
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx[0], d2N_dx2[0], U12)
            E6 = E_lin + E_nl
            E_warp = float(np.dot(B_lin_7[6], U_e))
            E = np.concatenate([E6, np.array([E_warp], dtype=np.float64)])
            S = D7 @ (E - E0_7)
            F_int += (B_tot.T @ S) * w_g * detJ
        return F_int

    def strain_at_gauss_points(self, U_e: np.ndarray) -> List[np.ndarray]:
        """
        Return strain E_lin + E_nl at each integration point (same order as gauss_data).

        Parameters
        ----------
        U_e : np.ndarray
            Element displacement vector, shape (12,) or (14,); only the first 12 entries affect EB strains.

        Returns
        -------
        List[np.ndarray]
            One strain vector per Gauss point, each shape (6,) or equivalent.
        """
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
            E = E_lin + E_nl
            result.append(np.asarray(E, dtype=np.float64))
        return result

    def element_stiffness_matrix(self):
        """
        Return ElementObject with K_e = initial tangent at U=0 (same as linear K_e).

        Caches gauss_data (B, D, shape_functions, shape_derivatives per Gauss point)
        and evaluate_shape_functions for post-processing. B2 shape-function coefficients
        are set (same monomial basis as linear EB) for save/load evaluation; see RESULTS_DESIGN.md.

        Returns
        -------
        ElementObject
            K_e = tangent_stiffness_matrix(0), gauss_data, integration_scheme, evaluate_shape_functions.
        """
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData

        self._assert_logging_ready()
        if self._n_dof == 12:
            U_zero = np.zeros(12, dtype=np.float64)
            K_e = self.tangent_stiffness_matrix(U_zero)
            xi, w = self.integration_points
            D = self.material_stiffness_operator.assembly_form()
            detJ = self.jacobian_determinant
            L_arr = np.array([[self.L]], dtype=np.float64)
            if self.logger_operator:
                self.logger_operator.log_text(
                    "stiffness",
                    f"\n=== Element {self.element_id} Stiffness (initial tangent at U=0) ===",
                )
                self.logger_operator.log_matrix("stiffness", L_arr, {"name": f"Element length  L  (1,1)"})
                self.logger_operator.log_matrix("stiffness", D, {"name": f"Material stiffness matrix  D  {D.shape}"})
            gauss_cache = []
            for g, (xi_g, w_g) in enumerate(zip(xi, w)):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)[0]
                Ke_contribution = B.T @ D @ B * w_g * detJ
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
                if self.logger_operator:
                    self._log_gauss_point_stiffness(g, xi_g, w_g, dN_dξ, d2N_dξ2, B, Ke_contribution)
            if self.logger_operator:
                self.logger_operator.log_matrix("stiffness", K_e, {"name": "Initial tangent K_e"})
                self.logger_operator.flush("stiffness")
            op = self.shape_function_operator
            evaluate_shape_functions = lambda xi_val: op.natural_coordinate_form(np.asarray(xi_val))
            N_coeffs, dN_coeffs, d2N_coeffs = op.build_shape_function_coefficients_b2()
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

        U_zero = np.zeros(_NL_W_DOF, dtype=np.float64)
        K_e = self.tangent_stiffness_matrix(U_zero)
        xi, w = self.integration_points
        assert self.warping_material_stiffness_operator is not None
        D = self.warping_material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        gauss_cache = []
        for g, (xi_g, w_g) in enumerate(zip(xi, w)):
            _N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            assert self.warping_strain_displacement_operator is not None
            B = self.warping_strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)[0]
            Ke_contribution = B.T @ D @ B * w_g * detJ
            gauss_cache.append(
                StiffnessGaussPointData(
                    xi=float(xi_g),
                    weight=float(w_g),
                    B_matrix=B.copy(),
                    D_matrix=D.copy(),
                    jacobian=float(detJ),
                    shape_functions=_N.copy(),
                    shape_derivatives=dN_dξ.copy(),
                )
            )
            if self.logger_operator:
                self._log_gauss_point_stiffness_warp(g, xi_g, w_g, dN_dξ, d2N_dξ2, B, Ke_contribution)
        if self.logger_operator:
            self.logger_operator.log_matrix("stiffness", K_e, {"name": "Initial tangent K_e (14×14)"})
            self.logger_operator.flush("stiffness")
        op = self.shape_function_operator
        evaluate_shape_functions = lambda xi_val: op.natural_coordinate_form(np.asarray(xi_val))
        N_coeffs, dN_coeffs, d2N_coeffs = op.build_shape_function_coefficients_b2()
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

    def element_force_vector(self):
        """
        Compute the element force vector (external loads): distributed and point loads.

        Returns
        -------
        ForceObject
            F_e and gauss_data; same convention as linear EB.

        Notes
        -----
        Combines distributed load ``F_dist += sum_g w_g * N.T @ q * detJ`` and point loads
        F_point = N(x_p)ᵀ P at load locations.
        """
        from pre_processing.element_library.gauss_point_data import ForceObject, ForceGaussPointData
        from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.interpolate_loads import LoadInterpolationOperator

        self._assert_logging_ready()
        if self.logger_operator:
            self.logger_operator.log_text(
                "force",
                f"\n=== Element {self.element_id} Force Vector Computation ===",
            )
        if self._n_dof == 12:
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
                if self.logger_operator:
                    self._log_distributed_loads(xi, w, N, q_gauss, Fe_dist)
            if self.point_load_array.size > 0:
                from pre_processing.element_library.point_load_utils import add_phased_increment, point_load_phase_rad

                for load in self.point_load_array:
                    x_p = float(load[0])
                    F_p = load[3:9].astype(np.float64)
                    if x_start <= x_p <= x_start + self.L:
                        xi_p = 2 * (x_p - x_start) / self.L - 1
                        N_p = self.shape_function_operator.natural_coordinate_form(np.array([xi_p]))[0][0]
                        Fe_trans = N_p[[0, 1, 2, 6, 7, 8], :3] @ F_p[:3]
                        Fe_rot = N_p[[3, 4, 5, 9, 10, 11], 3:] @ F_p[3:]
                        inc = np.zeros_like(Fe)
                        inc[[0, 1, 2, 6, 7, 8]] = Fe_trans
                        inc[[3, 4, 5, 9, 10, 11]] = Fe_rot
                        Fe = add_phased_increment(Fe, inc, point_load_phase_rad(load))
                        if self.logger_operator:
                            self._log_point_load(x_p, xi_p, F_p, N_p, Fe_trans, Fe_rot)
            if self.logger_operator:
                self.logger_operator.log_matrix("force", Fe.reshape(1, -1), {"name": "Final Force Vector"})
                self.logger_operator.flush("force")
            return ForceObject(
                element_id=self.element_id,
                element_type=self.element_type_name,
                F_e=Fe,
                gauss_data=gauss_cache,
                point_loads=self.point_load_array.copy() if self.point_load_array.size else None,
            )

        Fe = np.zeros(_NL_W_DOF, dtype=np.float64)
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
            Fe[:_NL_W_STD] += Fe_dist
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
            if self.logger_operator:
                self._log_distributed_loads(xi, w, N, q_gauss, Fe_dist)
        if self.point_load_array.size > 0:
            from pre_processing.element_library.point_load_utils import add_phased_increment, point_load_phase_rad

            for load in self.point_load_array:
                x_p = float(load[0])
                F_p = load[3:9].astype(np.float64)
                if x_start <= x_p <= x_start + self.L:
                    xi_p = 2 * (x_p - x_start) / self.L - 1
                    N_p = self.shape_function_operator.natural_coordinate_form(np.array([xi_p]))[0][0]
                    Fe_trans = N_p[[0, 1, 2, 6, 7, 8], :3] @ F_p[:3]
                    Fe_rot = N_p[[3, 4, 5, 9, 10, 11], 3:] @ F_p[3:]
                    inc = np.zeros_like(Fe)
                    inc[[0, 1, 2, 6, 7, 8]] = Fe_trans
                    inc[[3, 4, 5, 9, 10, 11]] = Fe_rot
                    Fe = add_phased_increment(Fe, inc, point_load_phase_rad(load))
                    if self.logger_operator:
                        self._log_point_load(x_p, xi_p, F_p, N_p, Fe_trans, Fe_rot)
        if self.logger_operator:
            self.logger_operator.log_matrix("force", Fe.reshape(1, -1), {"name": "Final Force Vector"})
            self.logger_operator.flush("force")
        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=gauss_cache,
            point_loads=self.point_load_array.copy() if self.point_load_array.size else None,
        )

    def _xi_to_x(self, xi: float) -> float:
        """Convert natural coordinate ξ to physical position x."""
        return (xi + 1) * (self.L / 2) + self.x_start

    def _log_gauss_point_stiffness(
        self,
        gp_idx: int,
        xi: float,
        weight: float,
        dN_dξ: np.ndarray,
        d2N_dξ2: np.ndarray,
        B: np.ndarray,
        contribution: np.ndarray,
    ) -> None:
        """Log detailed stiffness integration data for one Gauss point (same format as linear EB)."""
        if dN_dξ.shape[-2:] != (12, 6):
            raise ValueError(f"dN_dξ shape mismatch: {dN_dξ.shape} ≠ (12, 6)")
        if d2N_dξ2.shape[-2:] != (12, 6):
            raise ValueError(f"d2N_dξ2 shape mismatch: {d2N_dξ2.shape} ≠ (12, 6)")
        if B.shape[-2:] != (6, 12):
            raise ValueError(f"B-matrix shape mismatch: {B.shape} ≠ (*, 6, 12)")
        if contribution.shape != (12, 12):
            raise ValueError(f"Contribution shape mismatch: {contribution.shape} ≠ (12, 12)")
        metadata = {"name": f"GP {gp_idx + 1}", "precision": 6, "max_line_width": 120}
        self.logger_operator.log_text(
            "stiffness",
            f"\nGP {gp_idx + 1}/{self.quadrature_order}: "
            f"ξ = {xi:.6f},  x = {self._xi_to_x(xi):.6e},  w = {weight:.6e}",
        )
        self.logger_operator.log_matrix(
            "stiffness", dN_dξ, {**metadata, "name": f"Shape-function derivative  dN/dξ  {dN_dξ.shape}"}
        )
        self.logger_operator.log_matrix(
            "stiffness", d2N_dξ2, {**metadata, "name": f"Second derivative  d²N/dξ²  {d2N_dξ2.shape}"}
        )
        self.logger_operator.log_matrix(
            "stiffness", B, {**metadata, "name": f"Strain-displacement matrix  B  {B.shape}"}
        )
        self.logger_operator.log_matrix(
            "stiffness", contribution, {**metadata, "name": f"Gauss-point contribution  BᵀDB  {contribution.shape}"}
        )

    def _log_gauss_point_stiffness_warp(
        self,
        gp_idx: int,
        xi: float,
        weight: float,
        dN_dξ: np.ndarray,
        d2N_dξ2: np.ndarray,
        B: np.ndarray,
        contribution: np.ndarray,
    ) -> None:
        if dN_dξ.shape[-2:] != (12, 6):
            raise ValueError(f"dN_dξ shape mismatch: {dN_dξ.shape} ≠ (12, 6)")
        if d2N_dξ2.shape[-2:] != (12, 6):
            raise ValueError(f"d2N_dξ2 shape mismatch: {d2N_dξ2.shape} ≠ (12, 6)")
        if B.shape[-2:] != (7, 14):
            raise ValueError(f"B-matrix shape mismatch: {B.shape} ≠ (*, 7, 14)")
        if contribution.shape != (14, 14):
            raise ValueError(f"Contribution shape mismatch: {contribution.shape} ≠ (14, 14)")
        metadata = {"name": f"GP {gp_idx + 1}", "precision": 6, "max_line_width": 120}
        self.logger_operator.log_text(
            "stiffness",
            f"\nGP {gp_idx + 1}/{self.quadrature_order}: "
            f"ξ = {xi:.6f},  x = {self._xi_to_x(xi):.6e},  w = {weight:.6e}",
        )
        self.logger_operator.log_matrix(
            "stiffness", dN_dξ, {**metadata, "name": f"Shape-function derivative  dN/dξ  {dN_dξ.shape}"}
        )
        self.logger_operator.log_matrix(
            "stiffness", d2N_dξ2, {**metadata, "name": f"Second derivative  d²N/dξ²  {d2N_dξ2.shape}"}
        )
        self.logger_operator.log_matrix(
            "stiffness", B, {**metadata, "name": f"Strain-displacement matrix  B  {B.shape}"}
        )
        self.logger_operator.log_matrix(
            "stiffness", contribution, {**metadata, "name": f"Gauss-point contribution  BᵀDB  {contribution.shape}"}
        )

    def _log_distributed_loads(
        self,
        xi: np.ndarray,
        weights: np.ndarray,
        N: np.ndarray,
        q: np.ndarray,
        Fe: np.ndarray,
    ) -> None:
        """Log distributed load integration details (same format as linear EB)."""
        if not self.logger_operator:
            return
        if N.shape[1:] != (12, 6):
            raise ValueError(f"N shape mismatch: {N.shape} != (n_pts,12,6)")
        if q.shape[1] != 6:
            raise ValueError(f"Load vector shape: {q.shape} != (n_pts,6)")
        if Fe.shape != (12,):
            raise ValueError(f"Fe result shape: {Fe.shape} != (12,)")
        metadata = {"precision": 6, "max_line_width": 100}
        self.logger_operator.log_text("force", "\n=== Distributed Loads ===")
        for gp, (xi_g, w_g) in enumerate(zip(xi, weights)):
            gp_meta = {**metadata, "name": f"GP {gp+1}"}
            self.logger_operator.log_matrix("force", N[gp], {**gp_meta, "name": f"N {N[gp].shape}"})
            self.logger_operator.log_matrix("force", q[gp], {**gp_meta, "name": f"q {q[gp].shape}"})
        self.logger_operator.log_matrix("force", Fe.reshape(1, -1), {**metadata, "name": f"Total Fe {Fe.shape}"})

    def _log_point_load(
        self,
        x: float,
        xi: float,
        F: np.ndarray,
        N: np.ndarray,
        trans: np.ndarray,
        rot: np.ndarray,
    ) -> None:
        """Log point load application (same format as linear EB)."""
        if not self.logger_operator:
            return
        if N.shape != (12, 6):
            raise ValueError(f"N shape mismatch: {N.shape} != (12,6)")
        if trans.shape != (6,):
            raise ValueError(f"Translation vector shape: {trans.shape} != (6,)")
        if rot.shape != (6,):
            raise ValueError(f"Rotation vector shape: {rot.shape} != (6,)")
        metadata = {"precision": 6, "max_line_width": 120}
        self.logger_operator.log_text(
            "force",
            f"\n=== Point Load @ x={x:.6e} ===\n"
            f"Natural ξ={xi:.6f}, Element Range: {self.x_start:.6e}-{self.x_end:.6e}",
        )
        self.logger_operator.log_matrix("force", F.reshape(-1, 1), {**metadata, "name": "Force Vector [6×1]"})
        self.logger_operator.log_matrix("force", N, {**metadata, "name": f"Shape Functions {N.shape}"})
        self.logger_operator.log_matrix("force", trans.reshape(-1, 1), {**metadata, "name": "Translations [6×1]"})
        self.logger_operator.log_matrix("force", rot.reshape(-1, 1), {**metadata, "name": "Rotations [6×1]"})
