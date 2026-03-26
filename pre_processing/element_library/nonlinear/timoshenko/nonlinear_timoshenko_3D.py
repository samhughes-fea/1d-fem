# pre_processing/element_library/nonlinear/timoshenko/nonlinear_timoshenko_3D.py
"""
2-node 3D Timoshenko beam, Total Lagrangian geometric nonlinearity.

**Tensors:** ``U_e`` (12,); reference ``B`` (6, 12), ``D`` (6, 6) from linear Timoshenko (``kappa*G*A`` shear).
``E = E_lin + E_nl`` (6,) with ``include_shear=True`` in ``GreenLagrangeStrainOperator``; ``S = D @ E``.
``F_int += B.T @ S * w_g * detJ`` (``B`` = linear Timoshenko strain-displacement operator in physical coordinates);
``K_T = K_0 + K_sigma`` — ``K_0`` (12, 12) linear Timoshenko stiffness (selective integration); ``K_sigma`` (12, 12) from
``GeometricStiffnessOperator`` using ``N``, ``M_y``, ``M_z``. ``detJ = L/2``.

**Weak forms (Gauss, xi in [-1, 1]):** internal force and ``K_0`` as Gauss sums; ``K_sigma`` same kernel as nonlinear EB but with Timoshenko ``dN/dx``;
``F_dist``, ``M_e`` as linear Timoshenko.

**Quadrature:** ``max(quadrature_order, 2)`` default; matches linear Timoshenko shear needs.

**Public API:** See ``tangent_stiffness_matrix`` and internal force methods in source.

**See Also:** ``docs/element_library/total_lagrangian_beam_formulation.md``; ``FORMULATION_DOCSTRING_STANDARDS.md`` (TL utilities).
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


class NonlinearTimoshenkoBeamElement3D(Element1DBase):
    """
    TL Timoshenko: linear ``B``, ``D`` for ``K_0``; Green-Lagrange ``E_nl`` with shear; ``K_sigma`` from ``N``, ``M_y``, ``M_z``.

    Notes
    -----
    Module docstring lists tensor shapes and Gauss sums. ``GeometricStiffnessOperator`` is shared with nonlinear EB family (Timoshenko shapes).
    """

    element_type_name = "Timoshenko-3D-Nonlinear"

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
        """
        Initialize the nonlinear Timoshenko beam element.

        Parameters
        ----------
        element_id : int
            Element ID in the mesh.
        element_dictionary : dict
            Element connectivity and type data (ids, connectivity, types, etc.).
        grid_dictionary : dict
            Node coordinates (key "coordinates").
        section_dictionary : dict
            Cross-section properties (A, I_y, I_z, J_t, kappa if size >= 7).
        material_dictionary : dict
            Material properties (E, G, nu, rho).
        point_load_array : np.ndarray
            Point loads array (Nx9): x, y, z, Fx, Fy, Fz, Mx, My, Mz.
        distributed_load_array : np.ndarray
            Distributed loads array (Nx9).
        job_results_dir : str
            Directory for job logs and output.
        quadrature_order : int, optional
            Gauss–Legendre quadrature order; uses max(quadrature_order, 2) so shear terms
            are integrated (default 3).

        Notes
        -----
        x_start, x_end, x_global_start, and x_global_end are set from node_coords
        and grid_dictionary (same convention as linear elements).
        """
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
        """Validate critical element properties and log geometry."""
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
    def kappa(self) -> float:
        """Shear correction factor (default 5/6 if not in section)."""
        return float(self.section_array[5]) if self.section_array.size >= 7 else 5.0 / 6.0

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
        Reference-configuration consistent mass (same as linear Timoshenko) for modal/dynamic.
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
        """Material stiffness ``K_0 = sum_g B.T @ D @ B * w_g * detJ`` (linear Timoshenko ``B``, full quadrature).

        Returns
        -------
        np.ndarray
            Shape (12, 12); same as linear Timoshenko element stiffness at U=0.
        """
        if self._K_0 is not None:
            return self._K_0
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        K_0 = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            K_0 += B.T @ D @ B * w_g * detJ
        self._K_0 = K_0
        return K_0

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """
        Tangent stiffness K_T = K_0 + K_σ(U_e) in global (element) coordinates.

        Parameters
        ----------
        U_e : np.ndarray
            Element displacement vector, shape (12,): [u_x1, u_y1, u_z1, θ_x1, θ_y1, θ_z1, ...].

        Returns
        -------
        np.ndarray
            Tangent stiffness matrix, shape (12, 12).
        """
        K_0 = self._get_K_0()
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        dξ_dx = 2.0 / self.L
        n_g = len(xi)
        N_gp = np.zeros(n_g, dtype=np.float64)
        M_y_gp = np.zeros(n_g, dtype=np.float64)
        M_z_gp = np.zeros(n_g, dtype=np.float64)
        dN_dx_list = []
        for k, (xi_g, w_g) in enumerate(zip(xi, w)):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            dN_dx = dN_dξ.copy() * dξ_dx
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            E_lin = (B @ U_e).ravel()
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx[0], U_e)
            E = E_lin + E_nl
            N_i, M_y_i, M_z_i = self.stress_resultant_operator.section_forces_from_strain(E, D)
            N_gp[k] = N_i
            M_y_gp[k] = M_y_i
            M_z_gp[k] = M_z_i
            dN_dx_list.append(dN_dx[0])
        dN_dx_arr = np.stack(dN_dx_list, axis=0)
        K_sigma = self.geometric_stiffness_operator.assemble_K_sigma(
            N_gp, M_y_gp, M_z_gp, w, dN_dx_arr, detJ
        )
        return K_0 + K_sigma

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        """
        Internal (residual) force using linear Timoshenko ``B`` in physical coordinates.

        Accumulates ``F_int += B.T @ S * w_g * detJ`` over Gauss points with ``S = D @ E``,
        ``E = E_lin + E_nl``.

        Parameters
        ----------
        U_e : np.ndarray
            Element displacement vector, shape (12,).

        Returns
        -------
        np.ndarray
            Internal force vector, shape (12,).
        """
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        F_int = np.zeros(12, dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            E_lin = (B @ U_e).ravel()
            dξ_dx = 2.0 / self.L
            dN_dx = dN_dξ[0] * dξ_dx
            E_nl = self.green_lagrange_strain_operator.strain_nonlinear_part(dN_dx, U_e)
            E = E_lin + E_nl
            S = D @ E
            F_int += (B.T @ S) * w_g * detJ
        return F_int

    def strain_at_gauss_points(self, U_e: np.ndarray) -> List[np.ndarray]:
        """
        Return strain E_lin + E_nl at each integration point (same order as gauss_data).

        Parameters
        ----------
        U_e : np.ndarray
            Element displacement vector, shape (12,).

        Returns
        -------
        List[np.ndarray]
            One strain vector per Gauss point, each shape (6,) or equivalent.
        """
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
        """
        Return ElementObject with K_e = initial tangent at U=0 (same as linear Timoshenko K_e).

        Caches gauss_data (B, D, shape_functions, shape_derivatives per Gauss point)
        and evaluate_shape_functions for post-processing. B2 shape-function coefficients
        are set (linear Lagrange in ξ, same as linear Timoshenko) for save/load evaluation; see RESULTS_DESIGN.md.

        Returns
        -------
        ElementObject
            K_e = tangent_stiffness_matrix(0), gauss_data, integration_scheme, evaluate_shape_functions.
        """
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData
        self._assert_logging_ready()
        K_e = self.tangent_stiffness_matrix(np.zeros(12, dtype=np.float64))
        xi, w = self.integration_points
        D = self.material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        if self.logger_operator:
            self.logger_operator.log_text(
                "stiffness",
                f"\n=== Element {self.element_id} Stiffness (initial tangent) ===",
            )
            self.logger_operator.log_matrix("stiffness", np.array([[self.L]]), {"name": "Element length L"})
            self.logger_operator.log_matrix("stiffness", D, {"name": "Material matrix D"})
        gauss_cache = []
        for g, (xi_g, w_g) in enumerate(zip(xi, w)):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            Ke_contribution = B.T @ D @ B * w_g * detJ
            if self.logger_operator:
                self._log_gauss_point_stiffness(
                    g, float(xi_g), float(w_g),
                    dN_dξ[0], d2N_dξ2[0], B, Ke_contribution,
                )
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
            self.logger_operator.log_matrix("stiffness", K_e, {"name": "Final K_e (initial tangent)"})
            self.logger_operator.flush("stiffness")
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
        """Build B2 coefficient arrays (12, 6, 4) for N, dN/dξ, d²N/dξ² in monomial basis (linear Lagrange; same as linear Timoshenko)."""
        c = np.zeros((12, 6, 4), dtype=np.float64)
        dc = np.zeros((12, 6, 4), dtype=np.float64)
        d2c = np.zeros((12, 6, 4), dtype=np.float64)
        # N = 0.5(1±ξ) -> coeffs [0.5, ±0.5]; dN/dξ = ±0.5 (constant); d²N/dξ² = 0
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
        """
        Compute the element force vector (external loads): distributed and point loads.

        Returns
        -------
        ForceObject
            F_e and gauss_data; same convention as linear Timoshenko.

        Notes
        -----
        Combines distributed load ``F_dist += sum_g w_g * N.T @ q * detJ`` and point loads
        F_point = N(x_p)ᵀ P at load locations.
        """
        from pre_processing.element_library.gauss_point_data import ForceObject, ForceGaussPointData
        from pre_processing.element_library.linear.timoshenko.utilities.interpolate_loads import LoadInterpolationOperator
        self._assert_logging_ready()
        if self.logger_operator:
            self.logger_operator.log_text(
                "force",
                f"\n=== Element {self.element_id} Force Vector Computation ===",
            )
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
        """Log detailed stiffness integration data for one Gauss point (same format as linear elements)."""
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

    def _log_distributed_loads(
        self,
        xi: np.ndarray,
        weights: np.ndarray,
        N: np.ndarray,
        q: np.ndarray,
        Fe: np.ndarray,
    ) -> None:
        """Log distributed load integration details (same format as linear elements)."""
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
        """Log point load application (same format as linear elements)."""
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
