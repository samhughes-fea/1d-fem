# pre_processing/element_library/linear/timoshenko/linear_timoshenko_3D.py
"""
2-node 3D Timoshenko beam (shear-deformable).

**Tensors:** U_e (12,) node-major (u_x, u_y, u_z, θ_x, θ_y, θ_z) per node; K_e (12,12), F_e (12,);
per Gauss point B (6,12), D (6,6), ε (6,), S = D ε (6,) in the standard Voigt order.
See `docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`. detJ = L/2.

**Weak forms (Gauss, xi in [-1, 1]):** ``K_e += B.T @ D @ B * w_g * detJ`` summed over Gauss point sets
(full rule plus bending/shear block replacements — still weak-form quadrature);
``F_dist += w_g * N.T @ q * detJ``; ``F_point = N.T @ P`` at load station; ``M_e`` consistent mass per
``FORMULATION_DOCSTRING_STANDARDS.md``.

**Kinematics:** non-zero shear strains γ_xy, γ_xz; curvatures from rotations per theory. Local x along chord.

**Constitutive:** D includes EA, EI, κGA shear diagonal terms, and GJ_t (see `utilities/D_matrix.py`).

**Quadrature / selective integration:** ``element_stiffness_matrix`` calls ``assemble_timoshenko_K0`` with
``TimoshenkoQuadratureOrders`` resolved from ``element_array`` (axial, bending_y/z, shear_y/z for the full-rule count,
``shear_block`` default ``1`` for the shear-stiffness rows, torsion). Post-processing ``gauss_cache`` still uses the
max-order Gauss points only.

**Public API:** ``element_stiffness_matrix`` → ``ElementObject``; ``element_force_vector`` → ``ForceObject``;
``element_mass_matrix`` → ``MassObject``.
"""

import numpy as np
from typing import Tuple

# Import Element1DBase class
from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.beam_warping import (
    beam_warping_policy,
    enforce_strict_section_gamma,
    maybe_warn_timoshenko_default_kappa,
    mesh_uses_warping_dof,
    section_gamma_from_section_array,
    warn_if_degenerate_warping_stiffness,
)

# Import MaterialStiffnessOperator and StrainDisplacementOperator classes
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.shape_function_registry import get_shape_function_operator

# Import LoadInterpolationOperator class
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.interpolate_loads import LoadInterpolationOperator
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.k0_timoshenko import (
    assemble_timoshenko_K0,
    timoshenko_quadrature_orders_from_element_array,
)

# --- logging ----------------------------------------------
import logging
logger = logging.getLogger(__name__)
# --- element-level logger ----------------------------------
from pre_processing.element_library.base_logger_operator import BaseLoggerOperator

# Warping extension sizes (Vlasov χ at nodes 12–13; strain row 6)
_TS_W_N_STD = 12
_TS_W_N_STRAIN = 7
_TS_W_N_DOF = 14


class LinearTimoshenkoBeamElement3D(Element1DBase):
    """
    2-node straight beam, 6 DOF per node, 12 total; local ``x`` along the chord.

    Notes
    -----
    **B tensor:** ``(6, 12)`` with all six rows active, including shear rows
    ``gamma_xy = d(u_y)/dx - theta_z`` and ``gamma_xz = d(u_z)/dx - theta_y``.
    **D tensor:** ``(6, 6)`` with shear diagonals ``kappa*G*A`` (rows 3 and 4); optional coupling terms
    may appear if shear-centre offsets are set in the constitutive operator.
    **N tensor:** per-Gauss shape-function slice is ``(12, 6)``; shear rows in ``B`` explicitly depend on ``N`` terms.

    Kinematics: ``gamma_xy = d(u_y)/dx - theta_z``, ``gamma_xz = d(u_z)/dx - theta_y`` (``utilities/B_matrix.py``).
    Constitutive: ``D`` uses ``kappa * G * A`` on shear diagonals.
    Quadrature: orders from ``element_array``; shear columns default to 2 if zero; selective bending/shear assembly (module docstring).
    Loads: ``F_dist += w_g * N.T @ q * detJ`` like EB; optional logging.
    """
    
    # Element formulation identifier for tracking in multi-element meshes
    element_type_name = "Timoshenko-3D"
    
    def __init__(self,
                 *, 
                 element_id: int,
                 element_dictionary: dict,
                 grid_dictionary: dict,
                 section_dictionary: dict,
                 material_dictionary: dict,
                 point_load_array: np.ndarray,
                 distributed_load_array: np.ndarray,
                 job_results_dir: str,
                 quadrature_order: int = 3):
        
        """
        Parameters
        ----------
        element_id, element_dictionary, grid_dictionary, section_dictionary, material_dictionary
            Passed to ``Element1DBase``.
        point_load_array, distributed_load_array
            Point and distributed loads (see base class).
        job_results_dir
            Directory for element logs.
        quadrature_order
            If ``3`` or ``None`` (default path), ``self.quadrature_order`` is set from ``element_array``
            integration columns with shear at least 2. If set to another explicit int, that value is used.
        """
        idx = int(np.where(element_dictionary["ids"] == element_id)[0][0])
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

        # Mass / distributed-load quadrature: default from element_array via TimoshenkoQuadratureOrders.loop_order
        self._timoshenko_orders = timoshenko_quadrature_orders_from_element_array(self.element_array)
        if quadrature_order is None or quadrature_order == 3:
            self.quadrature_order = self._timoshenko_orders.loop_order
        else:
            self.quadrature_order = quadrature_order

        # Geometry
        self.node_coords = self.grid_array                         
        self.L = np.linalg.norm(self.node_coords[1] - self.node_coords[0])
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]

        self.x_global_start = grid_dictionary["coordinates"][:, 0].min()
        self.x_global_end = grid_dictionary["coordinates"][:, 0].max()

        self._warp_policy = beam_warping_policy(
            element_dictionary,
            idx,
            etype_str,
            section_gamma_from_section_array(self.section_array),
        )
        self._warp_mesh = self._warp_policy.mesh_allocates_chi_dof
        self._warp_stiff = self._warp_policy.warping_stiffness_on

        # Initialize element properties and validate
        self._validate_element_properties()
        self._assert_logging_ready()
    
        # Initialize operator classes (registry keys on baseline Timoshenko name)
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

        maybe_warn_timoshenko_default_kappa(self.section_array, element_id)

    def _validate_element_properties(self) -> None:
        """Validate critical element properties"""
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4:
            raise ValueError("Material array not properly initialised")
        if self.section_array.size not in (5, 7, 9, 10):
            raise ValueError("Material/section arrays not properly initialised")
    
        # quick access to node-id pair (n1, n2) from the slice array
        conn = tuple(self.element_array[1:3])

        logger.debug(        # ← use the module-level logger
            f"Element {self.element_id} geometry initialised\n"
            f"• Connectivity: {conn}\n"
            f"• Length: {self.L:.4e}\n"
            f"• Start/End X: {self.x_start:.4e} / {self.x_end:.4e}"
        )

    # Property definitions -----------------------------------------------------
    @property
    def A(self) -> float:
        """Cross-sectional area (m²)"""
        return self.section_array[0]
    
    @property
    def I_y(self) -> float:
        """Moment of inertia about y-axis (m⁴)"""
        return self.section_array[2]
    
    @property
    def I_z(self) -> float:
        """Moment of inertia about z-axis (m⁴)"""
        return self.section_array[3]
    
    @property
    def J_t(self) -> float:
        """Torsional constant (m⁴)"""
        return self.section_array[4]

    @property
    def kappa(self) -> float:
        """Shear correction factor κ (from section preprocessing if available, else 5/6)."""
        if self.section_array.size >= 7:
            return float(self.section_array[5])
        return 5.0 / 6.0  # default for rectangular section

    @property
    def y_sc(self) -> float:
        """Shear centre offset from centroid (y) [m]; 0 if not in section_array."""
        if self.section_array.size >= 9:
            return float(self.section_array[7])
        return 0.0

    @property
    def z_sc(self) -> float:
        """Shear centre offset from centroid (z) [m]; 0 if not in section_array."""
        if self.section_array.size >= 9:
            return float(self.section_array[8])
        return 0.0

    @property
    def Gamma(self) -> float:
        """Warping constant Γ [m⁶] when present (index 9)."""
        if self.section_array.size >= 10:
            return float(self.section_array[9])
        return 0.0

    @property
    def E(self) -> float:
        """Young's modulus (Pa)"""
        return self.material_array[0]

    @property
    def G(self) -> float:
        """Shear modulus (Pa)"""
        return self.material_array[1]

    @property
    def jacobian_determinant(self) -> float:
        """Shortcut to strain operator's Jacobian"""
        return self.L/2

    @property
    def integration_points(self) -> Tuple[np.ndarray, np.ndarray]:
        """Gauss quadrature points/weights"""
        return np.polynomial.legendre.leggauss(self.quadrature_order)

    def _B_warping_row_ts(self) -> np.ndarray:
        """Row 6 of B: φ_x′ = dθ_x/dx + dχ/dx (linear θ_x and χ)."""
        row = np.zeros(_TS_W_N_DOF, dtype=np.float64)
        row[3] = -1.0 / self.L
        row[9] = 1.0 / self.L
        row[12] = -1.0 / self.L
        row[13] = 1.0 / self.L
        return row

    def _B_matrix_7x14(
        self, dN_dξ: np.ndarray, d2N_dξ2: np.ndarray, N: np.ndarray
    ) -> np.ndarray:
        n_gauss = dN_dξ.shape[0]
        B_6x12 = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)
        B = np.zeros((n_gauss, _TS_W_N_STRAIN, _TS_W_N_DOF), dtype=np.float64)
        B[:, :6, :_TS_W_N_STD] = B_6x12
        wr = self._B_warping_row_ts()
        for g in range(n_gauss):
            B[g, 6, :] = wr
        return B

    def _D_matrix_7x7(self) -> np.ndarray:
        D6 = self.material_stiffness_operator.assembly_form()
        D = np.zeros((_TS_W_N_STRAIN, _TS_W_N_STRAIN), dtype=np.float64)
        D[:6, :6] = D6
        g_eff = self._warp_policy.gamma_effective
        D[6, 6] = self.E * g_eff
        return D

    def _N_14x6_at_xi(self, xi_g: float, N12: np.ndarray) -> np.ndarray:
        Nf = np.zeros((_TS_W_N_DOF, 6), dtype=np.float64)
        Nf[:_TS_W_N_STD, :] = N12
        xi = float(xi_g)
        Nf[12, 3] = 0.5 * (1.0 - xi)
        Nf[13, 3] = 0.5 * (1.0 + xi)
        return Nf

    def _precurvature_equivalent_load(self) -> np.ndarray:
        """``sum_g B.T @ D @ E_0 * w_g * detJ`` on the full-rule Gauss loop (matches stiffness cache)."""
        nd = self._n_dof
        if not np.any(self._E_0_voigt):
            return np.zeros(nd, dtype=np.float64)
        detJ = self.jacobian_determinant
        if nd == 12:
            D = self.material_stiffness_operator.assembly_form()
            max_order = self._timoshenko_orders.max_full_order
            xi, w = np.polynomial.legendre.leggauss(max_order)
            f_e = np.zeros(12, dtype=np.float64)
            for xi_g, w_g in zip(xi, w):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
                f_e += (B.T @ D @ self._E_0_voigt) * w_g * detJ
            return f_e
        E0 = np.concatenate([self._E_0_voigt, np.zeros(1, dtype=np.float64)])
        D = self._D_matrix_7x7()
        xi_full, w_full = np.polynomial.legendre.leggauss(self.quadrature_order)
        f_e = np.zeros(_TS_W_N_DOF, dtype=np.float64)
        for xi_g, w_g in zip(xi_full, w_full):
            N_12, dN_dξ_12, d2N_dξ2_12 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self._B_matrix_7x14(dN_dξ_12, d2N_dξ2_12, N_12)[0]
            f_e += (B.T @ D @ E0) * w_g * detJ
        return f_e

    # Operator-compatible formulation methods ----------------------------------
    def shape_functions(self, xi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Shape functions and natural derivatives for Kᵉ assembly."""
        return self.shape_function_operator.natural_coordinate_form(xi)

    def B_matrix(self, dN_dξ: np.ndarray, d2N_dξ2: np.ndarray, N: np.ndarray = None) -> np.ndarray:
        """``B_tilde`` in natural ``xi``; element sums ``K_e += B.T @ D @ B * w_g * detJ`` (physical ``B``)."""
        return self.strain_displacement_operator.natural_coordinate_form(dN_dξ, d2N_dξ2, N)

    def D_matrix(self) -> np.ndarray:
        """D-matrix for stiffness assembly (4×4)."""
        return self.material_stiffness_operator.assembly_form()
    
    # Ke tensor computations ---------------------------------------------------------
    def element_stiffness_matrix(self):
        """
        Compute stiffness matrix using operator classes with integrated logging.
        Uses selective integration: different quadrature orders for bending vs shear terms.
        
        Returns
        -------
        ElementObject
            Element formulation data with cached Gauss point information
        """
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData
        
        self._assert_logging_ready()

        if self._n_dof == 14:
            Ke = np.zeros((_TS_W_N_DOF, _TS_W_N_DOF), dtype=np.float64)
            D = self._D_matrix_7x7()
            xi_full, w_full = np.polynomial.legendre.leggauss(self.quadrature_order)
            gauss_cache = []
            for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
                N_12, dN_dξ_12, d2N_dξ2_12 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self._B_matrix_7x14(dN_dξ_12, d2N_dξ2_12, N_12)[0]
                detJ = self.jacobian_determinant
                Ke += B.T @ D @ B * w_g * detJ
                gauss_cache.append(StiffnessGaussPointData(
                    xi=float(xi_g),
                    weight=float(w_g),
                    B_matrix=B.copy(),
                    D_matrix=D.copy(),
                    jacobian=float(detJ),
                    shape_functions=N_12.copy(),
                    shape_derivatives=dN_dξ_12.copy(),
                ))
            if self.logger_operator:
                self.logger_operator.log_text(
                    "stiffness",
                    f"\n=== Element {self.element_id} Timoshenko (+warp) K_e (14,14) ===",
                )
                self.logger_operator.log_matrix("stiffness", Ke, {"name": "Element Stiffness Matrix"})
                self.logger_operator.flush("stiffness")
            return ElementObject(
                element_id=self.element_id,
                element_type=self.element_type_name,
                K_e=Ke,
                gauss_data=gauss_cache,
                integration_scheme="Gauss-Legendre",
            )

        L = np.array([[self.L]], dtype=np.float64)
        D = self.material_stiffness_operator.assembly_form()
        orders = self._timoshenko_orders
        max_order = orders.max_full_order
        detJ = self.jacobian_determinant

        gauss_cache = []

        if self.logger_operator:  # Modified logging block
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} Stiffness Matrix Computation (Selective Integration) ===")
            self.logger_operator.log_text(
                "stiffness",
                f"Integration orders: axial={orders.axial}, bending_y={orders.bending_y}, bending_z={orders.bending_z}, "
                f"shear_y={orders.shear_y}, shear_z={orders.shear_z}, torsion={orders.torsion}; "
                f"full rule order={max_order}; bending block={orders.bending_rule_order}; shear block={orders.shear_block}",
            )
            self.logger_operator.log_matrix("stiffness", L, {"name": f"Element length  L  {(1,1)}"})
            self.logger_operator.log_matrix("stiffness", D, {"name": f"Material stiffness matrix  D  {D.shape}"})

        xi_full, w_full = np.polynomial.legendre.leggauss(max_order)
        Ke_full_bending = np.zeros((12, 12), dtype=np.float64)
        Ke_full_shear = np.zeros((12, 12), dtype=np.float64)
        Ke_full = np.zeros((12, 12), dtype=np.float64)

        for g, (xi_g, w_g) in enumerate(zip(xi_full, w_full)):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            Ke_contrib = B.T @ D @ B * w_g * detJ
            Ke_full += Ke_contrib
            gauss_cache.append(StiffnessGaussPointData(
                xi=float(xi_g),
                weight=float(w_g),
                B_matrix=B.copy(),
                D_matrix=D.copy(),
                jacobian=float(detJ),
                shape_functions=N.copy(),
                shape_derivatives=dN_dξ.copy()
            ))
            if self.logger_operator:
                B_bending = B[[1, 2], :]
                D_bending_diag = np.array([D[1, 1], D[2, 2]])
                Ke_full_bending += B_bending.T @ np.diag(D_bending_diag) @ B_bending * w_g * detJ
                B_shear = B[[3, 4], :]
                D_shear_diag = np.array([D[3, 3], D[4, 4]])
                Ke_full_shear += B_shear.T @ np.diag(D_shear_diag) @ B_shear * w_g * detJ

        Ke = assemble_timoshenko_K0(
            self.shape_function_operator,
            self.strain_displacement_operator,
            D,
            detJ,
            orders,
        )

        if self.logger_operator:
            Ke_bending_block = np.zeros((12, 12), dtype=np.float64)
            xi_bending, w_bending = np.polynomial.legendre.leggauss(orders.bending_rule_order)
            for xi_g, w_g in zip(xi_bending, w_bending):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
                B_bending = B[[1, 2], :]
                D_bending_diag = np.array([D[1, 1], D[2, 2]])
                Ke_bending_block += B_bending.T @ np.diag(D_bending_diag) @ B_bending * w_g * detJ
            Ke_shear_block = np.zeros((12, 12), dtype=np.float64)
            xi_shear, w_shear = np.polynomial.legendre.leggauss(orders.shear_block)
            for xi_g, w_g in zip(xi_shear, w_shear):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
                B_shear = B[[3, 4], :]
                D_shear_diag = np.array([D[3, 3], D[4, 4]])
                Ke_shear_block += B_shear.T @ np.diag(D_shear_diag) @ B_shear * w_g * detJ
            self.logger_operator.log_matrix("stiffness", Ke_bending_block, {"name": "Bending contribution (selective)"})
            self.logger_operator.log_matrix("stiffness", Ke_shear_block, {"name": "Shear contribution (selective)"})
            self.logger_operator.log_matrix("stiffness", Ke_full - Ke_full_bending - Ke_full_shear, {"name": "Other terms contribution"})
            self.logger_operator.log_matrix("stiffness", Ke, {"name": "Final Element Stiffness Matrix"})
            self.logger_operator.flush("stiffness")

        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=Ke,
            gauss_data=gauss_cache,
            integration_scheme="Gauss-Legendre (Selective)"
        )

    # Fe tensor computations ---------------------------------------------------------
    def element_force_vector(self):
        """Compute the element force vector including distributed and point loads.

        Returns
        -------
        ForceObject
            Element force formulation data with cached Gauss point information

        Notes
        -----
        Combines contributions from:
        - Distributed loads: ``F_dist += sum_g w_g * N.T @ q * detJ``
        - Point loads: F_point = N(x_p)^T P
        """
        from pre_processing.element_library.gauss_point_data import ForceObject, ForceGaussPointData
        
        self._assert_logging_ready()

        Fe = np.zeros(self._n_dof, dtype=np.float64)
        gauss_cache = []
        
        if self.logger_operator:  # Modified logging call
            self.logger_operator.log_text("force", f"\n=== Element {self.element_id} Force Vector Computation ===")

        Fe += self._precurvature_equivalent_load()

        # Process distributed loads
        if self.distributed_load_array.size > 0:
            Fe_dist, dist_gauss_cache = self._compute_distributed_load_contribution()
            Fe[:_TS_W_N_STD] += Fe_dist
            gauss_cache = dist_gauss_cache

        # Process point loads
        point_loads_cache = None
        if self.point_load_array.size > 0:
            Fe[:_TS_W_N_STD] += self._compute_point_load_contribution()
            point_loads_cache = self.point_load_array.copy()

        if self.logger_operator:  # Modified logging block
            self.logger_operator.log_matrix("force", Fe.reshape(1, -1), {"name": "Final Force Vector"})
            self.logger_operator.flush("force")

        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=gauss_cache,
            point_loads=point_loads_cache
        )

    def linear_geometric_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """
        Geometric stiffness from linear strains :math:`\\varepsilon = B U_e`, :math:`S = D \\varepsilon`.
        **12-DOF** standard Timoshenko; **14-DOF** warping: same embedded **12×12** ``K_σ`` on the first
        twelve DOFs as ``NonlinearTimoshenkoBeamElement3D`` with warping (χ rows/cols zero in ``K_σ``).
        """
        U_e = np.asarray(U_e, dtype=np.float64).ravel()
        from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.geometric_stiffness import (
            GeometricStiffnessOperator,
        )
        from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.stress_resultant import (
            StressResultantOperator,
        )

        stress_op = StressResultantOperator()
        geo = GeometricStiffnessOperator(element_length=self.L)
        D6 = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        dξ_dx = 2.0 / self.L
        n_g = len(xi)
        N_gp = np.zeros(n_g, dtype=np.float64)
        M_y_gp = np.zeros(n_g, dtype=np.float64)
        M_z_gp = np.zeros(n_g, dtype=np.float64)
        dN_dx_list = []

        if U_e.size == 12:
            for k, (xi_g, w_g) in enumerate(zip(xi, w)):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
                eps = B @ U_e - self._E_0_voigt
                Ni, Myi, Mzi = stress_op.section_forces_from_strain(eps, D6)
                N_gp[k] = Ni
                M_y_gp[k] = Myi
                M_z_gp[k] = Mzi
                dN_dx_list.append(dN_dξ[0] * dξ_dx)
        elif U_e.size == 14 and self._n_dof == 14:
            E0 = np.concatenate([self._E_0_voigt, np.zeros(1, dtype=np.float64)])
            for k, (xi_g, w_g) in enumerate(zip(xi, w)):
                N12, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N12)[0]
                eps = B @ U_e - E0
                Ni, Myi, Mzi = stress_op.section_forces_from_strain(eps[:6], D6)
                N_gp[k] = Ni
                M_y_gp[k] = Myi
                M_z_gp[k] = Mzi
                dN_dx_list.append(dN_dξ[0] * dξ_dx)
        else:
            raise ValueError(
                f"linear_geometric_stiffness_matrix: expected 12 or (14 with warping element) DOFs, got {U_e.size}"
            )

        dN_dx_arr = np.stack(dN_dx_list, axis=0)
        K12 = geo.assemble_K_sigma(N_gp, M_y_gp, M_z_gp, w, dN_dx_arr, detJ)
        if U_e.size == 12:
            return K12
        K_sigma = np.zeros((_TS_W_N_DOF, _TS_W_N_DOF), dtype=np.float64)
        K_sigma[:_TS_W_N_STD, :_TS_W_N_STD] = K12
        return K_sigma

    def element_mass_matrix(self):
        """Consistent mass: translations rho*A; torsion rho*J_t; bending rotations rho*Iy / rho*Iz."""
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        if self._n_dof == 14:
            rho = float(self.material_array[3])
            mu = np.zeros(_TS_W_N_DOF, dtype=np.float64)
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
            M_e = np.zeros((_TS_W_N_DOF, _TS_W_N_DOF), dtype=np.float64)
            xi, w = self.integration_points
            detJ = self.jacobian_determinant
            for xi_g, w_g in zip(xi, w):
                N12, _, _ = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                Ng = self._N_14x6_at_xi(xi_g, N12[0])
                for i in range(_TS_W_N_DOF):
                    for j in range(_TS_W_N_DOF):
                        mij = 0.5 * (mu[i] + mu[j])
                        M_e[i, j] += mij * float(np.dot(Ng[i, :], Ng[j, :])) * w_g * detJ
            return MassObject(
                element_id=self.element_id,
                element_type=self.element_type_name,
                M_e=M_e,
            )

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



    # Fe - distributed load contribution ------------------------------------------------
    def _compute_distributed_load_contribution(self):
        """Compute distributed load contribution using Gauss quadrature.
        
        Returns
        -------
        tuple[np.ndarray, List[ForceGaussPointData]]
            Force vector and Gauss point cache
        """
        from pre_processing.element_library.gauss_point_data import ForceGaussPointData
        
        xi_gauss, weights = self.integration_points
        x_gauss = (xi_gauss + 1) * (self.L / 2) + self.x_start
        Fe_dist = np.zeros(12, dtype=np.float64)
        gauss_cache = []
        detJ = self.jacobian_determinant

        try:
            interpolator = LoadInterpolationOperator(
                distributed_loads_array=self.distributed_load_array,
                boundary_mode="error",
                interpolation_order="cubic",
                n_gauss_points=self.quadrature_order
            )
            q_gauss = interpolator.interpolate(x_gauss)
            N = np.stack([self.shape_function_operator.natural_coordinate_form(xi)[0][0]
                        for xi in xi_gauss])
            Fe_dist = np.einsum("gij,gj,g->i", N, q_gauss, weights) * (self.L / 2)
            
            # Cache Gauss point data
            for g, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
                gauss_cache.append(ForceGaussPointData(
                    xi=float(xi_g),
                    weight=float(w_g),
                    shape_functions=N[g].copy(),
                    jacobian=float(detJ),
                    distributed_load=q_gauss[g].copy() if g < len(q_gauss) else None
                ))

            if self.logger_operator:
                self._log_distributed_loads(xi_gauss, weights, N, q_gauss, Fe_dist)

        except Exception as e:
            logger.error(f"Distributed load error: {str(e)}")
            raise

        return Fe_dist, gauss_cache
    
    # Fe - point load contribution
    def _compute_point_load_contribution(self) -> np.ndarray:
        """Compute point load contributions using shape function evaluation."""
        from pre_processing.element_library.point_load_utils import add_phased_increment, point_load_phase_rad

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
            inc = np.zeros(12, dtype=np.float64)
            inc[[0, 1, 2, 6, 7, 8]] = Fe_trans
            inc[[3, 4, 5, 9, 10, 11]] = Fe_rot
            Fe_point = add_phased_increment(Fe_point, inc, point_load_phase_rad(load))

            if self.logger_operator:
                self._log_point_load(x_p, xi_p, F_p, N_p, Fe_trans, Fe_rot)

        return Fe_point

    # Stiffness logging helpers ----------------------------------------------------------
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
        """
        Log detailed stiffness-integration data for one Gauss point.

        All tensors are shape-checked before logging to guard against silent
        dimensional errors.
        """
        # --- shape validation ----------------------------------------------------
        if dN_dξ.shape[-2:] != (12, 6):
            raise ValueError(f"dN_dξ shape mismatch: {dN_dξ.shape} ≠ (12, 6)")
        if d2N_dξ2.shape[-2:] != (12, 6):
            raise ValueError(f"d2N_dξ2 shape mismatch: {d2N_dξ2.shape} ≠ (12, 6)")
        if B.shape[-2:] != (6, 12):                     # ← updated (was 4 × 12)
            raise ValueError(f"B-matrix shape mismatch: {B.shape} ≠ (*, 6, 12)")
        if contribution.shape != (12, 12):
            raise ValueError(f"Contribution shape mismatch: {contribution.shape} ≠ (12, 12)")

        # --- metadata for pretty print ------------------------------------------
        metadata = {
            "name": f"GP {gp_idx + 1}",
            "precision": 6,
            "max_line_width": 120,
        }

        # --- header --------------------------------------------------------------
        self.log_text(
            "stiffness",
            f"\nGP {gp_idx + 1}/{self.quadrature_order}: "
            f"ξ = {xi:.6f},  x = {self._xi_to_x(xi):.6e},  w = {weight:.6e}",
        )

        # --- derivative matrices -------------------------------------------------
        self.log_matrix(
            "stiffness", dN_dξ,
            {**metadata, "name": f"Shape-function derivative  dN/dξ  {dN_dξ.shape}"}
        )
        self.log_matrix(
            "stiffness", d2N_dξ2,
            {**metadata, "name": f"Second derivative  d²N/dξ²  {d2N_dξ2.shape}"}
        )

        # --- B-matrix & BᵀDB contribution --------------------------------------------
        self.log_matrix(
            "stiffness", B,
            {**metadata, "name": f"Strain-displacement matrix  B  {B.shape}"}
        )
        self.log_matrix(
            "stiffness", contribution,
            {**metadata, "name": f"Gauss-point contribution  BᵀDB  {contribution.shape}"}
        )

    # Force logging helpers ----------------------------------------------------------
    def _log_distributed_loads(self, xi: np.ndarray, weights: np.ndarray,
                            N: np.ndarray, q: np.ndarray, Fe: np.ndarray):
        """Log distributed load integration details with structural validation."""
        if not self.logger_operator:
            return

        # Validate tensor dimensions
        if N.shape[1:] != (12, 6):
            raise ValueError(f"N shape mismatch: {N.shape} != (n_pts,12,6)")
        if q.shape[1] != 6:
            raise ValueError(f"Load vector shape: {q.shape} != (n_pts,6)")
        if Fe.shape != (12,):  # Fixed shape validation
            raise ValueError(f"Fe result shape: {Fe.shape} != (12,)")

        metadata = {
            "precision": 6,
            "max_line_width": 100
        }
    
        self.logger_operator.log_text("force", "\n=== Distributed Loads ===")
    
        for gp, (xi_g, w_g) in enumerate(zip(xi, weights)):
            gp_meta = {**metadata, "name": f"GP {gp+1}"}
        
            # Validate per-GP matrices
            if N[gp].shape != (12, 6):
                raise ValueError(f"GP {gp} N shape: {N[gp].shape} != (12,6)")
            if q[gp].shape != (6,):
                raise ValueError(f"GP {gp} q shape: {q[gp].shape} != (6,)")

            self.logger_operator.log_matrix("force", N[gp],
                        {**gp_meta, "name": f"N {N[gp].shape}"})
            self.logger_operator.log_matrix("force", q[gp],
                        {**gp_meta, "name": f"q {q[gp].shape}"})

        self.logger_operator.log_matrix("force", Fe.reshape(1, -1),  # Ensure 2D for logging
                                    {**metadata, "name": f"Total Fe {Fe.shape}"})

    def _log_point_load(self, x: float, xi: float, F: np.ndarray,
                    N: np.ndarray, trans: np.ndarray, rot: np.ndarray):
        """Log point load application with mechanical validation."""
        if not self.logger_operator:
            return

        # Validate tensor dimensions
        if N.shape != (12, 6):
            raise ValueError(f"N shape mismatch: {N.shape} != (12,6)")
        if trans.shape != (6,):
            raise ValueError(f"Translation vector shape: {trans.shape} != (6,)")
        if rot.shape != (6,):
            raise ValueError(f"Rotation vector shape: {rot.shape} != (6,)")

        metadata = {
            "precision": 6,
            "max_line_width": 120
        }
    
        self.logger_operator.log_text("force",
            f"\n=== Point Load @ x={x:.6e} ===\n"
            f"Natural ξ={xi:.6f}, Element Range: {self.x_start:.6e}-{self.x_end:.6e}"
        )
    
        # Log with structural context
        self.logger_operator.log_matrix("force", F.reshape(-1, 1),  # Column vector
                                    {**metadata, "name": "Force Vector [6×1]"})
        self.logger_operator.log_matrix("force", N, 
                                    {**metadata, "name": f"Shape Functions {N.shape}"})
        self.logger_operator.log_matrix("force", trans.reshape(-1, 1),
                                    {**metadata, "name": "Translations [6×1]"})
        self.logger_operator.log_matrix("force", rot.reshape(-1, 1),
                                    {**metadata, "name": "Rotations [6×1]"})

    def _xi_to_x(self, xi: float) -> float:
        """Convert natural coordinate to physical position"""
        return (xi * self.L/2) + (self.x_start + self.L/2)

    # Utility methods ----------------------------------------------------------
    def _is_point_in_element(self, x: float) -> bool:
        """Check if point x is within element bounds. Uses half-open [x_start, x_end)
        for non-end elements so a point load at an interior node is assigned to
        exactly one element (the one to the right of the node). Avoids double-counting."""
        tol = 1e-12 * self.L
        if np.isclose(self.x_end, self.x_global_end):
            return (self.x_start - tol <= x <= self.x_end + tol)
        return (self.x_start - tol <= x < self.x_end)

    def __repr__(self) -> str:
        return (f"LinearTimoshenkoBeamElement3D(id={self.element_id}, L={self.L:.2e}m, "
                f"E={self.E:.1e}Pa, quad={self.quadrature_order})")