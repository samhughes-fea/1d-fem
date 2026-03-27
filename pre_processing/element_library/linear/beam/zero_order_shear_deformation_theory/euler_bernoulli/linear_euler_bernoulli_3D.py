# pre_processing/element_library/linear/euler_bernoulli/linear_euler_bernoulli_3D.py
"""
2-node 3D Euler–Bernoulli beam.

**Tensors:** U_e (12,) node-major DOFs (u_x, u_y, u_z, θ_x, θ_y, θ_z) per node;
K_e (12,12), F_e (12,); per Gauss point B (6,12), D (6,6), ε (6,), S = D ε (6,).
Voigt order follows `docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`. detJ = |J| = L/2.

**Weak forms (Gauss, ξ in [-1, 1]):** `K_e += B.T @ D @ B * w_g * detJ` summed over Gauss points;
`F_dist += w_g * N.T @ q * detJ`; `F_point = N.T @ P` at the load station; M_e via consistent mass.

**Kinematics:** γ_xy = γ_xz = 0; κ_y and κ_z from transverse displacement (Hermite); local x along chord.
See `utilities/B_matrix.py`.

**Constitutive:** D is diagonal in axial, bending, torsion; shear rows of D are zero (V_y, V_z from equilibrium, not D ε).

**Quadrature:** Gauss–Legendre order ``quadrature_order`` from argument or ``element_array`` (axial, bending_y, bending_z, torsion, load columns).

**Public API:** ``element_stiffness_matrix`` → ``ElementObject``; ``element_force_vector`` → ``ForceObject``;
``element_mass_matrix`` → ``MassObject`` (consistent mass).
"""

import numpy as np
from typing import Optional, Tuple

# Import Element1DBase class
from pre_processing.element_library.element_1D_base import Element1DBase

# Import MaterialStiffnessOperator and StrainDisplacementOperator classes
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.shape_function_registry import get_shape_function_operator

# Import LoadInterpolationOperator class
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.interpolate_loads import LoadInterpolationOperator

# --- logging ----------------------------------------------
import logging
logger = logging.getLogger(__name__)
# --- element-level logger ----------------------------------
from pre_processing.element_library.base_logger_operator import BaseLoggerOperator

class LinearEulerBernoulliBeamElement3D(Element1DBase):
    """
    2-node straight beam, 6 DOF per node, 12 total (see module docstring for ``U_e`` order).

    Notes
    -----
    **B tensor:** ``(6, 12)`` with active rows ``eps_x``, ``kappa_y``, ``kappa_z``, ``phi_x``;
    shear rows ``gamma_xy`` and ``gamma_xz`` are structural zeros in EB kinematics.
    **D tensor:** ``(6, 6)`` diagonal in axial/bending/torsion with shear rows/cols zero.
    **N tensor:** per-Gauss shape-function slice is ``(12, 6)`` in node-major DOF order used by loads and mass.

    Constitutive: ``D`` is diagonal in axial, bending, and torsion; rows 3–4 (shear) are zero, so
    ``V_y = V_z = 0`` from ``S = D @ eps``. Equilibrium gives shear ``V = dM/dx``, not ``D @ eps``.

    Loads: ``F_dist += w_g * N.T @ q * detJ``; point loads ``N.T @ P`` at the load station. See ``utilities/interpolate_loads.py``.

    Logging: optional per-element logs for stiffness and force matrices.
    """
    
    # Element formulation identifier for tracking in multi-element meshes
    element_type_name = "Euler-Bernoulli-3D"
    
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
                 quadrature_order: Optional[int] = None):
        
        """
        Parameters
        ----------
        element_id, element_dictionary, grid_dictionary, section_dictionary, material_dictionary
            Passed to ``Element1DBase`` (connectivity, coordinates, section and material arrays).
        point_load_array, distributed_load_array
            Loads; columns ``[x,y,z, Fx,Fy,Fz, Mx,My,Mz]`` for distributed table.
        job_results_dir
            Directory for element logs.
        quadrature_order
            Gauss–Legendre order. If ``None``, taken as the max of ``element_array`` fields
            ``[axial, bending_y, bending_z, torsion, load]`` (indices 3–5, 8–9), at least 2.
            See ``docs/element_library/shape_function_conventions.md``.
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

        # Quadrature order: from element_array when not provided (same 7-column convention as Timoshenko)
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

        # Geometry
        self.node_coords = self.grid_array                         
        self.L = np.linalg.norm(self.node_coords[1] - self.node_coords[0])
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]

        self.x_global_start = grid_dictionary["coordinates"][:, 0].min()
        self.x_global_end = grid_dictionary["coordinates"][:, 0].max()

        # Initialize element properties and validate
        self._validate_element_properties()
        self._assert_logging_ready()
    
        # Initialize operator classes (shape-function implementation from registry)
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

    def _validate_element_properties(self) -> None:
        """Validate critical element properties"""
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4 or self.section_array.size not in (5, 7):
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

    # Operator-compatible formulation methods ----------------------------------
    def shape_functions(self, xi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Shape functions and natural derivatives for Kᵉ assembly."""
        return self.shape_function_operator.natural_coordinate_form(xi)

    def B_matrix(self, dN_dξ: np.ndarray, d2N_dξ2: np.ndarray) -> np.ndarray:
        """``B_tilde`` in natural ``xi``; element sums ``K_e += B.T @ D @ B * w_g * detJ`` (physical ``B``)."""
        return self.strain_displacement_operator.natural_coordinate_form(dN_dξ, d2N_dξ2)

    def D_matrix(self) -> np.ndarray:
        """D-matrix for stiffness assembly (4×4)."""
        return self.material_stiffness_operator.assembly_form()

    # Ke tensor computations ---------------------------------------------------------
    def element_stiffness_matrix(self):
        """
        Compute stiffness matrix using operator classes with integrated logging.
        
        Returns
        -------
        ElementObject
            Element formulation data with cached Gauss point information
        """
        from pre_processing.element_library.gauss_point_data import ElementObject, StiffnessGaussPointData
        
        self._assert_logging_ready()

        Ke = np.zeros((12, 12), dtype=np.float64)
        L = np.array([[self.L]], dtype=np.float64)
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        
        gauss_cache = []

        if self.logger_operator:  # Modified logging block
            self.logger_operator.log_text("stiffness", f"\n=== Element {self.element_id} Stiffness Matrix Computation ===")
            self.logger_operator.log_matrix("stiffness", L, {"name": f"Element length  L  {(1,1)}"})
            self.logger_operator.log_matrix("stiffness", D, {"name": f"Material stiffness matrix  D  {D.shape}"})
        
        for g, (xi_g, w_g) in enumerate(zip(xi, w)):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)[0]
            detJ = self.jacobian_determinant
            Ke_contribution = B.T @ D @ B * w_g * detJ # note B is in cartesian form, w_g is not and requires scaling
            Ke += Ke_contribution
            
            # Cache Gauss point data
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
                self._log_gauss_point_stiffness(g, xi_g, w_g, dN_dξ, d2N_dξ2, B, Ke_contribution)

        if self.logger_operator:  # Modified logging block
            self.logger_operator.log_matrix("stiffness", Ke, {"name": "Final Element Stiffness Matrix"})
            self.logger_operator.flush("stiffness")

        # B1: attach callable so post-processing can evaluate N(ξ) at arbitrary ξ
        op = self.shape_function_operator
        evaluate_shape_functions = lambda xi: op.natural_coordinate_form(np.asarray(xi))

        # B2: polynomial coefficients for N, dN/dξ, d²N/dξ² (monomial basis ξ^0..ξ^3), shape (12, 6, 4)
        N_coeffs, dN_coeffs, d2N_coeffs = op.build_shape_function_coefficients_b2()

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

        Fe = np.zeros(12, dtype=np.float64)
        gauss_cache = []
        
        if self.logger_operator:  # Modified logging call
            self.logger_operator.log_text("force", f"\n=== Element {self.element_id} Force Vector Computation ===")

        # Process distributed loads
        if self.distributed_load_array.size > 0:
            Fe_dist, dist_gauss_cache = self._compute_distributed_load_contribution()
            Fe += Fe_dist
            gauss_cache = dist_gauss_cache

        # Process point loads
        point_loads_cache = None
        if self.point_load_array.size > 0:
            Fe += self._compute_point_load_contribution()
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

    def element_mass_matrix(self):
        """Consistent mass: translations rho*A; torsion rho*J_t; bending rotations rho*Iy / rho*Iz."""
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
        Fe_point = np.zeros(12, dtype=np.float64)
        
        for load in self.point_load_array:
            x_p = float(load[0])
            F_p = load[3:9].astype(np.float64)
            
            if not self._is_point_in_element(x_p):
                continue

            xi_p = 2 * (x_p - self.x_start) / self.L - 1
            N_p = self.shape_function_operator.natural_coordinate_form(np.array([xi_p]))[0][0]

            Fe_trans = np.einsum("ij,j->i", N_p[[0,1,2,6,7,8], :3], F_p[:3])
            Fe_rot = np.einsum("ij,j->i", N_p[[3,4,5,9,10,11], 3:], F_p[3:])
            Fe_point[[0,1,2,6,7,8]] += Fe_trans
            Fe_point[[3,4,5,9,10,11]] += Fe_rot

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
        return (f"LinearEulerBernoulliBeamElement3D(id={self.element_id}, L={self.L:.2e}m, "
                f"E={self.E:.1e}Pa, quad={self.quadrature_order})")