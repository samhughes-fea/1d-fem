# pre_processing/element_library/nonlinear/timoshenko/nonlinear_timoshenko_3D.py
"""
2-node 3D Timoshenko beam, Total Lagrangian geometric nonlinearity.

**Tensors:** ``U_e`` (12,) or **(14,)** when ``[warping]`` allocates χ DOFs; ``B`` (6, 12) or **(7, 14)**;
``D`` (6, 6) or **(7, 7)** with ``D[6,6]=E·Γ_eff``. ``E = E_lin + E_nl`` on the first six rows;
the warping strain row is **linear** in χ (no Green–Lagrange row 6), matching ``LinearTimoshenkoBeamElement3D`` + warping.
``F_int += B_tot.T @ S * w_g * detJ`` with ``B_tot = B_lin + B_nl`` (``B_nl`` padded to 7×14, seventh row zero).
``K_T = K_0 + K_delta + K_sigma`` — **12 DOF:** ``K_0`` from ``assemble_timoshenko_K0``; **14 DOF:** ``K_0`` from full Gauss on ``B_7x14`` (same as linear Timoshenko + warping).
``K_delta`` is ``Σ (B_totᵀ D B_tot − B_linᵀ D B_lin)``; ``K_sigma`` embeds 12×12 geometric stiffness on the first 12 DOFs.

**Weak forms (Gauss, xi in [-1, 1]):** TL loops use ``quadrature_order`` (``loop_order`` from element slice when unset).

**See Also:** ``LinearTimoshenkoBeamElement3D`` + ``[warping]``; ``NonlinearEulerBernoulliBeamElement3D`` 14-DOF policy.
"""

import logging
from typing import List, Tuple

import numpy as np

from pre_processing.element_library.beam_warping import (
    beam_warping_policy,
    enforce_strict_section_gamma,
    mesh_uses_warping_dof,
    section_gamma_from_section_array,
    warn_if_degenerate_warping_stiffness,
)
from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.nonlinear.timoshenko.tl_green_lagrange_voigt import (
    chord_frame_green_lagrange_voigt_timoshenko_12,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.B_matrix import (
    StrainDisplacementOperator,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.D_matrix import (
    MaterialStiffnessOperator,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.interpolate_loads import (
    LoadInterpolationOperator,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.k0_timoshenko import (
    assemble_timoshenko_K0,
    timoshenko_quadrature_orders_from_element_array,
)
from pre_processing.element_library.nonlinear.timoshenko.utilities.geometric_stiffness import GeometricStiffnessOperator
from pre_processing.element_library.nonlinear.timoshenko.utilities.green_lagrange_strain import GreenLagrangeStrainOperator
from pre_processing.element_library.nonlinear.timoshenko.utilities.stress_resultant import StressResultantOperator
from pre_processing.element_library.shape_function_registry import get_shape_function_operator

logger = logging.getLogger(__name__)

_NL_TS_W_STD = 12
_NL_TS_W_STRAIN = 7
_NL_TS_W_DOF = 14


def _B_nl_6x12_to_7x14(B_nl_6_12: np.ndarray) -> np.ndarray:
    """Pad TL nonlinear ``B`` from (6, 12) to (7, 14); warping row/cols remain zero."""
    B = np.zeros((_NL_TS_W_STRAIN, _NL_TS_W_DOF), dtype=np.float64)
    B[:6, :_NL_TS_W_STD] = B_nl_6_12
    return B


class NonlinearTimoshenkoBeamElement3D(Element1DBase):
    """
    TL Timoshenko: ``K_0`` (linear selective assembly); ``E_nl`` and ``B_tot`` for ``F_int`` and ``K_delta``;
    ``K_sigma`` from ``N``, ``M_y``, ``M_z``.

    **14 local DOFs** when the mesh uses warping (``[warping]`` / ``mesh_uses_warping_dof``): same ``B_7x14`` / ``D_7x7`` pattern as
    ``LinearTimoshenkoBeamElement3D``; Green–Lagrange only on the first six strain rows; ``K_sigma`` embeds 12×12 on the standard beam DOFs.

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
        quadrature_order: int | None = 3,
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
            If ``3`` or ``None``, TL loops use ``TimoshenkoQuadratureOrders.loop_order`` from
            ``element_array``. Otherwise this integer sets the TL loop order (still ≥ 2).

        Notes
        -----
        x_start, x_end, x_global_start, and x_global_end are set from node_coords
        and grid_dictionary (same convention as linear elements).
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
        self._mesh_element_dictionary = element_dictionary
        self._mesh_grid_dictionary = grid_dictionary
        self._mesh_material_dictionary = material_dictionary
        self._mesh_section_dictionary = section_dictionary
        self._timoshenko_orders = timoshenko_quadrature_orders_from_element_array(self.element_array)
        if quadrature_order is None or quadrature_order == 3:
            self.quadrature_order = self._timoshenko_orders.loop_order
        else:
            self.quadrature_order = max(int(quadrature_order), 2)
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
            y_sc=self.y_sc,
            z_sc=self.z_sc,
        )
        self.green_lagrange_strain_operator = GreenLagrangeStrainOperator(
            element_length=self.L,
            include_shear=True,
        )
        self.stress_resultant_operator = StressResultantOperator()
        self.geometric_stiffness_operator = GeometricStiffnessOperator(element_length=self.L)

        self._K_0: np.ndarray | None = None

        if self._n_dof == _NL_TS_W_DOF:
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
        elif self.section_array.size not in (5, 7, 9, 10):
            raise ValueError("Material/section arrays not properly initialised")
        conn = tuple(self.element_array[1:3])
        logger.debug(
            "Element %s geometry initialised\n"
            "• Connectivity: %s\n• Length: %.4e\n• Start/End X: %.4e / %.4e",
            self.element_id, conn, self.L, self.x_start, self.x_end,
        )

    # Property definitions -----------------------------------------------------
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
    def y_sc(self) -> float:
        """Shear centre offset y [m]; 0 if not in section_array."""
        if self.section_array.size >= 9:
            return float(self.section_array[7])
        return 0.0

    @property
    def z_sc(self) -> float:
        """Shear centre offset z [m]; 0 if not in section_array."""
        if self.section_array.size >= 9:
            return float(self.section_array[8])
        return 0.0

    @property
    def Gamma(self) -> float:
        """Warping constant Γ [m⁶] from ``section_array`` index 9 when present."""
        if self.section_array.size >= 10:
            return float(self.section_array[9])
        return 0.0

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

    @property
    def dxi_dx(self) -> float:
        """Natural to physical derivative :math:`\\mathrm{d}\\xi/\\mathrm{d}x = 2/L`."""
        return 2.0 / self.L

    # Operator-compatible formulation methods ----------------------------------
    def shape_functions(self, xi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Shape functions and natural derivatives for assembly (same contract as linear Timoshenko)."""
        return self.shape_function_operator.natural_coordinate_form(xi)

    def B_matrix(self, dN_dξ: np.ndarray, d2N_dξ2: np.ndarray, N: np.ndarray | None = None) -> np.ndarray:
        """``B`` in natural ``ξ`` (same API as linear); TL assembly uses ``physical_coordinate_form`` internally."""
        return self.strain_displacement_operator.natural_coordinate_form(dN_dξ, d2N_dξ2, N)

    def D_matrix(self) -> np.ndarray:
        """Material stiffness ``D`` (6, 6) for Timoshenko section law."""
        return self.material_stiffness_operator.assembly_form()

    def _B_warping_row_ts(self) -> np.ndarray:
        """Row 6 of B: φ_x′ = dθ_x/dx + dχ/dx (same as ``LinearTimoshenkoBeamElement3D``)."""
        row = np.zeros(_NL_TS_W_DOF, dtype=np.float64)
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
        B = np.zeros((n_gauss, _NL_TS_W_STRAIN, _NL_TS_W_DOF), dtype=np.float64)
        B[:, :6, :_NL_TS_W_STD] = B_6x12
        wr = self._B_warping_row_ts()
        for g in range(n_gauss):
            B[g, 6, :] = wr
        return B

    def _D_matrix_7x7(self) -> np.ndarray:
        D6 = self.material_stiffness_operator.assembly_form()
        D = np.zeros((_NL_TS_W_STRAIN, _NL_TS_W_STRAIN), dtype=np.float64)
        D[:6, :6] = D6
        g_eff = self._warp_policy.gamma_effective
        D[6, 6] = self.E * g_eff
        return D

    def _N_14x6_at_xi(self, xi_g: float, N12: np.ndarray) -> np.ndarray:
        Nf = np.zeros((_NL_TS_W_DOF, 6), dtype=np.float64)
        Nf[:_NL_TS_W_STD, :] = N12
        xi = float(xi_g)
        Nf[12, 3] = 0.5 * (1.0 - xi)
        Nf[13, 3] = 0.5 * (1.0 + xi)
        return Nf

    def _tl_voigt_strain_at_gauss(self, U_e: np.ndarray, xi_g: float) -> np.ndarray:
        """
        Six-component Voigt strain ``E`` at ``xi_g`` for the 12-DOF TL Timoshenko path.

        Subclasses (e.g. GESDB) may override for alternate strain definitions; ``internal_force_vector``
        and ``tangent_stiffness_matrix`` consume this through :meth:`_gp_tl_kinematics`.
        """
        return chord_frame_green_lagrange_voigt_timoshenko_12(self, U_e, xi_g)

    def _gp_tl_kinematics(
        self, U_e: np.ndarray, xi_g: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        One Gauss point: shapes, ``B``, total Lagrangian strain ``E``, and ``dN/dx`` first row (12, 6).

        Returns
        -------
        N, dN_dξ, d2N_dξ2, B, E, dN_dx_row
        """
        U_use = np.asarray(U_e, dtype=np.float64).ravel()
        if U_use.size >= _NL_TS_W_DOF:
            U_use = U_use[:_NL_TS_W_STD]
        N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
        B = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
        E = self._tl_voigt_strain_at_gauss(U_use, xi_g)
        dN_dx = dN_dξ[0] * self.dxi_dx
        return N, dN_dξ, d2N_dξ2, B, E, dN_dx

    # K_0 / tangent / nonlinear stress path -----------------------------------
    def _get_K_0(self) -> np.ndarray:
        """Material stiffness ``K_0``: selective Timoshenko (12) or full Gauss (14, warping)."""
        if self._K_0 is not None:
            return self._K_0
        if self._n_dof == _NL_TS_W_DOF:
            D = self._D_matrix_7x7()
            detJ = self.jacobian_determinant
            xi, w = np.polynomial.legendre.leggauss(self.quadrature_order)
            Ke = np.zeros((_NL_TS_W_DOF, _NL_TS_W_DOF), dtype=np.float64)
            for xi_g, w_g in zip(xi, w):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
                Ke += B.T @ D @ B * w_g * detJ
            self._K_0 = Ke
            return self._K_0
        D = self.material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        self._K_0 = assemble_timoshenko_K0(
            self.shape_function_operator,
            self.strain_displacement_operator,
            D,
            detJ,
            self._timoshenko_orders,
        )
        return self._K_0

    def linear_geometric_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """
        Geometric stiffness for **linear** buckling theory about displacement ``U_e``, using the same
        Gauss assembly as :class:`LinearTimoshenkoBeamElement3D` (small-strain ``B``, ``N``, ``M_y``, ``M_z``).

        Required when ``modal.buckling_prestress`` is ``nonlinear_static``: prestress comes from this TL
        element, but :math:`K_\\sigma` for the eigenproblem follows the linearised beam-column form.
        """
        from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
            LinearTimoshenkoBeamElement3D,
        )

        ed = dict(self._mesh_element_dictionary)
        n_t = len(ed["types"])
        ed["types"] = np.array(["LinearTimoshenkoBeamElement3D"] * n_t, dtype=object)
        lin = LinearTimoshenkoBeamElement3D(
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
        Tangent stiffness ``K_T = K_0 + K_delta(U_e) + K_sigma(U_e)``.

        ``K_0`` is the selective linear Timoshenko assembly; ``K_delta`` is the TL material
        stiffness increment ``Σ (B_totᵀ D B_tot − B_linᵀ D B_lin)`` on the TL Gauss loop so that
        at ``U_e = 0``, ``K_delta = 0`` and ``K_T`` matches the linear element stiffness.

        Parameters
        ----------
        U_e : np.ndarray
            Element displacement vector, shape (12,) or (14,) with warping χ DOFs.

        Returns
        -------
        np.ndarray
            Tangent stiffness matrix, shape (12, 12) or (14, 14).
        """
        if self._n_dof == _NL_TS_W_DOF:
            return self._tangent_stiffness_matrix_warping(U_e)

        K_0 = self._get_K_0()
        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        n_g = len(xi)
        N_gp = np.zeros(n_g, dtype=np.float64)
        M_y_gp = np.zeros(n_g, dtype=np.float64)
        M_z_gp = np.zeros(n_g, dtype=np.float64)
        dN_dx_list = []
        for k, (xi_g, w_g) in enumerate(zip(xi, w)):
            _, _, _, _, E, dN_dx_row = self._gp_tl_kinematics(U_e, xi_g)
            N_i, M_y_i, M_z_i = self.stress_resultant_operator.section_forces_from_strain(
                E - self._E_0_voigt, D
            )
            N_gp[k] = N_i
            M_y_gp[k] = M_y_i
            M_z_gp[k] = M_z_i
            dN_dx_list.append(dN_dx_row)
        dN_dx_arr = np.stack(dN_dx_list, axis=0)
        K_sigma = self.geometric_stiffness_operator.assemble_K_sigma(
            N_gp, M_y_gp, M_z_gp, w, dN_dx_arr, detJ
        )
        K_delta = np.zeros((12, 12), dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B_lin = self.strain_displacement_operator.physical_coordinate_form(dN_dξ, d2N_dξ2, N)[0]
            dN_dx = dN_dξ[0] * self.dxi_dx
            d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
            B_nl = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx, d2N_dx2, U_e, N[0]
            )
            B_tot = B_lin + B_nl
            K_delta += (B_tot.T @ D @ B_tot - B_lin.T @ D @ B_lin) * w_g * detJ
        return K_0 + K_delta + K_sigma

    def _tangent_stiffness_matrix_warping(self, U_e: np.ndarray) -> np.ndarray:
        """``K_T`` for 14 DOFs: ``K_0`` + ``K_delta`` (7×14) + embedded ``K_sigma`` (12×12 block)."""
        U_e = np.asarray(U_e, dtype=np.float64).reshape(_NL_TS_W_DOF)
        U12 = U_e[:_NL_TS_W_STD]
        K_0 = self._get_K_0()
        D7 = self._D_matrix_7x7()
        D6 = D7[:6, :6]
        xi, w = self.integration_points
        detJ = self.jacobian_determinant

        n_g = len(xi)
        N_gp = np.zeros(n_g, dtype=np.float64)
        M_y_gp = np.zeros(n_g, dtype=np.float64)
        M_z_gp = np.zeros(n_g, dtype=np.float64)
        dN_dx_list = []
        for k, (xi_g, _w_g) in enumerate(zip(xi, w)):
            _, _, _, _, E6, dN_dx_row = self._gp_tl_kinematics(U12, xi_g)
            N_i, M_y_i, M_z_i = self.stress_resultant_operator.section_forces_from_strain(
                E6 - self._E_0_voigt, D6
            )
            N_gp[k] = N_i
            M_y_gp[k] = M_y_i
            M_z_gp[k] = M_z_i
            dN_dx_list.append(dN_dx_row)
        dN_dx_arr = np.stack(dN_dx_list, axis=0)
        K_sigma_12 = self.geometric_stiffness_operator.assemble_K_sigma(
            N_gp, M_y_gp, M_z_gp, w, dN_dx_arr, detJ
        )
        K_sigma = np.zeros((_NL_TS_W_DOF, _NL_TS_W_DOF), dtype=np.float64)
        K_sigma[:_NL_TS_W_STD, :_NL_TS_W_STD] = K_sigma_12

        K_delta = np.zeros((_NL_TS_W_DOF, _NL_TS_W_DOF), dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
            B_lin_7 = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
            dN_dx = dN_dξ[0] * self.dxi_dx
            d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
            B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx, d2N_dx2, U12, N[0]
            )
            B_nl_7 = _B_nl_6x12_to_7x14(B_nl_12)
            B_tot = B_lin_7 + B_nl_7
            K_delta += (B_tot.T @ D7 @ B_tot - B_lin_7.T @ D7 @ B_lin_7) * w_g * detJ
        return K_0 + K_delta + K_sigma

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        """
        Internal force using ``B_tot = B_lin + B_nl`` with ``S = D @ E``, ``E = E_lin + E_nl``.

        Matches the nonlinear material tangent correction ``K_delta`` (same ``B_tot`` structure).
        """
        if self._n_dof == _NL_TS_W_DOF:
            return self._internal_force_vector_warping(U_e)

        D = self.material_stiffness_operator.assembly_form()
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        F_int = np.zeros(12, dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2, B_lin, E, dN_dx = self._gp_tl_kinematics(U_e, xi_g)
            d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
            B_nl = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx, d2N_dx2, U_e, N[0]
            )
            B_tot = B_lin + B_nl
            S = D @ (E - self._E_0_voigt)
            F_int += (B_tot.T @ S) * w_g * detJ
        return F_int

    def _internal_force_vector_warping(self, U_e: np.ndarray) -> np.ndarray:
        U_e = np.asarray(U_e, dtype=np.float64).reshape(_NL_TS_W_DOF)
        U12 = U_e[:_NL_TS_W_STD]
        D7 = self._D_matrix_7x7()
        E0_7 = np.concatenate([self._E_0_voigt, np.zeros(1, dtype=np.float64)])
        xi, w = self.integration_points
        detJ = self.jacobian_determinant
        F_int = np.zeros(_NL_TS_W_DOF, dtype=np.float64)
        for xi_g, w_g in zip(xi, w):
            N, dN_dξ, d2N_dξ2, B_lin, E6, dN_dx = self._gp_tl_kinematics(U12, xi_g)
            d2N_dx2 = d2N_dξ2[0] * (self.dxi_dx**2)
            B_nl_12 = self.green_lagrange_strain_operator.nonlinear_strain_displacement_gradient(
                dN_dx, d2N_dx2, U12, N[0]
            )
            B_lin_7 = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
            B_nl_7 = _B_nl_6x12_to_7x14(B_nl_12)
            B_tot = B_lin_7 + B_nl_7
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
            Element displacement vector, shape (12,) or (14,) with warping.

        Returns
        -------
        List[np.ndarray]
            One strain vector per Gauss point, shape (6,) or (7,) when warping.
        """
        U = np.asarray(U_e, dtype=np.float64).ravel()
        xi, w = self.integration_points
        result: List[np.ndarray] = []
        if U.size == _NL_TS_W_DOF:
            for xi_g, w_g in zip(xi, w):
                N, dN_dξ, d2N_dξ2, _, E6, _ = self._gp_tl_kinematics(U[:_NL_TS_W_STD], xi_g)
                B_lin_7 = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
                E_warp = float(np.dot(B_lin_7[6], U))
                E = np.concatenate([E6, np.array([E_warp], dtype=np.float64)])
                result.append(E)
            return result
        for xi_g, w_g in zip(xi, w):
            _, _, _, _, E, _ = self._gp_tl_kinematics(U_e, xi_g)
            result.append(np.asarray(E, dtype=np.float64))
        return result

    # Initial tangent / ElementObject cache ------------------------------------
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
        if self._n_dof == _NL_TS_W_DOF:
            U_zero = np.zeros(_NL_TS_W_DOF, dtype=np.float64)
            K_e = self.tangent_stiffness_matrix(U_zero)
            xi, w = self.integration_points
            D = self._D_matrix_7x7()
            detJ = self.jacobian_determinant
            if self.logger_operator:
                self.log_text(
                    "stiffness",
                    f"\n=== Element {self.element_id} Stiffness (initial tangent 14×14) ===",
                )
                self.log_matrix("stiffness", D, {"name": "Material matrix D (7×7)"})
            gauss_cache = []
            for g, (xi_g, w_g) in enumerate(zip(xi, w)):
                N, dN_dξ, d2N_dξ2 = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                B = self._B_matrix_7x14(dN_dξ, d2N_dξ2, N)[0]
                Ke_contribution = B.T @ D @ B * w_g * detJ
                gauss_cache.append(
                    StiffnessGaussPointData(
                        xi=float(xi_g),
                        weight=float(w_g),
                        B_matrix=B.copy(),
                        D_matrix=D.copy(),
                        jacobian=float(detJ),
                        shape_functions=None,
                        shape_derivatives=None,
                    )
                )
                if self.logger_operator:
                    self.log_text("stiffness", f"\nGP {g + 1}: ξ = {xi_g:.6f}, w = {w_g:.6e}")
                    self.log_matrix("stiffness", B, {"name": f"B (7,14) {B.shape}"})
                    self.log_matrix("stiffness", Ke_contribution, {"name": "BᵀDB"})
            if self.logger_operator:
                self.log_matrix("stiffness", K_e, {"name": "Final K_e (initial tangent 14×14)"})
                self.flush_logs("stiffness")
            op = self.shape_function_operator
            evaluate_shape_functions = lambda xi_val: op.natural_coordinate_form(np.asarray(xi_val))
            N_coeffs, dN_coeffs, d2N_coeffs = op.natural_coordinate_form_coefficients()
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

        K_e = self.tangent_stiffness_matrix(np.zeros(12, dtype=np.float64))
        xi, w = self.integration_points
        D = self.material_stiffness_operator.assembly_form()
        detJ = self.jacobian_determinant
        if self.logger_operator:
            self.log_text(
                "stiffness",
                f"\n=== Element {self.element_id} Stiffness (initial tangent) ===",
            )
            self.log_matrix("stiffness", np.array([[self.L]]), {"name": "Element length L"})
            self.log_matrix("stiffness", D, {"name": "Material matrix D"})
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
            self.log_matrix("stiffness", K_e, {"name": "Final K_e (initial tangent)"})
            self.flush_logs("stiffness")
        op = self.shape_function_operator
        evaluate_shape_functions = lambda xi_val: op.natural_coordinate_form(np.asarray(xi_val))
        N_coeffs, dN_coeffs, d2N_coeffs = op.natural_coordinate_form_coefficients()
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

    # Fe tensor computations ----------------------------------------------------
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
        from pre_processing.element_library.gauss_point_data import ForceObject

        self._assert_logging_ready()
        if self.logger_operator:
            self.log_text(
                "force",
                f"\n=== Element {self.element_id} Force Vector Computation ===",
            )
        Fe = np.zeros(self._n_dof, dtype=np.float64)
        gauss_cache = []
        if self.distributed_load_array.size > 0:
            Fe_dist, dist_gauss_cache = self._compute_distributed_load_contribution()
            Fe[:_NL_TS_W_STD] += Fe_dist
            gauss_cache = dist_gauss_cache
        point_loads_cache = None
        if self.point_load_array.size > 0:
            Fe[:_NL_TS_W_STD] += self._compute_point_load_contribution()
            point_loads_cache = self.point_load_array.copy()
        if self.logger_operator:
            self.log_matrix("force", Fe.reshape(1, -1), {"name": "Final Force Vector"})
            self.flush_logs("force")
        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=Fe,
            gauss_data=gauss_cache,
            point_loads=point_loads_cache,
        )

    def _compute_distributed_load_contribution(self):
        """Distributed load contribution using Gauss quadrature (same pattern as linear Timoshenko)."""
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
                n_gauss_points=self.quadrature_order,
            )
            q_gauss = interpolator.interpolate(x_gauss)
            N = np.stack([
                self.shape_function_operator.natural_coordinate_form(np.array([xi]))[0][0]
                for xi in xi_gauss
            ])
            Fe_dist = np.einsum("gij,gj,g->i", N, q_gauss, weights) * (self.L / 2)
            for g, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
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
                self._log_distributed_loads(xi_gauss, weights, N, q_gauss, Fe_dist)
        except Exception as e:
            logger.error("Distributed load error: %s", str(e))
            raise
        return Fe_dist, gauss_cache

    def _compute_point_load_contribution(self) -> np.ndarray:
        """Point load contributions (half-open element interval matches linear Timoshenko)."""
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
            if self.logger_operator:
                self._log_point_load(x_p, xi_p, F_p, N_p, Fe_trans, Fe_rot)
        return Fe_point

    def element_mass_matrix(self):
        """Reference-configuration consistent mass (same as linear Timoshenko) for modal/dynamic."""
        from pre_processing.element_library.gauss_point_data import MassObject

        self._assert_logging_ready()
        if self._n_dof == _NL_TS_W_DOF:
            rho = float(self.material_array[3])
            mu = np.zeros(_NL_TS_W_DOF, dtype=np.float64)
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
            M_e = np.zeros((_NL_TS_W_DOF, _NL_TS_W_DOF), dtype=np.float64)
            xi, w = self.integration_points
            detJ = self.jacobian_determinant
            for xi_g, w_g in zip(xi, w):
                N12, _, _ = self.shape_function_operator.natural_coordinate_form(np.array([xi_g]))
                Ng = self._N_14x6_at_xi(xi_g, N12[0])
                for i in range(_NL_TS_W_DOF):
                    for j in range(_NL_TS_W_DOF):
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

    # Stiffness logging helpers -------------------------------------------------
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
        """Log detailed stiffness-integration data for one Gauss point (same format as linear elements)."""
        if dN_dξ.shape[-2:] != (12, 6):
            raise ValueError(f"dN_dξ shape mismatch: {dN_dξ.shape} ≠ (12, 6)")
        if d2N_dξ2.shape[-2:] != (12, 6):
            raise ValueError(f"d2N_dξ2 shape mismatch: {d2N_dξ2.shape} ≠ (12, 6)")
        if B.shape[-2:] != (6, 12):
            raise ValueError(f"B-matrix shape mismatch: {B.shape} ≠ (*, 6, 12)")
        if contribution.shape != (12, 12):
            raise ValueError(f"Contribution shape mismatch: {contribution.shape} ≠ (12, 12)")
        metadata = {"name": f"GP {gp_idx + 1}", "precision": 6, "max_line_width": 120}
        self.log_text(
            "stiffness",
            f"\nGP {gp_idx + 1}/{self.quadrature_order}: "
            f"ξ = {xi:.6f},  x = {self._xi_to_x(xi):.6e},  w = {weight:.6e}",
        )
        self.log_matrix(
            "stiffness", dN_dξ, {**metadata, "name": f"Shape-function derivative  dN/dξ  {dN_dξ.shape}"}
        )
        self.log_matrix(
            "stiffness", d2N_dξ2, {**metadata, "name": f"Second derivative  d²N/dξ²  {d2N_dξ2.shape}"}
        )
        self.log_matrix(
            "stiffness", B, {**metadata, "name": f"Strain-displacement matrix  B  {B.shape}"}
        )
        self.log_matrix(
            "stiffness", contribution, {**metadata, "name": f"Gauss-point contribution  BᵀDB  {contribution.shape}"}
        )

    # Force logging helpers -----------------------------------------------------
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
        self.log_text("force", "\n=== Distributed Loads ===")
        for gp, (xi_g, w_g) in enumerate(zip(xi, weights)):
            gp_meta = {**metadata, "name": f"GP {gp+1}"}
            self.log_matrix("force", N[gp], {**gp_meta, "name": f"N {N[gp].shape}"})
            self.log_matrix("force", q[gp], {**gp_meta, "name": f"q {q[gp].shape}"})
        self.log_matrix("force", Fe.reshape(1, -1), {**metadata, "name": f"Total Fe {Fe.shape}"})

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
        self.log_text(
            "force",
            f"\n=== Point Load @ x={x:.6e} ===\n"
            f"Natural ξ={xi:.6f}, Element Range: {self.x_start:.6e}-{self.x_end:.6e}",
        )
        self.log_matrix("force", F.reshape(-1, 1), {**metadata, "name": "Force Vector [6×1]"})
        self.log_matrix("force", N, {**metadata, "name": f"Shape Functions {N.shape}"})
        self.log_matrix("force", trans.reshape(-1, 1), {**metadata, "name": "Translations [6×1]"})
        self.log_matrix("force", rot.reshape(-1, 1), {**metadata, "name": "Rotations [6×1]"})

    # Utility methods -----------------------------------------------------------
    def _xi_to_x(self, xi: float) -> float:
        """Convert natural coordinate to physical position (same mapping as linear Timoshenko)."""
        return (xi * self.L / 2) + (self.x_start + self.L / 2)

    def _is_point_in_element(self, x: float) -> bool:
        """Half-open [x_start, x_end) for interior nodes; closed at global right boundary."""
        tol = 1e-12 * self.L
        if np.isclose(self.x_end, self.x_global_end):
            return (self.x_start - tol <= x <= self.x_end + tol)
        return (self.x_start - tol <= x < self.x_end)

    def __repr__(self) -> str:
        return (
            f"NonlinearTimoshenkoBeamElement3D(id={self.element_id}, L={self.L:.2e}m, "
            f"E={self.E:.1e}Pa, quad={self.quadrature_order}, dof={self._n_dof})"
        )
