# pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_3D.py
"""
Co-rotational 3D beam element (12 DOF).

Small-strain linear Timoshenko or Euler–Bernoulli material stiffness is formed in a **corotated**
basis whose first axis follows the **current chord** between deformed nodes; the result is mapped
to global Cartesian DOFs with a block-orthogonal transform **T**(U_e) built from four 3×3 blocks **R**.T
(two translational and two rotational triplets).

**Internal force:** F_int = T.T @ K_local(L) @ (T @ U_e), with K_local assembled like the linear beam at length L.

**Tangent:** The consistent Jacobian follows from differentiating ``F_int`` w.r.t. ``U_e``. By default
``tangent_stiffness_matrix`` uses **central finite differences** on ``internal_force_vector``
(``O(12)`` force evaluations per call). Alternatively ``tangent_stiffness_mode="elastic_material"``
returns only the symmetric **elastic material** block ``Tᵀ K_local(L) T`` at the current chord; this
**omits spin stiffness** from ∂**T**/∂**U** and ∂**L**/∂**U** (use for diagnostics or modified Newton,
not full equivalence to the FD tangent at finite rotation).

**Linear analysis:** ``element_stiffness_matrix`` uses the reference chord (undeformed geometry).

See ``docs/element_library/large_rotation_vs_total_lagrangian.md`` for contrasts with Total Lagrangian beams.

See Also
--------
pre_processing.element_library.nonlinear.timoshenko.nonlinear_timoshenko_3D
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
"""

from __future__ import annotations

import logging
from typing import Literal, Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.gauss_point_data import ElementObject, ForceObject, MassObject
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.B_matrix import (
    StrainDisplacementOperator as TimoshenkoStrain,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.D_matrix import (
    MaterialStiffnessOperator as TimoshenkoMaterial,
)
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.k0_timoshenko import (
    assemble_timoshenko_K0,
    timoshenko_quadrature_orders_from_element_array,
)
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.B_matrix import (
    StrainDisplacementOperator as EulerBernoulliStrain,
)
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.D_matrix import (
    MaterialStiffnessOperator as EulerBernoulliMaterial,
)
from pre_processing.element_library.shape_function_registry import get_shape_function_operator

logger = logging.getLogger(__name__)


def _orthogonal_basis_from_chord(d: np.ndarray) -> np.ndarray:
    """
    Rotation matrix R (3, 3) with columns e1, e2, e3 = orthonormal basis;
    e1 aligns with chord direction ``d`` (global Cartesian).
    """
    L = float(np.linalg.norm(d))
    if L <= 1e-18:
        raise ValueError("Corotational chord length is degenerate.")
    e1 = d / L
    a = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    if abs(float(np.dot(e1, a))) > 0.95:
        a = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    e2 = np.cross(e1, a)
    e2 = e2 / np.linalg.norm(e2)
    e3 = np.cross(e1, e2)
    return np.column_stack((e1, e2, e3))


def _block_twelve(R: np.ndarray) -> np.ndarray:
    """Map global nodal vectors to corotated axes: four R.T blocks (trans/rot at each node)."""
    Rt = R.T
    return np.kron(np.eye(4, dtype=np.float64), Rt)


class CorotationalBeamElement3D(Element1DBase):
    """
    Co-rotational 3D beam (``kernel`` = ``timoshenko`` or ``euler_bernoulli``).

    Notes
    -----
    Uses linear ``K_\\mathrm{local}`` at ``L`` with beam operators aligned with corotated axis ``e_1``.
    """

    element_type_name = "Corotational-3D"

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
        quadrature_order: int | None = None,
        kernel: str = "timoshenko",
        tangent_stiffness_mode: Literal["finite_difference", "elastic_material"] = "finite_difference",
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

        self._element_dictionary = element_dictionary
        self._grid_dictionary = grid_dictionary
        self._section_dictionary = section_dictionary
        self._material_dictionary = material_dictionary

        if kernel not in ("euler_bernoulli", "timoshenko"):
            raise ValueError(
                f"kernel must be 'euler_bernoulli' or 'timoshenko', got {kernel!r}"
            )
        self.kernel = kernel
        if tangent_stiffness_mode not in ("finite_difference", "elastic_material"):
            raise ValueError(
                "tangent_stiffness_mode must be 'finite_difference' or 'elastic_material', "
                f"got {tangent_stiffness_mode!r}"
            )
        self.tangent_stiffness_mode = tangent_stiffness_mode

        if quadrature_order is not None:
            self.quadrature_order = int(quadrature_order)
        else:
            axial_order = int(self.element_array[3])
            bending_y_order = int(self.element_array[4])
            bending_z_order = int(self.element_array[5])
            shear_y_order = int(self.element_array[6])
            shear_z_order = int(self.element_array[7])
            torsion_order = int(self.element_array[8])
            load_order = int(self.element_array[9])
            self.quadrature_order = max(
                axial_order,
                bending_y_order,
                bending_z_order,
                shear_y_order,
                shear_z_order,
                torsion_order,
                load_order,
                2,
            )

        self.node_coords = self.grid_array
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = float(grid_dictionary["coordinates"][:, 0].min())
        self.x_global_end = float(grid_dictionary["coordinates"][:, 0].max())

        self._validate_element_properties()
        self._assert_logging_ready()

    def _validate_element_properties(self) -> None:
        if self.L <= 0:
            raise ValueError(f"Invalid element length {self.L:.2e} for element {self.element_id}")
        if self.material_array.size != 4 or self.section_array.size not in (5, 7, 9, 10):
            raise ValueError("Material/section arrays not properly initialised")

    @property
    def A(self) -> float:
        return float(self.section_array[0])

    @property
    def I_y(self) -> float:
        return float(self.section_array[2])

    @property
    def I_z(self) -> float:
        return float(self.section_array[3])

    @property
    def J_t(self) -> float:
        return float(self.section_array[4])

    @property
    def kappa(self) -> float:
        return float(self.section_array[5]) if self.section_array.size >= 7 else 5.0 / 6.0

    @property
    def y_sc(self) -> float:
        return float(self.section_array[7]) if self.section_array.size >= 9 else 0.0

    @property
    def z_sc(self) -> float:
        return float(self.section_array[8]) if self.section_array.size >= 9 else 0.0

    @property
    def E(self) -> float:
        return float(self.material_array[0])

    @property
    def G(self) -> float:
        return float(self.material_array[1])

    def _linear_reference_element(self) -> Element1DBase:
        """Undeformed linear beam element for consistent loads and mass (reference chord)."""
        kwargs = dict(
            element_id=self.element_id,
            element_dictionary=self._element_dictionary,
            grid_dictionary=self._grid_dictionary,
            section_dictionary=self._section_dictionary,
            material_dictionary=self._material_dictionary,
            point_load_array=self.point_load_array,
            distributed_load_array=self.distributed_load_array,
            job_results_dir=self.job_results_dir,
            quadrature_order=self.quadrature_order,
        )
        if self.kernel == "timoshenko":
            from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
                LinearTimoshenkoBeamElement3D,
            )

            return LinearTimoshenkoBeamElement3D(**kwargs)
        from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
            LinearEulerBernoulliBeamElement3D,
        )

        return LinearEulerBernoulliBeamElement3D(**kwargs)

    def _timoshenko_K_local(self, L: float) -> np.ndarray:
        orders = timoshenko_quadrature_orders_from_element_array(self.element_array)
        sf = get_shape_function_operator("LinearTimoshenkoBeamElement3D", L)
        strain = TimoshenkoStrain(element_length=L)
        mat = TimoshenkoMaterial(
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
        D = mat.assembly_form()
        detJ = L / 2.0
        return assemble_timoshenko_K0(sf, strain, D, detJ, orders)

    def _euler_bernoulli_K_local(self, L: float) -> np.ndarray:
        axial = max(int(self.element_array[3]), 1)
        bend_y = max(int(self.element_array[4]), 2)
        bend_z = max(int(self.element_array[5]), 2)
        torsion = max(int(self.element_array[8]), 1)
        n_gauss = max(axial, bend_y, bend_z, torsion, 2)
        sf = get_shape_function_operator("LinearEulerBernoulliBeamElement3D", L)
        strain = EulerBernoulliStrain(element_length=L)
        mat = EulerBernoulliMaterial(
            youngs_modulus=self.E,
            shear_modulus=self.G,
            cross_section_area=self.A,
            moment_inertia_y=self.I_y,
            moment_inertia_z=self.I_z,
            torsion_constant=self.J_t,
        )
        D = mat.assembly_form()
        xi, w = np.polynomial.legendre.leggauss(n_gauss)
        Ke = np.zeros((12, 12), dtype=np.float64)
        detJ = L / 2.0
        for xi_g, w_g in zip(xi, w):
            _, dN_dξ, d2N_dξ2 = sf.natural_coordinate_form(np.array([xi_g]))
            B = strain.physical_coordinate_form(dN_dξ, d2N_dξ2)[0]
            Ke += B.T @ D @ B * w_g * detJ
        return Ke

    def _K_local(self, L: float) -> np.ndarray:
        if self.kernel == "timoshenko":
            return self._timoshenko_K_local(L)
        return self._euler_bernoulli_K_local(L)

    def _Tm_reference(self) -> np.ndarray:
        d0 = self.node_coords[1] - self.node_coords[0]
        R = _orthogonal_basis_from_chord(d0)
        return _block_twelve(R)

    def _Tm_current(self, U_e: np.ndarray) -> Tuple[np.ndarray, float]:
        u = np.asarray(U_e, dtype=np.float64).reshape(12)
        x1 = self.node_coords[0] + u[0:3]
        x2 = self.node_coords[1] + u[6:9]
        d = x2 - x1
        L = float(np.linalg.norm(d))
        if L <= 1e-18:
            raise ValueError("Corotational chord collapsed — check displacements.")
        R = _orthogonal_basis_from_chord(d)
        return _block_twelve(R), L

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        Tm, Lcurr = self._Tm_current(U_e)
        Kl = self._K_local(Lcurr)
        u = np.asarray(U_e, dtype=np.float64).reshape(12)
        return Tm.T @ Kl @ (Tm @ u)

    def _tangent_elastic_material(self, U_e: np.ndarray) -> np.ndarray:
        """Symmetric elastic material stiffness ``Tᵀ K_local T`` at current chord (no spin stiffness)."""
        Tm, Lcurr = self._Tm_current(U_e)
        Kl = self._K_local(Lcurr)
        K = Tm.T @ Kl @ Tm
        return 0.5 * (K + K.T)

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        """
        Tangent stiffness for Newton iterations.

        ``finite_difference`` (default): consistent Jacobian via central differences on
        ``internal_force_vector`` (``O(12)`` force evaluations).

        ``elastic_material``: analytic ``Tᵀ K_local(L) T`` only; omits corotational spin terms
        (see module docstring).
        """
        if self.tangent_stiffness_mode == "elastic_material":
            return self._tangent_elastic_material(U_e)

        u = np.asarray(U_e, dtype=np.float64).reshape(12)
        scale = np.sqrt(np.finfo(float).eps) * (1.0 + np.linalg.norm(u))
        if scale <= 0:
            scale = np.sqrt(np.finfo(float).eps)
        Kt = np.zeros((12, 12), dtype=np.float64)
        for j in range(12):
            du = np.zeros(12, dtype=np.float64)
            du[j] = scale
            fp = self.internal_force_vector(u + du)
            fm = self.internal_force_vector(u - du)
            Kt[:, j] = (fp - fm) / (2.0 * scale)
        return 0.5 * (Kt + Kt.T)

    def element_stiffness_matrix(self) -> ElementObject:
        self._assert_logging_ready()
        Tm0 = self._Tm_reference()
        Kl0 = self._K_local(self.L)
        Ke = Tm0.T @ Kl0 @ Tm0
        ref_el = self._linear_reference_element()
        eo = ref_el.element_stiffness_matrix()
        if self.logger_operator:
            self.logger_operator.log_text(
                "stiffness",
                f"\n=== Element {self.element_id} Corotational K_e (reference chord) ===",
            )
            self.logger_operator.log_matrix("stiffness", Ke, {"name": "Element Stiffness Matrix"})
            self.logger_operator.flush("stiffness")
        return ElementObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            K_e=Ke,
            gauss_data=[],
            integration_scheme="Gauss-Legendre",
            evaluate_shape_functions=eo.evaluate_shape_functions,
            shape_function_N_coefficients=eo.shape_function_N_coefficients,
            shape_function_dN_dxi_coefficients=eo.shape_function_dN_dxi_coefficients,
            shape_function_d2N_dxi2_coefficients=eo.shape_function_d2N_dxi2_coefficients,
        )

    def element_force_vector(self) -> ForceObject:
        fo = self._linear_reference_element().element_force_vector()
        return ForceObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            F_e=fo.F_e,
            gauss_data=fo.gauss_data,
            point_loads=fo.point_loads,
        )

    def element_mass_matrix(self) -> MassObject:
        mo = self._linear_reference_element().element_mass_matrix()
        Tm0 = self._Tm_reference()
        M_e = Tm0.T @ mo.M_e @ Tm0
        return MassObject(
            element_id=self.element_id,
            element_type=self.element_type_name,
            M_e=M_e,
        )
