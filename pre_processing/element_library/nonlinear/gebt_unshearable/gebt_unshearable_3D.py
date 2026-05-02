# pre_processing/element_library/nonlinear/gebt_unshearable/gebt_unshearable_3D.py
"""
2-node 3D GEBT unshearable (Euler-Bernoulli) beam element (12 DOF) — STUB.

This module sits beside ``NonlinearTimoshenkoBeamElement3D`` (TL shear-deformable) and
``NonlinearEulerBernoulliBeamElement3D`` (Total-Lagrangian Euler-Bernoulli).  It is
the **unshearable** variant of GEBT: the geometrically-exact stack with the linear
Euler-Bernoulli kernel as ``K_0``; ``K_sigma`` is built from the EB nonlinear
stress resultants.

Tensors
-------
- ``U_e`` : (12,) — element nodal vector.
- ``B_lin`` : (4, 12) — linear Euler-Bernoulli ``StrainDisplacementOperator`` row
  layout (axial / bending-y / bending-z / torsion).
- ``D`` : (4, 4) — Euler-Bernoulli ``MaterialStiffnessOperator``
  (``EA``, ``EI_y``, ``EI_z``, ``GJ_t``).
- ``S = D @ E`` : (4,) — section forces (axial, M_y, M_z, T).
- ``K_0`` : (12, 12) — linear EB material stiffness with the standard EB selective
  integration (axial / torsion at axial order, bending at bending order).
- ``K_sigma`` : (12, 12) — EB geometric stiffness from
  ``GeometricStiffnessOperator(include_shear=False)``.
- ``K_T = K_0 + K_sigma``.

Weak forms (Gauss, ``xi`` in [-1, 1])
-------------------------------------
- ``F_int += B.T @ S * w_g * detJ`` with ``S = D @ E``
- ``K_T = K_0 + K_sigma`` (TL stack on linear EB)
- ``M_e`` : reference-configuration consistent mass (same shape functions as
  linear EB) for modal/dynamic analyses.
- ``detJ = L / 2``

Kinematics
----------
Geometrically-exact unshearable beam: shear strains ``γ_xy`` and ``γ_xz`` are
constrained out (Euler-Bernoulli kinematic constraint ``θ = -dw/dx``).  Tangent
stiffness sums material part (``K_0``, EB kernel) and geometric part (``K_sigma``).

Constitutive
------------
Linear-elastic Euler-Bernoulli ``MaterialStiffnessOperator``
(``EA``, ``EI_y``, ``EI_z``, ``GJ_t``).

Quadrature
----------
Gauss-Legendre.  ``quadrature_order`` defaults to ``max`` of the integration order
columns (axial / bending / torsion / load), minimum 2.

Public API
----------
- ``element_stiffness_matrix() -> ElementObject``
- ``element_force_vector() -> ForceObject``
- ``element_mass_matrix() -> MassObject``
- ``tangent_stiffness_matrix(U_e) -> np.ndarray`` (12, 12)
- ``internal_force_vector(U_e) -> np.ndarray`` (12,)

Limitations
-----------
This is a stub; the ``_get_K_0`` selective-integration routine (analogous to the
``NonlinearTimoshenkoBeamElement3D`` selective Timoshenko integration but with the EB row layout and
without shear blocks) is not yet implemented; assembly methods raise.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/nonlinear/timoshenko/nonlinear_timoshenko_3D.py
pre_processing/element_library/nonlinear/euler_bernoulli/nonlinear_euler_bernoulli_3D.py
pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/linear_euler_bernoulli_3D.py
"""

from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase

logger = logging.getLogger(__name__)


class GEBTUnshearableBeamElement3D(Element1DBase):
    """
    GEBT unshearable (Euler-Bernoulli kernel) 3D beam element (12 DOF) — STUB.

    Identity
    --------
    Two-node 3D beam with ``dof_per_node = 6``;
    ``element_type_name = "GEBTUnshearable-3D"``.

    Contract vs baseline
    --------------------
    Sits beside ``NonlinearTimoshenkoBeamElement3D`` and
    ``NonlinearEulerBernoulliBeamElement3D`` (TL EB).  Same TL pipeline idea as
    TL Timoshenko but using the linear EB ``B``, ``D``, the EB
    ``GreenLagrangeStrainOperator`` and ``GeometricStiffnessOperator`` with
    ``include_shear=False``.

    Notes
    -----
    Once implemented:

    - ``K_0`` cached from the linear EB kernel via ``_get_K_0`` (selective
      integration of axial, bending, torsion blocks; no shear blocks).
    - ``K_sigma`` is assembled from EB section forces (``N``, ``M_y``, ``M_z``).
    - At ``U_e = 0``, ``K_T == K_e`` of the linear EB element.

    See Also
    --------
    NonlinearTimoshenkoBeamElement3D
    NonlinearEulerBernoulliBeamElement3D
    LinearEulerBernoulliBeamElement3D
    """

    element_type_name = "GEBTUnshearable-3D"

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

        if quadrature_order is not None:
            self.quadrature_order = max(int(quadrature_order), 2)
        else:
            axial_order = int(self.element_array[3])
            bending_y_order = int(self.element_array[4])
            bending_z_order = int(self.element_array[5])
            torsion_order = int(self.element_array[8])
            load_order = int(self.element_array[9])
            self.quadrature_order = max(
                axial_order, bending_y_order, bending_z_order,
                torsion_order, load_order, 2
            )

        self.node_coords = self.grid_array
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = float(grid_dictionary["coordinates"][:, 0].min())
        self.x_global_end = float(grid_dictionary["coordinates"][:, 0].max())

        self._validate_element_properties()
        self._assert_logging_ready()

        self._K_0: np.ndarray | None = None

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
    def E(self) -> float:
        return float(self.material_array[0])

    @property
    def G(self) -> float:
        return float(self.material_array[1])

    @property
    def jacobian_determinant(self) -> float:
        return self.L / 2.0

    @property
    def integration_points(self) -> Tuple[np.ndarray, np.ndarray]:
        return np.polynomial.legendre.leggauss(self.quadrature_order)

    def _get_K_0(self) -> np.ndarray:
        """Material stiffness ``K_0`` with EB selective integration (no shear blocks).

        Marker method: in the GEBT unshearable pipeline ``K_0`` is the linear EB
        ``K_e``; body raises until the GEBT unshearable utilities tree lands.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__}._get_K_0: stub — implementation pending."
        )

    def element_stiffness_matrix(self):
        raise NotImplementedError(
            f"{self.__class__.__name__}: element_stiffness_matrix stub — implementation pending."
        )

    def element_force_vector(self):
        raise NotImplementedError(
            f"{self.__class__.__name__}: element_force_vector stub — implementation pending."
        )

    def element_mass_matrix(self):
        raise NotImplementedError(
            f"{self.__class__.__name__}: element_mass_matrix stub — implementation pending."
        )

    def tangent_stiffness_matrix(self, U_e: np.ndarray) -> np.ndarray:
        raise NotImplementedError(
            f"{self.__class__.__name__}: tangent_stiffness_matrix stub — implementation pending."
        )

    def internal_force_vector(self, U_e: np.ndarray) -> np.ndarray:
        raise NotImplementedError(
            f"{self.__class__.__name__}: internal_force_vector stub — implementation pending."
        )
