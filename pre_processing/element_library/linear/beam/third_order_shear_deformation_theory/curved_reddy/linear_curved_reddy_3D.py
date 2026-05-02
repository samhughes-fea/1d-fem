# pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/curved_reddy/linear_curved_reddy_3D.py
"""
2-node 3D linear curved Reddy 3rd-order beam (12 DOF), constant ``kappa0`` — STUB.

This module ports the curved-beam structure used by ``LinearCurvedTimoshenkoBeamElement3D``
to Reddy 3rd-order kinematics.  The class scaffolding and operator wiring are
documented exactly as for the straight ``LinearReddyBeamElement3D``, but every
assembly method raises ``NotImplementedError`` until the curved-Reddy ``B`` body is
filled in.

Tensors
-------
- ``U_e`` : (12,) — standard 3D beam DOFs (``dof_per_node = 6``).
- ``B`` : (6, 12) — strain-displacement per Gauss point from the new
  ``CurvedStrainDisplacementOperator`` (this module's ``utilities/B_matrix.py``); same
  Voigt order as Reddy with curvature coupling in rows 0 and 3.
- ``D`` : (6, 6) — straight Reddy ``MaterialStiffnessOperator`` (with the
  alpha-coefficient cross-terms).
- ``K_e`` : (12, 12), ``F_e`` : (12,).

Weak forms (Gauss, ``xi`` in [-1, 1])
-------------------------------------
- ``K_e += B.T @ D @ B * w_g * detJ``
- ``F_dist += w_g * N.T @ q * detJ``
- ``F_point = N.T @ P``
- ``detJ = L / 2``

Kinematics
----------
Reddy third-order shear deformation theory on a chord-frame curve of constant initial
curvature ``kappa0`` (``element_dictionary["curvature"]``).  ``kappa0 -> 0`` recovers
``LinearReddyBeamElement3D``.

Constitutive
------------
Straight Reddy ``MaterialStiffnessOperator``; reused unchanged.

Quadrature
----------
Gauss-Legendre.  ``quadrature_order`` defaults to ``max`` of integration order columns
(shear at least ``3``).

Public API
----------
- ``element_stiffness_matrix() -> ElementObject``
- ``element_force_vector() -> ForceObject``
- ``element_mass_matrix() -> MassObject``

Limitations
-----------
The ``CurvedStrainDisplacementOperator`` body raises until the curved-Reddy ``B`` is
filled in.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/curved_beam/linear_curved_timoshenko_3D.py
pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/reddy/linear_reddy_3D.py
"""

from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.curved_reddy.utilities.B_matrix import (
    CurvedStrainDisplacementOperator,
)

logger = logging.getLogger(__name__)


class LinearCurvedReddyBeamElement3D(Element1DBase):
    """
    Linear curved Reddy 3rd-order beam element (12 DOF) — STUB.

    Identity
    --------
    Two-node 3D beam with ``dof_per_node = 6``;
    ``element_type_name = "CurvedReddy-3D"``.

    Contract vs baseline
    --------------------
    Composes ``LinearReddyBeamElement3D``'s ``MaterialStiffnessOperator`` and shape
    function registry with a new ``CurvedStrainDisplacementOperator`` that adds
    curvature coupling on top of Reddy's alpha-corrected bending/shear.  ``kappa0 -> 0``
    recovers ``LinearReddyBeamElement3D``.

    Notes
    -----
    Once implemented, ``B`` will be 6 x 12 with rows 0 and 3 coupled by ``kappa0``.

    See Also
    --------
    LinearReddyBeamElement3D
    LinearCurvedTimoshenkoBeamElement3D
    """

    element_type_name = "CurvedReddy-3D"

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

        idx = int(np.where(element_dictionary["ids"] == element_id)[0][0])
        self.curvature = float(
            element_dictionary.get(
                "curvature",
                np.zeros(len(element_dictionary["ids"]), dtype=np.float64),
            )[idx]
        )

        if quadrature_order is not None:
            self.quadrature_order = int(quadrature_order)
        else:
            axial_order = int(self.element_array[3])
            bending_y_order = int(self.element_array[4])
            bending_z_order = int(self.element_array[5])
            shear_y_order = int(self.element_array[6]) or 3
            shear_z_order = int(self.element_array[7]) or 3
            torsion_order = int(self.element_array[8])
            load_order = int(self.element_array[9])
            self.quadrature_order = max(
                axial_order, bending_y_order, bending_z_order,
                shear_y_order, shear_z_order, torsion_order, load_order, 3
            )

        self.node_coords = self.grid_array
        self.L = float(np.linalg.norm(self.node_coords[1] - self.node_coords[0]))
        self.x_start, *_, self.x_end = self.node_coords[[0, 1], 0]
        self.x_global_start = float(grid_dictionary["coordinates"][:, 0].min())
        self.x_global_end = float(grid_dictionary["coordinates"][:, 0].max())

        self._validate_element_properties()
        self._assert_logging_ready()

        self.strain_displacement_operator = CurvedStrainDisplacementOperator(
            element_length=self.L,
            curvature=self.curvature,
            alpha_coeff=self.alpha_coeff,
        )

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
    def alpha_coeff(self) -> float:
        if self.section_array.size >= 7:
            return float(self.section_array[6])
        return float(self.I_z) / float(self.A) if self.A > 0 else 0.0

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
