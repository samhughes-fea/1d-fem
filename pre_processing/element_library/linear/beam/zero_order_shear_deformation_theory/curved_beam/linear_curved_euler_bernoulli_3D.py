# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/curved_beam/linear_curved_euler_bernoulli_3D.py
"""
2-node 3D linear curved Euler-Bernoulli beam, constant ``kappa0`` (initial curvature in
the chord plane) — STUB.

This module ports the curved-beam structure used by ``LinearCurvedTimoshenkoBeamElement3D``
to the zero-order shear deformation theory (no transverse shear DOFs).  The class
scaffolding and operator wiring are documented exactly as for that analogue, but every
assembly method raises ``NotImplementedError`` until the curved-EB ``B`` and ``D``
expressions are filled in.

Tensors
-------
- ``U_e`` : (12,) — standard 3D beam DOFs (``dof_per_node = 6``).
- ``B`` : (4, 12) — strain-displacement per Gauss point, Voigt rows
  ``[eps_s, kappa_y, kappa_z, phi_x]`` from the new ZOSDT
  ``CurvedStrainDisplacementOperator`` (this module's ``utilities/B_matrix.py``).
- ``D`` : (4, 4) — straight EB ``MaterialStiffnessOperator`` reduced to
  ``diag(EA, EI_y, EI_z, GJ_t)``.
- ``K_e`` : (12, 12), ``F_e`` : (12,).

Weak forms (Gauss, ``xi`` in [-1, 1])
-------------------------------------
- ``K_e += B.T @ D @ B * w_g * detJ``
- ``F_dist += w_g * N.T @ q * detJ`` (load mapping per
  ``FORMULATION_DOCSTRING_STANDARDS.md``)
- ``F_point = N.T @ P``
- ``detJ = L / 2`` (chord map ``d/dx = (2/L) d/dxi``)

Kinematics
----------
ZOSDT (Euler-Bernoulli; rotations slaved to ``d(u_y)/dx`` and ``d(u_z)/dx``).  Constant
initial curvature ``kappa0`` in the chord plane couples row 0 of ``B`` (axial) with
``u_y`` and modifies row 2 (in-plane bending) with ``d(u_x)/dx``.  ``kappa0 -> 0``
recovers ``LinearEulerBernoulliBeamElement3D``.

Constitutive
------------
Straight EB ``MaterialStiffnessOperator`` (``EA``, ``EI_y``, ``EI_z``, ``GJ_t``); reuses
the existing ZOSDT ``D`` operator.

Quadrature
----------
Gauss-Legendre.  ``quadrature_order`` defaults to
``max(axial, bending_y, bending_z, torsion, load, 2)`` from the element integration
order block, mirroring ``LinearEulerBernoulliBeamElement3D``.

Public API
----------
- ``element_stiffness_matrix() -> ElementObject``
- ``element_force_vector() -> ForceObject``
- ``element_mass_matrix() -> MassObject``

Limitations
-----------
The body of ``CurvedStrainDisplacementOperator.physical_coordinate_form`` raises
``NotImplementedError`` (see ``utilities/B_matrix.py``); all assembly methods on this
class therefore raise as well.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/curved_beam/linear_curved_timoshenko_3D.py
pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/linear_euler_bernoulli_3D.py
"""

from __future__ import annotations

import logging
import warnings
from typing import Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.curved_beam.utilities.B_matrix import (
    CurvedStrainDisplacementOperator,
)

logger = logging.getLogger(__name__)


class LinearCurvedEulerBernoulliBeamElement3D(Element1DBase):
    """
    Linear curved Euler-Bernoulli element (12 DOF) — STUB.

    Identity
    --------
    Two-node 3D beam with ``dof_per_node = 6``;
    ``element_type_name = "CurvedEulerBernoulli-3D"``.

    Contract vs baseline
    --------------------
    Mirrors ``LinearCurvedTimoshenkoBeamElement3D`` but on ZOSDT (no transverse shear);
    ``kappa0 -> 0`` recovers ``LinearEulerBernoulliBeamElement3D``.

    Notes
    -----
    Once implemented:

    - ``B`` is provided by ``CurvedStrainDisplacementOperator`` (this module's
      ``utilities/B_matrix.py``); body currently raises ``NotImplementedError``.
    - ``D`` is the straight EB ``MaterialStiffnessOperator``.
    - ``N`` is supplied by the standard EB shape function registry.

    See Also
    --------
    LinearCurvedTimoshenkoBeamElement3D
    LinearEulerBernoulliBeamElement3D
    """

    element_type_name = "CurvedEulerBernoulli-3D"

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
        warnings.warn(
            "LinearCurvedEulerBernoulliBeamElement3D is deprecated; use "
            "LinearEulerBernoulliBeamElement3D with precurvature.txt (reference k_x0, k_y0, k_z0).",
            DeprecationWarning,
            stacklevel=2,
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

        self._validate_element_properties()
        self._assert_logging_ready()

        self.strain_displacement_operator = CurvedStrainDisplacementOperator(
            element_length=self.L, curvature=self.curvature
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
