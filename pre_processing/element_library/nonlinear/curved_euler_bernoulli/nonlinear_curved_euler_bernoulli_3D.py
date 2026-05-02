# pre_processing/element_library/nonlinear/curved_euler_bernoulli/nonlinear_curved_euler_bernoulli_3D.py
"""
Nonlinear (Total Lagrangian) 3D curved Euler-Bernoulli beam element (12 DOF) — STUB.

This module is a STUB.  It bolts the Total-Lagrangian pipeline used by
``NonlinearEulerBernoulliBeamElement3D`` onto the curved-EB ``B`` operator introduced
under ``linear/beam/zero_order_shear_deformation_theory/curved_beam`` and raises
``NotImplementedError`` from every assembly method until the curved-EB ``B`` body is
filled in.

Tensors
-------
- ``U_e`` : (12,) — standard 3D beam DOFs (``dof_per_node = 6``).
- ``B_lin`` : (4, 12) — linear curved-EB strain-displacement from the new
  ``CurvedStrainDisplacementOperator`` (ZOSDT); rows ``[eps_s, kappa_y, kappa_z, phi_x]``.
- ``B_nl`` : (4, 12) — nonlinear (Green-Lagrange) increment; reuses TL EB
  ``GreenLagrangeStrainOperator.nonlinear_strain_displacement_gradient`` on the chord
  derivatives.
- ``D`` : (4, 4) — straight EB ``MaterialStiffnessOperator`` (curvature does not change
  the constitutive law).
- ``S`` : (4,) — section forces ``S = D @ E`` (axial N, two bending moments, torsion T).
- ``E`` : (4,) — Green-Lagrange strain ``E = E_lin + E_nl``.
- ``K_0`` : (12, 12), ``K_sigma`` : (12, 12), ``K_T = K_0 + K_sigma``.

Weak forms (Gauss, ``xi`` in [-1, 1])
-------------------------------------
- ``K_e += B.T @ D @ B * w_g * detJ`` (linear cache: ``B = B_lin``)
- ``F_int += B_lin.T @ S * w_g * detJ`` with ``S = D @ E``
- ``K_T = K_0 + K_sigma`` or ``K_T = K_mat + K_sigma`` (full nonlinear curvature)
- ``detJ = L / 2`` (chord map)

Kinematics
----------
ZOSDT (Euler-Bernoulli) on a chord-frame curve of constant initial curvature ``kappa0``.
``kappa0`` is supplied through ``element_dictionary["curvature"]`` and feeds
``CurvedStrainDisplacementOperator`` (this stub's parent linear class).  Strains are
referred to the **initial** (undeformed) configuration.  ``kappa0 -> 0`` recovers
``NonlinearEulerBernoulliBeamElement3D``.

Constitutive
------------
Straight EB ``MaterialStiffnessOperator`` (``EA``, ``EI_y``, ``EI_z``, ``GJ_t``).

Quadrature
----------
Gauss-Legendre.  ``quadrature_order`` defaults to
``max(axial, bending_y, bending_z, torsion, load, 2)`` from the element integration
order block.

Public API
----------
- ``element_stiffness_matrix() -> ElementObject``
- ``element_force_vector() -> ForceObject``
- ``element_mass_matrix() -> MassObject``
- ``tangent_stiffness_matrix(U_e) -> np.ndarray`` (12, 12)
- ``internal_force_vector(U_e) -> np.ndarray`` (12,)

Limitations
-----------
Depends on the linear curved-EB ``CurvedStrainDisplacementOperator`` (whose body raises);
all assembly methods raise as a result.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/curved_beam/linear_curved_euler_bernoulli_3D.py
pre_processing/element_library/nonlinear/euler_bernoulli/nonlinear_euler_bernoulli_3D.py
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


class NonlinearCurvedEulerBernoulliBeamElement3D(Element1DBase):
    """
    Total-Lagrangian 3D curved Euler-Bernoulli beam element (12 DOF) — STUB.

    Identity
    --------
    Two-node 3D beam with ``dof_per_node = 6``;
    ``element_type_name = "CurvedEulerBernoulli-3D-Nonlinear"``.

    Contract vs baseline
    --------------------
    Embeds ``LinearCurvedEulerBernoulliBeamElement3D`` (linear curved ``B``, ``D``,
    ``element_mass_matrix``) and overlays the Total-Lagrangian pipeline used by
    ``NonlinearEulerBernoulliBeamElement3D``.  ``kappa0 -> 0`` recovers
    ``NonlinearEulerBernoulliBeamElement3D``; ``U_e = 0`` recovers
    ``LinearCurvedEulerBernoulliBeamElement3D``.

    Notes
    -----
    Once implemented:

    - ``K_0`` will equal the cached ``element_stiffness_matrix()`` result of the linear
      curved-EB element.
    - ``K_sigma`` will be the geometric stiffness from
      ``GeometricStiffnessOperator.assemble_K_sigma`` evaluated with ``dN/dx`` along the
      chord and the curvature-aware ``B`` for the bending block.
    - ``B_lin`` is from the curved-EB strain operator (4 x 12).
    - ``B_nl`` is the TL EB nonlinear increment (axial GL term).
    - ``S`` is ``D @ (E_lin + E_nl)``.

    See Also
    --------
    LinearCurvedEulerBernoulliBeamElement3D
    NonlinearEulerBernoulliBeamElement3D
    """

    element_type_name = "CurvedEulerBernoulli-3D-Nonlinear"

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
            "NonlinearCurvedEulerBernoulliBeamElement3D is deprecated; use "
            "NonlinearEulerBernoulliBeamElement3D with precurvature.txt (reference k_x0, k_y0, k_z0).",
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
