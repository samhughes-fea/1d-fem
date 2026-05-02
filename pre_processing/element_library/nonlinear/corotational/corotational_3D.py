# pre_processing/element_library/nonlinear/corotational/corotational_3D.py
"""
Co-rotational 3D beam element (12 DOF) — STUB.

This module sits beside ``NonlinearTimoshenkoBeamElement3D`` and the Total/Updated-Lagrangian
elements.  The framework is **co-rotational**: a moving local frame ``T(U_e)`` is
extracted from the current displacements every step, the small-strain linear EB or
Timoshenko kernel is applied **in that frame**, and the result is rotated back to
the global element frame with a frame-rotation derivative term.

Tensors
-------
- ``U_e`` : (12,) — element nodal vector in the global element frame.
- ``T(U_e)`` : (12, 12) — block-diagonal frame transformation (per-node ``3 x 3``
  rotations on translations and rotations).
- ``K_local`` : (12, 12) — small-strain linear EB or Timoshenko stiffness in the
  current local frame.
- ``K_geom_frame`` : (12, 12) — geometric stiffness contribution from the frame
  rotation derivatives.
- ``K_T = T.T @ K_local @ T + K_geom_frame``.

Weak forms (Gauss, ``xi`` in [-1, 1])
-------------------------------------
- ``K_local += B.T @ D @ B * w_g * detJ`` (small-strain kernel)
- ``F_int = T.T @ F_local`` with ``F_local = K_local @ U_local``
- ``K_T = T.T @ K_local @ T + K_geom_frame``
- ``detJ = L_curr / 2``

Kinematics
----------
Co-rotational large-displacement / small-strain.  At each step the frame ``T`` is
obtained from a chord and node-rotation extraction; the strain inside ``T`` is the
small-strain measure (Euler-Bernoulli or Timoshenko).  Distinct from TL (fixed
reference) and UL (last converged reference): co-rotational separates the rigid
rotation cleanly.

Constitutive
------------
Whichever small-strain kernel is selected by ``kernel`` (default ``"timoshenko"``);
the linear EB or Timoshenko ``MaterialStiffnessOperator`` provides ``D``.

Quadrature
----------
Gauss-Legendre.  ``quadrature_order`` defaults to ``max`` of the integration order
columns.

Public API
----------
- ``element_stiffness_matrix() -> ElementObject``
- ``element_force_vector() -> ForceObject``
- ``element_mass_matrix() -> MassObject``
- ``tangent_stiffness_matrix(U_e) -> np.ndarray`` (12, 12)
- ``internal_force_vector(U_e) -> np.ndarray`` (12,)

Limitations
-----------
The co-rotational frame algebra (``T(U_e)``, frame-rotation derivatives,
``K_geom_frame``) is not yet implemented; assembly methods raise.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/nonlinear/timoshenko/nonlinear_timoshenko_3D.py
pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/linear_timoshenko_3D.py
pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/linear_euler_bernoulli_3D.py
"""

from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase

logger = logging.getLogger(__name__)


class CorotationalBeamElement3D(Element1DBase):
    """
    Co-rotational 3D beam element (12 DOF) — STUB.

    Identity
    --------
    Two-node 3D beam with ``dof_per_node = 6``;
    ``element_type_name = "Corotational-3D"``.

    Contract vs baseline
    --------------------
    Sits beside TL beam elements.  Uses the linear EB or Timoshenko kernel
    in a moving local frame; ``kernel`` constructor argument switches between them
    (default ``"timoshenko"``).

    Notes
    -----
    Once implemented:

    - ``T(U_e)`` is the per-node frame transform; ``K_local`` is the linear kernel's
      ``K_e`` evaluated on the current chord.
    - ``K_geom_frame`` captures the derivative of ``T`` w.r.t. ``U_e``.
    - ``F_int = T.T @ K_local @ U_local``.

    See Also
    --------
    NonlinearTimoshenkoBeamElement3D
    LinearTimoshenkoBeamElement3D
    LinearEulerBernoulliBeamElement3D
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

        if kernel not in ("euler_bernoulli", "timoshenko"):
            raise ValueError(
                f"kernel must be 'euler_bernoulli' or 'timoshenko', got {kernel!r}"
            )
        self.kernel = kernel

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
                axial_order, bending_y_order, bending_z_order,
                shear_y_order, shear_z_order, torsion_order, load_order, 2
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
