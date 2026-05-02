# pre_processing/element_library/nonlinear/euler_bernoulli/updated_lagrangian_euler_bernoulli_3D.py
"""
Updated-Lagrangian (UL) 3D Euler-Bernoulli beam element (12 DOF) — STUB.

This module sits beside ``NonlinearEulerBernoulliBeamElement3D`` (Total Lagrangian).
The reference configuration of every Gauss point integral is the **last converged
configuration** — *not* the initial geometry — so ``B``, ``D``, and the geometric
stiffness ``K_sigma`` are recomputed about the current state on every step.  See
``docs/element_library/total_lagrangian_beam_formulation.md`` for the contrasting
TL setting.

**Strain/stress (planned implementation):** Not the TL Green–Lagrange / reference-\\(X\\) pipeline in ``nonlinear/.../utilities/green_lagrange_strain.py``; UL uses small-strain theory on the **updated** configuration with stresses conjugate to that choice (see ``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md``).

Tensors
-------
- ``U_e`` : (12,) — element nodal vector in the current configuration's local frame.
- ``B_lin`` : (4, 12) — linear EB strain-displacement evaluated about the **current**
  geometry (Almansi or Cauchy, both reduce to the small-strain ``B`` once the chord
  derivatives are taken on the deformed length).
- ``D`` : (4, 4) — straight EB ``MaterialStiffnessOperator`` (small-strain elastic
  modulus).
- ``S`` : (4,) — current section forces ``S = D @ E``.
- ``K_mat`` : (12, 12) — material stiffness about the current state.
- ``K_sigma`` : (12, 12) — geometric stiffness on the current state.
- ``K_T = K_mat + K_sigma`` (UL).

Weak forms (Gauss, ``xi`` in [-1, 1])
-------------------------------------
- ``K_e += B.T @ D @ B * w_g * detJ`` evaluated on the **current** chord
- ``F_int += B.T @ S * w_g * detJ`` with ``S = D @ E``
- ``K_T = K_mat + K_sigma`` (UL; reference is the last converged configuration)
- ``detJ = L_curr / 2`` (chord map on the current chord)

Kinematics
----------
Updated-Lagrangian.  After each load increment the converged configuration becomes
the new reference.  Strains are the small (Almansi or Cauchy) strain on the updated
reference; the constitutive law is unchanged.  This contrasts with the TL counterpart
``NonlinearEulerBernoulliBeamElement3D`` whose reference is fixed at the initial
geometry.

Constitutive
------------
Straight EB ``MaterialStiffnessOperator`` (``EA``, ``EI_y``, ``EI_z``, ``GJ_t``).

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
The Updated-Lagrangian utilities (Almansi or Cauchy strain operator,
current-configuration ``K_sigma``, configuration update API) are not yet provided;
``_get_K_0_current()`` is a docstring-only marker that the material/geometric tangents
must be evaluated about the current configuration (see ``Notes``).  All assembly
methods raise ``NotImplementedError`` until that lands.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
docs/element_library/total_lagrangian_beam_formulation.md
pre_processing/element_library/nonlinear/euler_bernoulli/nonlinear_euler_bernoulli_3D.py
"""

from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase

logger = logging.getLogger(__name__)


class UpdatedLagrangianEulerBernoulliBeamElement3D(Element1DBase):
    """
    Updated-Lagrangian 3D Euler-Bernoulli beam element (12 DOF) — STUB.

    Identity
    --------
    Two-node 3D beam with ``dof_per_node = 6``;
    ``element_type_name = "EulerBernoulli-3D-UpdatedLagrangian"``.

    Contract vs baseline
    --------------------
    Sits beside ``NonlinearEulerBernoulliBeamElement3D`` (TL counterpart).  Both share
    the linear EB ``B`` and ``D``; the difference is the **reference**:
    UL evaluates everything about the last converged configuration, TL keeps the
    reference fixed at the initial geometry.

    Notes
    -----
    Once implemented:

    - ``_get_K_0_current()`` returns ``K_mat`` evaluated on the current chord.
    - ``K_sigma`` is built on the current chord with ``dN/dx_curr``.
    - ``S`` is the small-strain Cauchy stress times ``D`` on the updated reference.

    See Also
    --------
    NonlinearEulerBernoulliBeamElement3D
    LinearEulerBernoulliBeamElement3D
    """

    element_type_name = "EulerBernoulli-3D-UpdatedLagrangian"

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

    def _get_K_0_current(self) -> np.ndarray:
        """Material stiffness about the **current** configuration (UL).

        Marker method: in the UL pipeline ``K_mat`` is rebuilt on the converged chord
        every step, so there is no fixed initial cache.  Body raises until the UL
        utilities tree lands.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__}._get_K_0_current: stub — implementation pending."
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
