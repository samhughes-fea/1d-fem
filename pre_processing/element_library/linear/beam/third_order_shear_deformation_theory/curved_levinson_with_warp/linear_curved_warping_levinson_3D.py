# pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/curved_levinson_with_warp/linear_curved_warping_levinson_3D.py
"""
2-node 3D linear curved + warping Levinson 3rd-order beam (14 DOF) — STUB.

This module combines two extensions on top of ``LinearLevinsonBeamElement3D``:

1. Initial curvature ``kappa0`` from ``element_dictionary["curvature"]`` (rows 0 and 3
   of ``B`` coupled, see ``curved_levinson/utilities/B_matrix.py``).
2. Vlasov warping intensity ``chi`` carried as a seventh DOF per node
   (``dof_per_node = 7``).

Tensors
-------
- ``U_e`` : (14,) — element nodal vector with ``chi`` as the seventh DOF per node.
- ``B`` : (7, 14) — strain-displacement; rows 0..5 reuse the curved Levinson ``B``
  (with ``kappa0`` coupling) on columns 0..11; row 6 is the Vlasov warping rate
  ``d(chi)/dx``.
- ``D`` : (7, 7) — Levinson ``MaterialStiffnessOperator`` extended with
  ``D[6, 6] = E * Gamma``.
- ``K_e`` : (14, 14), ``F_e`` : (14,).

Weak forms (Gauss, ``xi`` in [-1, 1])
-------------------------------------
- ``K_e += B.T @ D @ B * w_g * detJ``
- ``F_dist += w_g * N.T @ q * detJ``
- ``F_point = N.T @ P``
- ``detJ = L / 2``

Kinematics
----------
Levinson 3rd-order kinematics on a chord-frame curve of constant initial curvature
``kappa0``, plus a Vlasov warping intensity ``chi``.  Limit cases:

- ``kappa0 = 0``         → ``LinearWarpingLevinsonBeamElement3D``
- ``Gamma = 0``          → ``LinearCurvedLevinsonBeamElement3D``
- both zero              → ``LinearLevinsonBeamElement3D``

Constitutive
------------
Levinson ``MaterialStiffnessOperator`` extended with ``D[6, 6] = E * Gamma``.

Quadrature
----------
Gauss-Legendre.  ``quadrature_order`` defaults to ``max`` of the integration order
columns (shear at least ``3``).

Public API
----------
- ``element_stiffness_matrix() -> ElementObject``
- ``element_force_vector() -> ForceObject``
- ``element_mass_matrix() -> MassObject``

Limitations
-----------
The curved-Levinson ``CurvedStrainDisplacementOperator`` body raises; the warping
extension to ``B`` is not yet wired in this stub.  All assembly methods raise.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/curved_levinson/linear_curved_levinson_3D.py
pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/levinson_with_warp/linear_warping_levinson_3D.py
"""

from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

from pre_processing.element_library.element_1D_base import Element1DBase

logger = logging.getLogger(__name__)


class LinearCurvedWarpingLevinsonBeamElement3D(Element1DBase):
    """
    Linear curved Levinson 3rd-order beam with Vlasov warping (14 DOF) — STUB.

    Identity
    --------
    Two-node 3D beam with ``dof_per_node = 7``;
    ``element_type_name = "CurvedWarpingLevinson-3D"``.

    Contract vs baseline
    --------------------
    Combines ``LinearCurvedLevinsonBeamElement3D`` and
    ``LinearWarpingLevinsonBeamElement3D``.  Both reductions
    (``kappa0 = 0`` or ``Gamma = 0``) recover the corresponding single-extension
    element; both zero recovers ``LinearLevinsonBeamElement3D``.

    Notes
    -----
    Once implemented, ``B`` (7, 14) is the curved-Levinson 6-row matrix promoted to
    width 14 with the Vlasov warping rate as row 6; ``D`` (7, 7) is the Levinson
    section law plus ``D[6, 6] = E * Gamma``.

    See Also
    --------
    LinearLevinsonBeamElement3D
    LinearCurvedLevinsonBeamElement3D
    LinearWarpingLevinsonBeamElement3D
    """

    element_type_name = "CurvedWarpingLevinson-3D"

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
            dof_per_node=7,
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
    def Gamma(self) -> float:
        if self.section_array.size >= 10:
            return float(self.section_array[9])
        return 0.0

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
