# pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/curved_reddy/utilities/B_matrix.py
"""
Curved Reddy ``B`` (6, 12) per Gauss point in physical coordinates along the chord map
(``d/dx = (2/L) d/dxi``, ``detJ = L/2``) — STUB.

Voigt rows ``[eps_s, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]``.  ``kappa0`` couples
row 0 (axial-bending coupling ``-kappa0 * u_y``) and row 3 (shear-axial coupling
``+kappa0 * u_x``) on top of Reddy's alpha-corrected bending/shear strain.
``kappa0 -> 0`` recovers ``LinearReddyBeamElement3D``'s straight ``B``.

This file currently provides only the dataclass signature; ``physical_coordinate_form``
raises ``NotImplementedError`` until the curved-Reddy strain operator is filled in.

Parent element
--------------
``linear_curved_reddy_3D.LinearCurvedReddyBeamElement3D``:
``K_e += B.T @ D @ B * w_g * detJ``.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/curved_beam/utilities/B_matrix.py
pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/reddy/utilities/B_matrix.py
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CurvedStrainDisplacementOperator:
    """
    Builds curved-Reddy ``B`` (6, 12) for constant initial curvature ``kappa0`` — STUB.

    Parameters
    ----------
    element_length : float
        Chord length ``L`` (m).
    curvature : float
        Initial curvature ``kappa0 = 1/R`` (1/m); ``0`` for straight.
    alpha_coeff : float, optional
        Reddy alpha-coefficient (defaults to ``0`` until the linear element passes its
        own value); used by the curved-Reddy ``B`` rows once implemented.

    Notes
    -----
    Body of ``physical_coordinate_form`` raises ``NotImplementedError``; ``D`` is the
    straight Reddy section law (see ``linear_reddy_3D``).

    See Also
    --------
    pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.curved_beam.utilities.B_matrix.CurvedStrainDisplacementOperator
    pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.reddy.utilities.B_matrix
    """

    element_length: float
    curvature: float
    alpha_coeff: float = 0.0

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        object.__setattr__(self, "_jacobian", self.element_length / 2.0)
        object.__setattr__(self, "_dxi_dx", 2.0 / self.element_length)
        object.__setattr__(self, "_d2xi_dx2", 4.0 / self.element_length ** 2)

    @property
    def jacobian(self) -> float:
        return self._jacobian

    @property
    def dxi_dx(self) -> float:
        return self._dxi_dx

    @property
    def d2xi_dx2(self) -> float:
        return self._d2xi_dx2

    def physical_coordinate_form(
        self,
        dN_dxi: np.ndarray,
        d2N_dxi2: np.ndarray,
        N: np.ndarray | None = None,
    ) -> np.ndarray:
        """Return ``B`` (n_gauss, 6, 12) in physical coordinates — STUB body."""
        raise NotImplementedError(
            "CurvedStrainDisplacementOperator.physical_coordinate_form (curved Reddy): "
            "stub — implementation pending."
        )
