# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/curved_beam/utilities/B_matrix.py
"""
Curved Euler-Bernoulli ``B`` (4, 12) per Gauss point in physical coordinates along the chord
map (``d/dx = (2/L) d/dxi``, ``detJ = L/2``) — STUB.

Voigt order
-----------
``eps`` rows: [eps_s, kappa_y, kappa_z, phi_x] (no transverse shear; ZOSDT).
- ``kappa0`` couples row 0 (``eps_s = d(u_x)/dx - kappa0 * u_y``) and modifies row 2
  (``kappa_z = d^2(u_y)/dx^2 + kappa0 * d(u_x)/dx``) for an in-plane circular arc.
- Rows 1, 3 unchanged from straight EB.
- ``kappa0 -> 0`` recovers straight EB ``B``.

This file currently provides only the dataclass signature; ``physical_coordinate_form``
raises ``NotImplementedError`` until the curved-EB strain operator is filled in.

Parent element
--------------
``linear/beam/zero_order_shear_deformation_theory/curved_beam/linear_curved_euler_bernoulli_3D.py``:
``K_e += B.T @ D @ B * w_g * detJ``.

See Also
--------
docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/curved_beam/utilities/B_matrix.py
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CurvedStrainDisplacementOperator:
    """
    Builds ``B`` (4, 12) for constant initial curvature ``kappa0`` in the chord plane.

    This is a STUB.  The constructor mirrors the FOSDT curved-beam analogue but
    ``physical_coordinate_form`` raises ``NotImplementedError`` until the curved-EB row
    expressions are populated.

    Parameters
    ----------
    element_length : float
        Chord length ``L`` (m).
    curvature : float
        Initial curvature ``kappa0 = 1/R`` (1/m); ``0`` for straight.

    Notes
    -----
    Weak-form linkage: ``LinearCurvedEulerBernoulliBeamElement3D``; ``D`` is straight EB
    section law (``EA``, ``EI_y``, ``EI_z``, ``GJ_t``).  Body raises until implemented.

    See Also
    --------
    pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.curved_beam.utilities.B_matrix.CurvedStrainDisplacementOperator
    """

    element_length: float
    curvature: float

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
        """Return ``B`` (n_gauss, 4, 12) in physical coordinates — STUB body."""
        raise NotImplementedError(
            "CurvedStrainDisplacementOperator.physical_coordinate_form (curved EB): "
            "stub — implementation pending."
        )
