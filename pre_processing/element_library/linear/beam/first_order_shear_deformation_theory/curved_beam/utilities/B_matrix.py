# pre_processing/element_library/linear/curved_beam/utilities/B_matrix.py
"""
Curved Timoshenko ``B`` (6, 12) per Gauss point in physical coordinates along the chord map
(``d/dx = (2/L) d/dxi``, ``detJ = L/2``).

``eps`` rows: [eps_s, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x] — ``kappa0`` couples row 0
(``eps_s = d(u_x)/dx - kappa0 * u_y``) and row 3 (``gamma_xy = d(u_y)/dx + kappa0 * u_x - theta_z``).
Rows 1,2,4,5 match straight Timoshenko. ``kappa0 = 0`` recovers straight ``B``.

Parent element: ``K_e += B.T @ D @ B * w_g * detJ`` (``linear_curved_timoshenko_3D.py``).
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class CurvedStrainDisplacementOperator:
    """
    Builds ``B`` (6, 12) for constant initial curvature ``kappa0``; see module docstring for row definitions.

    ``physical_coordinate_form`` expects ``dN_dxi``, ``d2N_dxi2``, and ``N`` (for ``kappa0 != 0`` coupling terms).

    Notes
    -----
    Weak-form linkage: ``linear_curved_timoshenko_3D.LinearCurvedTimoshenkoBeamElement3D``; ``D`` is straight Timoshenko section law.

    See Also
    --------
    linear_curved_timoshenko_3D.LinearCurvedTimoshenkoBeamElement3D
    """

    element_length: float
    curvature: float  # κ0 = 1/R (1/m); 0 for straight

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        object.__setattr__(self, "_jacobian", self.element_length / 2)
        object.__setattr__(self, "_dξ_dx", 2 / self.element_length)
        object.__setattr__(self, "_d2ξ_dx2", 4 / self.element_length ** 2)

    @property
    def jacobian(self) -> float:
        return self._jacobian

    @property
    def dξ_dx(self) -> float:
        return self._dξ_dx

    @property
    def d2ξ_dx2(self) -> float:
        return self._d2ξ_dx2

    def physical_coordinate_form(
        self,
        dN_dξ: np.ndarray,
        d2N_dξ2: np.ndarray,
        N: np.ndarray = None,
    ) -> np.ndarray:
        """B in physical coords (6, 12) per GP. Curvature κ0 couples axial with u_y and shear with u_x."""
        n_gauss = dN_dξ.shape[0]
        B = np.zeros((n_gauss, 6, 12))
        dξ_dx = self.dξ_dx
        κ0 = self.curvature

        # Axial: ε_s = du_x/dx - κ0*u_y
        B[:, 0, [0, 6]] = dN_dξ[:, [0, 6], 0] * dξ_dx
        if κ0 != 0 and N is not None:
            B[:, 0, [1, 7]] = -κ0 * N[:, [1, 7], 1]

        # Bending (unchanged)
        B[:, 1, [4, 10]] = dN_dξ[:, [4, 10], 4] * dξ_dx
        B[:, 2, [5, 11]] = dN_dξ[:, [5, 11], 5] * dξ_dx

        # Shear γ_xy = du_y/dx + κ0*u_x - θ_z
        if N is not None:
            B[:, 3, [1, 7]] = dN_dξ[:, [1, 7], 1] * dξ_dx
            B[:, 3, [5, 11]] = -N[:, [5, 11], 5]
            if κ0 != 0:
                B[:, 3, [0, 6]] = κ0 * N[:, [0, 6], 0]

        # Shear γ_xz (unchanged)
        if N is not None:
            B[:, 4, [2, 8]] = dN_dξ[:, [2, 8], 2] * dξ_dx
            B[:, 4, [4, 10]] = -N[:, [4, 10], 4]

        # Torsion
        B[:, 5, [3, 9]] = dN_dξ[:, [3, 9], 3] * dξ_dx

        return B
