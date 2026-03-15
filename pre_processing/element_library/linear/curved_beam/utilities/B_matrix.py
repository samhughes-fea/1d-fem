# pre_processing/element_library/linear/curved_beam/utilities/B_matrix.py
"""Strain-displacement matrix B for 2-node 3-D curved Timoshenko beam. Constant curvature κ0 in x-y plane; ε = [ε_s, κ_y, κ_z, γ_xy, γ_xz, φ_x] with curvature coupling in axial and shear. B shape (6, 12) per GP. When κ0=0 reduces to straight Timoshenko."""

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class CurvedStrainDisplacementOperator:
    """
    B-matrix for 2-node curved Timoshenko beam with constant initial curvature κ0 (in x-y plane).
    Strain: ε_s = du_x/ds - κ0*u_y; κ_y, κ_z unchanged; γ_xy = du_y/ds + κ0*u_x - θ_z; γ_xz, φ_x unchanged.
    When κ0=0 identical to straight Timoshenko.
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
