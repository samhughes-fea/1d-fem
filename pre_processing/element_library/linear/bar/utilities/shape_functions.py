# pre_processing/element_library/linear/bar/utilities/shape_functions.py
"""
Bar shape functions: ``natural_coordinate_form(xi_vec)`` returns ``N``, ``dN_dxi``, ``d2N_dxi2`` with shape
(n_gp, 12, 6). Distributed loads: ``F_dist += w_g * N.T @ q * detJ``, ``detJ = L/2``, ``xi in [-1, 1]``.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class ShapeFunctionOperator:
    """
    Operator for evaluating 3D Bar element shape functions and their derivatives.
    Bar has axial and torsion only; shape functions are linear Lagrange in ξ.

    Mathematical Formulation
    -----------------------
    - Axial displacement: Linear Lagrange N₁(ξ) = 0.5(1−ξ), N₂(ξ) = 0.5(1+ξ) for DOF 0, 6.
    - Torsional rotation: Same linear Lagrange for DOF 3, 9.

    Coordinate Transformation:
    - Physical to natural: ξ = (2x - L)/L
    - ∂N/∂x = (∂N/∂ξ)(∂ξ/∂x) = (∂N/∂ξ)(2/L)

    Parameters
    ----------
    element_length : float
        Physical length of element (x ∈ [0,L], L > 0)

    Attributes
    ----------
    dξ_dx : float
        First derivative transform (∂ξ/∂x = 2/L)
    d2ξ_dx2 : float
        Second derivative transform (∂²ξ/∂x² = 4/L²); zero used for bar (no bending).
    """

    element_length: float

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        object.__setattr__(self, '_dξ_dx', 2 / self.element_length)
        object.__setattr__(self, '_d2ξ_dx2', 4 / (self.element_length ** 2))

    @property
    def dξ_dx(self) -> float:
        """First derivative transform ∂ξ/∂x = 2/L (unitless)"""
        return self._dξ_dx

    @property
    def d2ξ_dx2(self) -> float:
        """Second derivative transform ∂²ξ/∂x² = 4/L² (1/m²)"""
        return self._d2ξ_dx2

    def natural_coordinate_form(self, ξ: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate shape functions and derivatives in natural coordinates (ξ-space).
        Only axial (DOF 0, 6) and torsion (DOF 3, 9) are non-zero.

        Parameters
        ----------
        ξ : np.ndarray
            Natural coordinates ∈ [-1, 1], shape (n_points,)

        Returns
        -------
        N : np.ndarray
            Shape function matrix [n_points, 12, 6]
        dN_dξ : np.ndarray
            First derivatives ∂N/∂ξ [n_points, 12, 6]
        d2N_dξ2 : np.ndarray
            Second derivatives (zeros for bar) [n_points, 12, 6]
        """
        ξ = np.asarray(ξ, dtype=np.float64)
        n_points = ξ.size
        ξ = ξ.reshape(-1, 1, 1)
        ξ_flat = ξ.squeeze()

        N = np.zeros((n_points, 12, 6), dtype=np.float64)
        dN_dξ = np.zeros((n_points, 12, 6), dtype=np.float64)
        d2N_dξ2 = np.zeros((n_points, 12, 6), dtype=np.float64)

        # Axial: DOF 0 (node1 u_x), 6 (node2 u_x)
        N[:, [0, 6], 0] = 0.5 * np.array([1 - ξ_flat, 1 + ξ_flat]).T
        dN_dξ[:, [0, 6], 0] = 0.5 * np.array([-1, 1])

        # Torsion: DOF 3 (node1 θ_x), 9 (node2 θ_x)
        N[:, [3, 9], 3] = 0.5 * np.array([1 - ξ_flat, 1 + ξ_flat]).T
        dN_dξ[:, [3, 9], 3] = 0.5 * np.array([-1, 1])

        return N, dN_dξ, d2N_dξ2

    def physical_coordinate_form(self, ξ: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluate shape functions and derivatives in physical coordinates (x-space).
        ∂N/∂x = (∂N/∂ξ)(∂ξ/∂x).
        """
        N, dN_dξ, d2N_dξ2 = self.natural_coordinate_form(ξ)
        dN_dx = dN_dξ * self.dξ_dx
        d2N_dx2 = d2N_dξ2 * self.d2ξ_dx2
        return N, dN_dx, d2N_dx2
