# pre_processing/element_library/linear/truss/utilities/shape_functions.py
"""
Truss shape functions: ``natural_coordinate_form`` → **``N``**, **``dN_dξ``**, **``d2N_dξ2``** with shape **(n_gp, 12, 6)**.
Weak-form distributed loads: **``F_dist = ∫ Nᵀ q |J| dξ``**, **``|J| = L/2``**.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class ShapeFunctionOperator:
    """
    Operator for evaluating 3D Truss element shape functions and their derivatives.
    Truss has axial, transverse, and torsion; all use linear Lagrange in ξ.

    Linear Lagrange for: axial (DOF 0, 6), transverse (DOF 1, 2, 7, 8), torsion (DOF 3, 9).

    Parameters
    ----------
    element_length : float
        Physical length of element (x ∈ [0,L], L > 0)

    Attributes
    ----------
    dξ_dx : float
        First derivative transform (∂ξ/∂x = 2/L)
    d2ξ_dx2 : float
        Second derivative transform (∂²ξ/∂x² = 4/L²)
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
        Linear Lagrange for u_x (0,6), u_y (1,7), u_z (2,8), θ_x (3,9).

        Returns
        -------
        N : np.ndarray [n_points, 12, 6]
        dN_dξ : np.ndarray [n_points, 12, 6]
        d2N_dξ2 : np.ndarray [n_points, 12, 6] (zeros)
        """
        ξ = np.asarray(ξ, dtype=np.float64)
        n_points = ξ.size
        ξ = ξ.reshape(-1, 1, 1)
        ξ_flat = ξ.squeeze()

        N = np.zeros((n_points, 12, 6), dtype=np.float64)
        dN_dξ = np.zeros((n_points, 12, 6), dtype=np.float64)
        d2N_dξ2 = np.zeros((n_points, 12, 6), dtype=np.float64)

        # Axial u_x: DOF 0, 6
        N[:, [0, 6], 0] = 0.5 * np.array([1 - ξ_flat, 1 + ξ_flat]).T
        dN_dξ[:, [0, 6], 0] = 0.5 * np.array([-1, 1])

        # Transverse u_y, u_z: DOF 1,7 and 2,8
        N[:, [1, 7], 1] = 0.5 * np.array([1 - ξ_flat, 1 + ξ_flat]).T
        dN_dξ[:, [1, 7], 1] = 0.5 * np.array([-1, 1])
        N[:, [2, 8], 2] = 0.5 * np.array([1 - ξ_flat, 1 + ξ_flat]).T
        dN_dξ[:, [2, 8], 2] = 0.5 * np.array([-1, 1])

        # Torsion θ_x: DOF 3, 9
        N[:, [3, 9], 3] = 0.5 * np.array([1 - ξ_flat, 1 + ξ_flat]).T
        dN_dξ[:, [3, 9], 3] = 0.5 * np.array([-1, 1])

        return N, dN_dξ, d2N_dξ2

    def physical_coordinate_form(self, ξ: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate shape functions and derivatives in physical coordinates."""
        N, dN_dξ, d2N_dξ2 = self.natural_coordinate_form(ξ)
        dN_dx = dN_dξ * self.dξ_dx
        d2N_dx2 = d2N_dξ2 * self.d2ξ_dx2
        return N, dN_dx, d2N_dx2
