# pre_processing/element_library/linear/truss/utilities/B_matrix.py
"""Strain-displacement matrix B for 2-node 3-D Truss (axial + transverse + torsion). ε = B @ u_e; B shape (3, 12) per GP."""

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class StrainDisplacementOperator:
    """
    Strain–displacement matrix **B** for a 2-node 3-D Truss element (axial + transverse + torsion).

    Strain vector: ε = [ ε_axial  γ_transverse  φ_torsion ]ᵀ
    - ε_axial = ∂(axial displacement)/∂x
    - γ_transverse = ∂(transverse displacement)/∂x
    - φ_torsion = ∂θ_x/∂x

    B is constant in ξ.

    Parameters
    ----------
    element_length : float
        Length L of the element (must be > 0).
    axial : np.ndarray, shape (3,)
        Unit vector along the element.
    transverse : np.ndarray, shape (3,)
        Unit vector for transverse direction (local y).
    """

    element_length: float
    axial: np.ndarray
    transverse: np.ndarray

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        axial = np.asarray(self.axial, dtype=np.float64)
        transverse = np.asarray(self.transverse, dtype=np.float64)
        if axial.shape != (3,) or np.abs(np.linalg.norm(axial) - 1.0) > 1e-9:
            raise ValueError("axial must be a unit vector of shape (3,)")
        if transverse.shape != (3,) or np.abs(np.linalg.norm(transverse) - 1.0) > 1e-9:
            raise ValueError("transverse must be a unit vector of shape (3,)")
        object.__setattr__(self, '_axial', axial)
        object.__setattr__(self, '_transverse', transverse)
        object.__setattr__(self, '_jacobian', self.element_length / 2)
        object.__setattr__(self, '_dξ_dx', 2 / self.element_length)

    @property
    def jacobian(self) -> float:
        """Jacobian dx/dξ = L/2"""
        return self._jacobian

    @property
    def dξ_dx(self) -> float:
        """∂ξ/∂x = 2/L"""
        return self._dξ_dx

    def natural_coordinate_form(self, dN_dξ: np.ndarray) -> np.ndarray:
        """
        Construct B in natural coordinates. Returns (n_gauss, 3, 12).
        """
        n_gauss = dN_dξ.shape[0]
        B = np.zeros((n_gauss, 3, 12), dtype=np.float64)
        L = self.element_length
        # Row 0: ε_axial
        B[:, 0, 0:3] = -1.0 / L * self._axial
        B[:, 0, 6:9] = 1.0 / L * self._axial
        # Row 1: γ_transverse (transverse component of displacement derivative)
        B[:, 1, 0:3] = -1.0 / L * self._transverse
        B[:, 1, 6:9] = 1.0 / L * self._transverse
        # Row 2: φ_torsion
        B[:, 2, 3] = -1.0 / L
        B[:, 2, 9] = 1.0 / L
        return B

    def physical_coordinate_form(self, dN_dξ: np.ndarray) -> np.ndarray:
        """Construct B in physical coordinates. Same as natural for truss."""
        return self.natural_coordinate_form(dN_dξ)

    def verify_coordinate_transforms(self, tol: float = 1e-12) -> Tuple[bool, str]:
        """Validate jacobian and dξ_dx."""
        if abs(self.jacobian - self.element_length / 2) > tol:
            return False, f"Jacobian error > {tol}"
        if abs(self.dξ_dx - 2 / self.element_length) > tol:
            return False, f"dξ_dx error > {tol}"
        return True, "All transforms valid"
