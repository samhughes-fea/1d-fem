# pre_processing\element_library\bar\utilities\B_matrix.py

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class StrainDisplacementOperator:
    """
    Strain–displacement matrix **B** for a 2-node 3-D Bar element (axial + torsion).

    Strain vector: ε = [ ε_axial  φ_torsion ]ᵀ
    - ε_axial = ∂(axial displacement)/∂x
    - φ_torsion = ∂θ_x/∂x

    B is constant in ξ (linear shape functions). One Gauss point suffices for exact integration.

    Parameters
    ----------
    element_length : float
        Length L of the element (must be > 0).
    axial : np.ndarray, shape (3,)
        Unit vector (direction cosines) along the element.

    Attributes
    ----------
    jacobian : float
        Jacobian dx/dξ = L/2
    dξ_dx : float
        ∂ξ/∂x = 2/L
    """

    element_length: float
    axial: np.ndarray

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        axial = np.asarray(self.axial, dtype=np.float64)
        if axial.shape != (3,) or np.abs(np.linalg.norm(axial) - 1.0) > 1e-9:
            raise ValueError("axial must be a unit vector of shape (3,)")
        object.__setattr__(self, '_axial', axial)
        object.__setattr__(self, '_jacobian', self.element_length / 2)
        object.__setattr__(self, '_dξ_dx', 2 / self.element_length)

    @property
    def jacobian(self) -> float:
        """Jacobian of isoparametric mapping (dx/dξ = L/2)"""
        return self._jacobian

    @property
    def dξ_dx(self) -> float:
        """First derivative ∂ξ/∂x = 2/L"""
        return self._dξ_dx

    def natural_coordinate_form(self, dN_dξ: np.ndarray) -> np.ndarray:
        """
        Construct strain-displacement matrix B in natural coordinates.
        For bar, B is constant; dN_dξ is used only for shape compatibility.
        Returns B of shape (n_gauss, 2, 12) with n_gauss = dN_dξ.shape[0].
        """
        n_gauss = dN_dξ.shape[0]
        B = np.zeros((n_gauss, 2, 12), dtype=np.float64)
        L = self.element_length
        cx, cy, cz = self._axial
        # Row 0: ε_axial = (1/L)(-u1_axial + u2_axial) in global
        B[:, 0, 0:3] = -1.0 / L * self._axial
        B[:, 0, 6:9] = 1.0 / L * self._axial
        # Row 1: φ_torsion = (1/L)(-θx1 + θx2)
        B[:, 1, 3] = -1.0 / L
        B[:, 1, 9] = 1.0 / L
        return B

    def physical_coordinate_form(self, dN_dξ: np.ndarray) -> np.ndarray:
        """
        Construct strain-displacement matrix B in physical coordinates.
        For bar, same as natural (B already encodes ∂/∂x via 1/L).
        Returns (n_gauss, 2, 12).
        """
        return self.natural_coordinate_form(dN_dξ)

    def verify_coordinate_transforms(self, tol: float = 1e-12) -> Tuple[bool, str]:
        """Validate jacobian and dξ_dx."""
        if abs(self.jacobian - self.element_length / 2) > tol:
            return False, f"Jacobian error > {tol}"
        if abs(self.dξ_dx - 2 / self.element_length) > tol:
            return False, f"dξ_dx error > {tol}"
        return True, "All transforms valid"
