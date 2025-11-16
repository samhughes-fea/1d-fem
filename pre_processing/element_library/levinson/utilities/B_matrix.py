# pre_processing\element_library\levinson\utilities\B_matrix.py

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class StrainDisplacementOperator:
    """
    Constructs the strain-displacement matrix `B` for a 3D Levinson beam element.

    The operator transforms first and second derivatives of shape functions with respect 
    to the natural coordinate ξ ∈ [-1, 1] into physical strain measures in x ∈ [0, L]. 

    Strain vector:
        ε = [εₓ, κ_z, κ_y, γ_xy, γ_xz, φₓ]ᵀ

    where:
        - εₓ  = ∂uₓ/∂x          (axial strain)
        - κ_z = ∂θ_y/∂x         (curvature due to bending in x–y plane) - rotation-based
        - κ_y = ∂θ_z/∂x         (curvature due to bending in x–z plane) - rotation-based
        - γ_xy = ∂u_y/∂x - θ_z + α(∂²θ_z/∂x²)  (higher-order shear in x-y plane)
        - γ_xz = ∂u_z/∂x - θ_y + α(∂²θ_y/∂x²)  (higher-order shear in x-z plane)
        - φₓ  = ∂θₓ/∂x          (torsional strain)
        
    The α coefficient for higher-order shear is typically h²/12 for rectangular sections,
    where h is the beam depth. For general sections, it may be computed from section properties.

    Coordinate mapping:
        - x(ξ) = ((1 - ξ) / 2) * x₁ + ((1 + ξ) / 2) * x₂
        - dx/dξ = L/2 ⇒ ∂ξ/∂x = 2/L
        - ∂²ξ/∂x² = 4 / L²

    Parameters
    ----------
    element_length : float
        Length `L` of the beam element in the global x-direction (must be > 0).

    Attributes
    ----------
    jacobian : float
        Determinant of the isoparametric mapping: dx/dξ = L / 2

    dξ_dx : float
        First derivative of ξ with respect to x: ∂ξ/∂x = 2 / L

    d2ξ_dx2 : float
        Second derivative of ξ with respect to x: ∂²ξ/∂x² = 4 / L²
    """

    element_length: float
    alpha_coefficient: float = 0.0  # Higher-order shear coefficient α (default 0, should be set from section properties)

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        object.__setattr__(self, '_jacobian', self.element_length / 2)
        object.__setattr__(self, '_dξ_dx', 2 / self.element_length)
        object.__setattr__(self, '_d2ξ_dx2', 4 / self.element_length ** 2)

    @property
    def jacobian(self) -> float:
        """float: Jacobian of isoparametric mapping (dx/dξ = L / 2)"""
        return self._jacobian

    @property
    def dξ_dx(self) -> float:
        """float: First derivative ∂ξ/∂x = 2 / L"""
        return self._dξ_dx

    @property
    def d2ξ_dx2(self) -> float:
        """float: Second derivative ∂²ξ/∂x² = 4 / L²"""
        return self._d2ξ_dx2

    def natural_coordinate_form(self,
                                dN_dξ: np.ndarray,
                                d2N_dξ2: np.ndarray,
                                N: np.ndarray = None) -> np.ndarray:
        """
        Construct strain-displacement matrix `B̃` in natural coordinates (ξ-space).

        Parameters
        ----------
        dN_dξ : ndarray of shape (n_gauss, 12, 6)
            First derivatives ∂N/∂ξ of shape functions with respect to ξ.
        d2N_dξ2 : ndarray of shape (n_gauss, 12, 6)
            Second derivatives ∂²N/∂ξ² of shape functions with respect to ξ.
        N : ndarray of shape (n_gauss, 12, 6), optional
            Shape functions (required for Levinson shear terms)

        Returns
        -------
        B : ndarray of shape (n_gauss, 6, 12)
            Strain-displacement matrix in ξ-space, used before transformation to physical space.

        Notes
        -----
        This form is used for symbolic verification and internal consistency checks.
        Levinson theory includes higher-order shear terms with α coefficient.
        """
        B = np.zeros((dN_dξ.shape[0], 6, 12))

        # Axial strain εₓ = ∂uₓ/∂ξ
        B[:, 0, [0, 6]] = dN_dξ[:, [0, 6], 0]

        # Bending curvature κ_z = ∂θ_y/∂ξ (Levinson: rotation-based)
        B[:, 1, [4, 10]] = dN_dξ[:, [4, 10], 4]

        # Bending curvature κ_y = ∂θ_z/∂ξ (Levinson: rotation-based)
        B[:, 2, [5, 11]] = dN_dξ[:, [5, 11], 5]

        # Higher-order shear strain: γ_xy = ∂u_y/∂ξ - θ_z + α(∂²θ_z/∂ξ²)
        if N is not None:
            B[:, 3, [1, 7]] = dN_dξ[:, [1, 7], 1]   # du_y/dx term
            B[:, 3, [5, 11]] = -N[:, [5, 11], 5]   # -θ_z term
            if self.alpha_coefficient != 0.0:
                B[:, 3, [5, 11]] += self.alpha_coefficient * d2N_dξ2[:, [5, 11], 5] * (4 / self.element_length**2)  # α(d²θ_z/dx²) term

        # Higher-order shear strain: γ_xz = ∂u_z/∂ξ - θ_y + α(∂²θ_y/∂x²)
        if N is not None:
            B[:, 4, [2, 8]] = dN_dξ[:, [2, 8], 2]   # du_z/dx term
            B[:, 4, [4, 10]] = -N[:, [4, 10], 4]   # -θ_y term
            if self.alpha_coefficient != 0.0:
                B[:, 4, [4, 10]] += self.alpha_coefficient * d2N_dξ2[:, [4, 10], 4] * (4 / self.element_length**2)  # α(d²θ_y/dx²) term

        # Torsional strain φₓ = ∂θₓ/∂ξ
        B[:, 5, [3, 9]] = dN_dξ[:, [3, 9], 3]

        return B

    def physical_coordinate_form(self,
                                 dN_dξ: np.ndarray,
                                 d2N_dξ2: np.ndarray,
                                 N: np.ndarray = None) -> np.ndarray:
        """
        Construct strain-displacement matrix `B` in physical coordinates (x-space).

        Parameters
        ----------
        dN_dξ : ndarray of shape (n_gauss, 12, 6)
            First derivatives ∂N/∂ξ of shape functions with respect to ξ.
        d2N_dξ2 : ndarray of shape (n_gauss, 12, 6)
            Second derivatives ∂²N/∂ξ² of shape functions with respect to ξ.
        N : ndarray of shape (n_gauss, 12, 6), optional
            Shape functions (required for Levinson shear terms)

        Returns
        -------
        B : ndarray of shape (n_gauss, 6, 12)
            Physical strain-displacement matrix such that ε = B @ u_e

        Notes
        -----
        - The coordinate transformation is handled internally.
        - Levinson theory includes higher-order shear terms with α coefficient.
        """
        B = np.zeros((dN_dξ.shape[0], 6, 12))

        # εₓ = ∂uₓ/∂x = ∂uₓ/∂ξ * ∂ξ/∂x
        B[:, 0, [0, 6]] = dN_dξ[:, [0, 6], 0] * self.dξ_dx

        # κ_z = ∂θ_y/∂x (Levinson: rotation-based)
        B[:, 1, [4, 10]] = dN_dξ[:, [4, 10], 4] * self.dξ_dx

        # κ_y = ∂θ_z/∂x (Levinson: rotation-based)
        B[:, 2, [5, 11]] = dN_dξ[:, [5, 11], 5] * self.dξ_dx

        # Higher-order shear strain: γ_xy = ∂u_y/∂x - θ_z + α(∂²θ_z/∂x²)
        if N is not None:
            B[:, 3, [1, 7]] = dN_dξ[:, [1, 7], 1] * self.dξ_dx      # du_y/dx term
            B[:, 3, [5, 11]] = -N[:, [5, 11], 5]                    # -θ_z term (no coordinate transform)
            if self.alpha_coefficient != 0.0:
                B[:, 3, [5, 11]] += self.alpha_coefficient * d2N_dξ2[:, [5, 11], 5] * self.d2ξ_dx2  # α(d²θ_z/dx²) term

        # Higher-order shear strain: γ_xz = ∂u_z/∂x - θ_y + α(∂²θ_y/∂x²)
        if N is not None:
            B[:, 4, [2, 8]] = dN_dξ[:, [2, 8], 2] * self.dξ_dx      # du_z/dx term
            B[:, 4, [4, 10]] = -N[:, [4, 10], 4]                    # -θ_y term (no coordinate transform)
            if self.alpha_coefficient != 0.0:
                B[:, 4, [4, 10]] += self.alpha_coefficient * d2N_dξ2[:, [4, 10], 4] * self.d2ξ_dx2  # α(d²θ_y/dx²) term

        # φₓ = ∂θₓ/∂x
        B[:, 5, [3, 9]] = dN_dξ[:, [3, 9], 3] * self.dξ_dx

        return B

    def verify_coordinate_transforms(self, tol: float = 1e-12) -> Tuple[bool, str]:
        """
        Check analytical coordinate transform identities within tolerance.

        Parameters
        ----------
        tol : float, optional
            Numerical tolerance for validation. Default is 1e-12.

        Returns
        -------
        Tuple[bool, str]
            (True, message) if valid; otherwise (False, error message).
        """
        checks = [
            ("Jacobian", abs(self.jacobian - self.element_length / 2)),
            ("First derivative", abs(self.dξ_dx - 2 / self.element_length)),
            ("Second derivative", abs(self.d2ξ_dx2 - 4 / self.element_length ** 2))
        ]
        for name, error in checks:
            if error > tol:
                return False, f"{name} transform error: {error:.2e} > {tol}"
        return True, "All coordinate transforms valid"