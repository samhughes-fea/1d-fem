# pre_processing/element_library/linear/timoshenko/utilities/B_matrix.py
"""Strain-displacement B (6, 12) per Gauss point for 2-node 3-D Timoshenko beam.

ε = B U_e with Voigt order [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x].
Shear terms are γ_xy = ∂u_y/∂x - θ_z and γ_xz = ∂u_z/∂x - θ_y.
Parent assembly uses `K_e += B.T @ D @ B * w_g * detJ`.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class StrainDisplacementOperator:
    """
    Strain-displacement tensor B ∈ ℝ^{6×12} for a 2-node 3-D Timoshenko beam.

    B is a rank-2 tensor defined at each Gauss point such that ε = B U_e,
    where ε ∈ ℝ^6 is the generalised strain vector and U_e ∈ ℝ^{12} is the
    element displacement vector in node-major order:

        U_e = [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x², u_y², u_z², θ_x², θ_y², θ_z²]^T

    **Kinematic equations (Timoshenko first-order shear theory)**

    The six Voigt strain components are:

        ε_x  = ∂u_x/∂x               (axial extension)
        κ_y  = ∂θ_y/∂x               (curvature about y, rotation-based)
        κ_z  = ∂θ_z/∂x               (curvature about z, rotation-based)
        γ_xy = ∂u_y/∂x − θ_z         (shear strain in XY plane; non-zero in Timoshenko)
        γ_xz = ∂u_z/∂x − θ_y         (shear strain in XZ plane; non-zero in Timoshenko)
        φ_x  = ∂θ_x/∂x               (twist rate)

    Unlike Euler-Bernoulli, curvatures are first derivatives of independent rotations
    (not second derivatives of transverse displacements). Shear strains γ_xy and γ_xz
    couple displacement and rotation DOFs and enter the constitutive relation S = D ε
    with D[3,3] = D[4,4] = κ·G·A (shear correction factor κ explicit).

    Parameters
    ----------
    element_length : float
        Length L of the beam element (must be > 0).

    Attributes
    ----------
    jacobian : float
        Jacobian of isoparametric mapping, ∂x/∂ξ = L/2.
    dξ_dx : float
        First coordinate transform factor, ∂ξ/∂x = 2/L.
    d2ξ_dx2 : float
        Second coordinate transform factor, ∂²ξ/∂x² = 4/L².

    Notes
    -----
    **Isoparametric mapping and shape function basis**

    x ∈ [0, L] maps to ξ ∈ [−1, 1] via x(ξ) = L(1 + ξ)/2, giving
    ∂ξ/∂x = 2/L and ∂²ξ/∂x² = 4/L².

    The same Hermite cubics and linear Lagrange functions used in EB apply here;
    however, B uses only first derivatives of rotation functions for bending rows
    (not second derivatives of displacement functions):

        L₁(ξ) = ½(1 − ξ),   dL₁/dξ = −½
        L₂(ξ) = ½(1 + ξ),   dL₂/dξ = +½
        H₁(ξ) = ¼(1 − ξ)²(2 + ξ)        dH₁/dξ = −¾(1 − ξ²)
        H₂(ξ) = (L/8)(1 − ξ)²(1 + ξ)    dH₂/dξ = (L/8)(3ξ² − 2ξ − 1)
        H₃(ξ) = ¼(1 + ξ)²(2 − ξ)        dH₃/dξ = ¾(1 − ξ²)
        H₄(ξ) = −(L/8)(1 + ξ)²(1 − ξ)   dH₄/dξ = −(L/8)(1 − 2ξ − 3ξ²)

    **Sparsity structure of B (single Gauss point)**

    ```text
    ε = B U_e,   ε ∈ ℝ^6,   U_e ∈ ℝ^{12}
    DOF cols: [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x², u_y², u_z², θ_x², θ_y², θ_z²]
    B =
    [ b1,1  0     0     0     0     0    b1,7  0     0     0      0      0   ]  # ε_x
    [ 0     0     0     0    b2,5   0     0    0     0     0    b2,11    0   ]  # κ_y
    [ 0     0     0     0     0    b3,6   0    0     0     0      0    b3,12 ]  # κ_z
    [ 0    b4,2   0     0     0    b4,6   0   b4,8   0     0      0    b4,12 ]  # γ_xy (non-zero)
    [ 0     0    b5,3   0    b5,5   0     0    0    b5,9   0    b5,11    0   ]  # γ_xz (non-zero)
    [ 0     0     0    b6,4   0     0     0    0     0    b6,10   0      0   ]  # φ_x
    ```

    **Non-zero entries of B in physical coordinates**

    ```text
    ε_x row (row 0) — axial:
      B[0,0]  = (dL₁/dξ)(2/L) = −1/L          (u_x, node 1)
      B[0,6]  = (dL₂/dξ)(2/L) = +1/L          (u_x, node 2)

    κ_y row (row 1) — bending about y, first derivative of θ_y:
      B[1,4]  = (d(−H₂)/dξ)(2/L) = −(L/8)(3ξ² − 2ξ − 1) · (2/L)   (θ_y, node 1)
      B[1,10] = (d(−H₄)/dξ)(2/L) = (L/8)(1 − 2ξ − 3ξ²) · (2/L)    (θ_y, node 2)

    κ_z row (row 2) — bending about z, first derivative of θ_z:
      B[2,5]  = (dH₂/dξ)(2/L) = (L/8)(3ξ² − 2ξ − 1) · (2/L)       (θ_z, node 1)
      B[2,11] = (dH₄/dξ)(2/L) = −(L/8)(1 − 2ξ − 3ξ²) · (2/L)      (θ_z, node 2)

    γ_xy row (row 3) — shear in XY plane:
      B[3,1]  = (dH₁/dξ)(2/L) = −¾(1 − ξ²) · (2/L)     (∂u_y/∂x, node 1)
      B[3,5]  = −H₂(ξ) = −(L/8)(1 − ξ)²(1 + ξ)         (−θ_z, node 1)
      B[3,7]  = (dH₃/dξ)(2/L) =  ¾(1 − ξ²) · (2/L)     (∂u_y/∂x, node 2)
      B[3,11] = −H₄(ξ) = (L/8)(1 + ξ)²(1 − ξ)          (−θ_z, node 2)

    γ_xz row (row 4) — shear in XZ plane:
      B[4,2]  = (dH₁/dξ)(2/L) = −¾(1 − ξ²) · (2/L)     (∂u_z/∂x, node 1)
      B[4,4]  = −N[4,4] = H₂(ξ)                         (−θ_y, node 1; N_θy = −H₂, so −N_θy = H₂)
      B[4,8]  = (dH₃/dξ)(2/L) =  ¾(1 − ξ²) · (2/L)     (∂u_z/∂x, node 2)
      B[4,10] = −N[10,4] = H₄(ξ)                        (−θ_y, node 2; N_θy = −H₄, so −N_θy = H₄)

    φ_x row (row 5) — torsion:
      B[5,3]  = (dL₁/dξ)(2/L) = −1/L          (θ_x, node 1)
      B[5,9]  = (dL₂/dξ)(2/L) = +1/L          (θ_x, node 2)
    ```

    Weak-form assembly: `K_e += B.T @ D @ B * w_g * detJ` with `detJ = L/2`.
    `physical_coordinate_form` applies the chain-rule factors; `natural_coordinate_form`
    omits them. Shear rows require N (function values), not only dN_dξ.

    See Also
    --------
    linear_timoshenko_3D.LinearTimoshenkoBeamElement3D
    """

    element_length: float

    def __post_init__(self):
        if self.element_length <= 0:
            raise ValueError(f"Element length must be positive, got {self.element_length}")
        object.__setattr__(self, '_jacobian', self.element_length / 2)
        object.__setattr__(self, '_dξ_dx', 2 / self.element_length)
        object.__setattr__(self, '_d2ξ_dx2', 4 / self.element_length ** 2)

    @property
    def jacobian(self) -> float:
        """float: Jacobian of isoparametric mapping (dx/dξ = L/2)"""
        return self._jacobian

    @property
    def dξ_dx(self) -> float:
        """float: First derivative ∂ξ/∂x = 2/L"""
        return self._dξ_dx

    @property
    def d2ξ_dx2(self) -> float:
        """float: Second derivative ∂²ξ/∂x² = 4/L²"""
        return self._d2ξ_dx2

    def natural_coordinate_form(self,
                                dN_dξ: np.ndarray,
                                d2N_dξ2: np.ndarray,
                                N: np.ndarray = None) -> np.ndarray:
        """
        Construct strain-displacement matrix `B̃` in natural coordinates (ξ-space).

        Parameters
        ----------
        dN_dξ : np.ndarray (n_gauss, 12, 6)
            First derivatives of shape functions
        d2N_dξ2 : np.ndarray (n_gauss, 12, 6)
            Second derivatives of shape functions (not used for Timoshenko bending)
        N : np.ndarray (n_gauss, 12, 6), optional
            Shape functions (required for Timoshenko shear terms)

        Returns
        -------
        B : np.ndarray (n_gauss, 6, 12)
            Strain-displacement matrix in ξ-space
        """
        n_gauss = dN_dξ.shape[0]
        B = np.zeros((n_gauss, 6, 12))

        # Axial strain: ε_x = ∂u_x/∂ξ
        B[:, 0, [0, 6]] = dN_dξ[:, [0, 6], 0]       # u_x

        # Bending about y-axis: κ_y = ∂θ_y/∂ξ (Timoshenko: rotation-based, not displacement-based)
        B[:, 1, [4, 10]] = dN_dξ[:, [4, 10], 4]     # θ_y

        # Bending about z-axis: κ_z = ∂θ_z/∂ξ (Timoshenko: rotation-based, not displacement-based)
        B[:, 2, [5, 11]] = dN_dξ[:, [5, 11], 5]     # θ_z

        # Shear strain: γ_xy = ∂u_y/∂ξ - θ_z (Timoshenko includes shear)
        if N is not None:
            B[:, 3, [1, 7]] = dN_dξ[:, [1, 7], 1]   # du_y/dx term
            B[:, 3, [5, 11]] = -N[:, [5, 11], 5]   # -θ_z term

        # Shear strain: γ_xz = ∂u_z/∂ξ - θ_y (Timoshenko includes shear)
        if N is not None:
            B[:, 4, [2, 8]] = dN_dξ[:, [2, 8], 2]   # du_z/dx term
            B[:, 4, [4, 10]] = -N[:, [4, 10], 4]   # -θ_y term

        # Torsional strain: φ_x = ∂θ_x/∂ξ
        B[:, 5, [3, 9]] = dN_dξ[:, [3, 9], 3]       # θ_x

        return B

    def physical_coordinate_form(self,
                                 dN_dξ: np.ndarray,
                                 d2N_dξ2: np.ndarray,
                                 N: np.ndarray = None) -> np.ndarray:
        """
        Construct strain-displacement matrix `B` in physical coordinates (x-space).

        Parameters
        ----------
        dN_dξ : np.ndarray (n_gauss, 12, 6)
            First derivatives of shape functions
        d2N_dξ2 : np.ndarray (n_gauss, 12, 6)
            Second derivatives of shape functions (not used for Timoshenko bending)
        N : np.ndarray (n_gauss, 12, 6), optional
            Shape functions (required for Timoshenko shear terms)

        Returns
        -------
        B : np.ndarray (n_gauss, 6, 12)
            Physical strain-displacement matrix (ε = B @ u_e)
        """
        n_gauss = dN_dξ.shape[0]
        B = np.zeros((n_gauss, 6, 12))

        # Axial strain: ε_x = ∂u_x/∂x
        B[:, 0, [0, 6]] = dN_dξ[:, [0, 6], 0] * self.dξ_dx           # u_x

        # Bending about y-axis: κ_y = ∂θ_y/∂x (Timoshenko: rotation-based)
        B[:, 1, [4, 10]] = dN_dξ[:, [4, 10], 4] * self.dξ_dx         # θ_y

        # Bending about z-axis: κ_z = ∂θ_z/∂x (Timoshenko: rotation-based)
        B[:, 2, [5, 11]] = dN_dξ[:, [5, 11], 5] * self.dξ_dx         # θ_z

        # Shear strain: γ_xy = ∂u_y/∂x - θ_z (Timoshenko includes shear)
        if N is not None:
            B[:, 3, [1, 7]] = dN_dξ[:, [1, 7], 1] * self.dξ_dx      # du_y/dx term
            B[:, 3, [5, 11]] = -N[:, [5, 11], 5]                    # -θ_z term (no coordinate transform needed)

        # Shear strain: γ_xz = ∂u_z/∂x - θ_y (Timoshenko includes shear)
        if N is not None:
            B[:, 4, [2, 8]] = dN_dξ[:, [2, 8], 2] * self.dξ_dx      # du_z/dx term
            B[:, 4, [4, 10]] = -N[:, [4, 10], 4]                    # -θ_y term (no coordinate transform needed)

        # Torsional strain: φ_x = ∂θ_x/∂x
        B[:, 5, [3, 9]] = dN_dξ[:, [3, 9], 3] * self.dξ_dx           # θ_x

        return B

    def verify_coordinate_transforms(self, tol: float = 1e-12) -> Tuple[bool, str]:
        """
        Validate coordinate transformation parameters.

        Parameters
        ----------
        tol : float
            Numerical tolerance for validation

        Returns
        -------
        Tuple[bool, str]
            Validation status and message
        """
        checks = [
            ("Jacobian", abs(self.jacobian - self.element_length/2)),
            ("∂ξ/∂x", abs(self.dξ_dx - 2/self.element_length)),
            ("∂²ξ/∂x²", abs(self.d2ξ_dx2 - 4/self.element_length**2))
        ]
        for name, error in checks:
            if error > tol:
                return False, f"{name} error: {error:.2e} > {tol}"
        return True, "All transforms valid"