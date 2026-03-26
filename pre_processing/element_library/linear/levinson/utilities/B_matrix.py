# pre_processing/element_library/linear/levinson/utilities/B_matrix.py
"""Strain-displacement B (6, 12) per Gauss point for 2-node 3-D Levinson beam.

Voigt order is [ε_x, κ_z, κ_y, γ_xy, γ_xz, φ_x];
shear rows include α ∂²θ/∂x² terms.
Parent assembly uses `K_e += B.T @ D @ B * w_g * detJ`
with selective bending/shear rules in `linear_levinson_3D.py`.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class StrainDisplacementOperator:
    """
    Strain-displacement tensor B ∈ ℝ^{6×12} for a 2-node 3-D Levinson beam.

    B is a rank-2 tensor defined at each Gauss point such that ε = B U_e,
    where ε ∈ ℝ^6 is the generalised strain vector and U_e ∈ ℝ^{12} is the
    element displacement vector in node-major order:

        U_e = [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x², u_y², u_z², θ_x², θ_y², θ_z²]^T

    **Note on Voigt row order (Levinson)**

    The Levinson implementation uses a non-standard Voigt ordering with κ_z before κ_y,
    matching the internal implementation:

        ε = [ε_x, κ_z, κ_y, γ_xy, γ_xz, φ_x]^T
        S = [N,   M_z, M_y, V_y,  V_z,  T  ]^T

    **Kinematic equations (Levinson higher-order shear theory)**

        ε_x  = ∂u_x/∂x                            (axial extension)
        κ_z  = ∂θ_z/∂x                            (curvature about z, rotation-based)
        κ_y  = ∂θ_y/∂x                            (curvature about y, rotation-based)
        γ_xy = ∂u_y/∂x − θ_z + α ∂²θ_z/∂x²       (higher-order shear in XY; α = section property)
        γ_xz = ∂u_z/∂x − θ_y + α ∂²θ_y/∂x²       (higher-order shear in XZ)
        φ_x  = ∂θ_x/∂x                            (twist rate)

    The coefficient α encodes the cross-section geometry (e.g. α = h²/12 for a
    rectangular section of height h). Unlike Timoshenko, no shear correction factor κ
    is used in D; the higher-order correction enters kinematically through B via the
    α ∂²θ/∂x² terms, which require second derivatives of rotation shape functions.

    Parameters
    ----------
    element_length : float
        Length L of the beam element (must be > 0).
    alpha_coefficient : float, optional
        Higher-order shear coefficient α (default 0; set from section properties).

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

        L₁(ξ) = ½(1 − ξ),   dL₁/dξ = −½
        L₂(ξ) = ½(1 + ξ),   dL₂/dξ = +½
        H₁(ξ) = ¼(1 − ξ)²(2 + ξ)        dH₁/dξ = −¾(1 − ξ²)            d²H₁/dξ² = (3/2)ξ
        H₂(ξ) = (L/8)(1 − ξ)²(1 + ξ)    dH₂/dξ = (L/8)(3ξ² − 2ξ − 1)   d²H₂/dξ² = (L/8)(6ξ − 2)
        H₃(ξ) = ¼(1 + ξ)²(2 − ξ)        dH₃/dξ = ¾(1 − ξ²)             d²H₃/dξ² = −(3/2)ξ
        H₄(ξ) = −(L/8)(1 + ξ)²(1 − ξ)   dH₄/dξ = −(L/8)(1 − 2ξ − 3ξ²)  d²H₄/dξ² = (L/8)(6ξ + 2)

    **Sparsity structure of B (single Gauss point, Levinson row order)**

    ```text
    ε = B U_e,   ε ∈ ℝ^6,   U_e ∈ ℝ^{12}
    DOF cols: [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x², u_y², u_z², θ_x², θ_y², θ_z²]
    B =
    [ b1,1  0     0     0     0     0    b1,7  0     0     0      0      0   ]  # ε_x
    [ 0     0     0     0     0    b2,6   0    0     0     0      0    b2,12 ]  # κ_z
    [ 0     0     0     0    b3,5   0     0    0     0     0    b3,11    0   ]  # κ_y
    [ 0    b4,2   0     0     0    b4,6   0   b4,8   0     0      0    b4,12 ]  # γ_xy (non-zero, +α terms)
    [ 0     0    b5,3   0    b5,5   0     0    0    b5,9   0    b5,11    0   ]  # γ_xz (non-zero, +α terms)
    [ 0     0     0    b6,4   0     0     0    0     0    b6,10   0      0   ]  # φ_x
    ```

    **Non-zero entries of B in physical coordinates**

    ```text
    ε_x row (row 0) — axial:
      B[0,0]  = (dL₁/dξ)(2/L) = −1/L          (u_x, node 1)
      B[0,6]  = (dL₂/dξ)(2/L) = +1/L          (u_x, node 2)

    κ_z row (row 1) — curvature about z, first derivative of θ_z:
      B[1,5]  = (dH₂/dξ)(2/L)                 (θ_z, node 1)
      B[1,11] = (dH₄/dξ)(2/L)                 (θ_z, node 2)

    κ_y row (row 2) — curvature about y, first derivative of θ_y:
      B[2,4]  = (d(−H₂)/dξ)(2/L)              (θ_y, node 1; sign from N_θy = −H₂)
      B[2,10] = (d(−H₄)/dξ)(2/L)              (θ_y, node 2)

    γ_xy row (row 3) — higher-order shear in XY, with α ∂²θ_z/∂x²:
      B[3,1]  = (dH₁/dξ)(2/L)                          (∂u_y/∂x, node 1)
      B[3,5]  = −H₂(ξ) + α · (d²H₂/dξ²)(4/L²)         (−θ_z + α κ_z correction, node 1)
      B[3,7]  = (dH₃/dξ)(2/L)                          (∂u_y/∂x, node 2)
      B[3,11] = −H₄(ξ) + α · (d²H₄/dξ²)(4/L²)         (−θ_z + α κ_z correction, node 2)

    γ_xz row (row 4) — higher-order shear in XZ, with α ∂²θ_y/∂x²:
      B[4,2]  = (dH₁/dξ)(2/L)                          (∂u_z/∂x, node 1)
      B[4,4]  = H₂(ξ) + α · (d²(−H₂)/dξ²)(4/L²)       (−θ_y + α κ_y correction, node 1)
      B[4,8]  = (dH₃/dξ)(2/L)                          (∂u_z/∂x, node 2)
      B[4,10] = H₄(ξ) + α · (d²(−H₄)/dξ²)(4/L²)       (−θ_y + α κ_y correction, node 2)

    φ_x row (row 5) — torsion:
      B[5,3]  = (dL₁/dξ)(2/L) = −1/L          (θ_x, node 1)
      B[5,9]  = (dL₂/dξ)(2/L) = +1/L          (θ_x, node 2)
    ```

    When α = 0 (default), the γ_xy and γ_xz rows reduce to the standard Timoshenko
    form (without κ factor). The full α expression requires second derivatives of the
    Hermite rotation functions.

    Weak-form assembly: `K_e += B.T @ D @ B * w_g * detJ` with `detJ = L/2`.
    Shear rows require N (function values) and d²N/dξ²; selective integration applies.

    See Also
    --------
    linear_levinson_3D.LinearLevinsonBeamElement3D
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

        # Bending curvature κ_z = ∂θ_z/∂ξ (x–y plane, bending about z)
        B[:, 1, [5, 11]] = dN_dξ[:, [5, 11], 5]

        # Bending curvature κ_y = ∂θ_y/∂ξ (x–z plane, bending about y)
        B[:, 2, [4, 10]] = dN_dξ[:, [4, 10], 4]

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

        # κ_z = ∂θ_z/∂x (x–y plane, bending about z)
        B[:, 1, [5, 11]] = dN_dξ[:, [5, 11], 5] * self.dξ_dx

        # κ_y = ∂θ_y/∂x (x–z plane, bending about y)
        B[:, 2, [4, 10]] = dN_dξ[:, [4, 10], 4] * self.dξ_dx

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