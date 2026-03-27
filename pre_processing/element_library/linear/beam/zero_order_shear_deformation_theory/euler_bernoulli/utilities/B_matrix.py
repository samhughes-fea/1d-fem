# pre_processing/element_library/linear/euler_bernoulli/utilities/B_matrix.py
"""Strain-displacement B (6, 12) per Gauss point for 2-node 3-D Euler-Bernoulli beam.

ε = B U_e with ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]; γ_xy = γ_xz = 0.
Parent element uses `K_e += B.T @ D @ B * w_g * detJ` with `detJ = L/2`.
Voigt order follows `FORMULATION_DOCSTRING_STANDARDS.md`.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class StrainDisplacementOperator:
    """
    Strain-displacement tensor B ∈ ℝ^{6×12} for a 2-node 3-D Euler-Bernoulli beam.

    B is a rank-2 tensor defined at each Gauss point such that ε = B U_e,
    where ε ∈ ℝ^6 is the generalised strain vector and U_e ∈ ℝ^{12} is the
    element displacement vector in node-major order:

        U_e = [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x², u_y², u_z², θ_x², θ_y², θ_z²]^T

    **Kinematic equations (Euler-Bernoulli theory)**

    The six Voigt strain components and their physical interpretations are:

        ε_x  = ∂u_x/∂x               (axial extension)
        κ_y  = ∂²u_z/∂x²             (curvature about y; EB: θ_y = −∂u_z/∂x)
        κ_z  = ∂²u_y/∂x²             (curvature about z; EB: θ_z = ∂u_y/∂x)
        γ_xy = 0                      (shear-inextensibility, Euler-Bernoulli)
        γ_xz = 0                      (shear-inextensibility, Euler-Bernoulli)
        φ_x  = ∂θ_x/∂x               (twist rate, St. Venant torsion)

    The constitutive relation is S = D ε with S = [N, M_y, M_z, V_y, V_z, T]^T.
    The EB shear resultants V_y and V_z are zero from the constitutive path; shear
    forces are recovered from equilibrium V = dM/dx, not from D ε.

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
    **Isoparametric mapping**

    The physical coordinate x ∈ [0, L] maps to ξ ∈ [−1, 1] via:

        x(ξ) = L(1 + ξ)/2,   ∂x/∂ξ = L/2,   ∂ξ/∂x = 2/L,   ∂²ξ/∂x² = 4/L².

    Chain rule: ∂N/∂x = (∂N/∂ξ)(2/L),  ∂²N/∂x² = (∂²N/∂ξ²)(4/L²).

    **Shape function basis**

    Linear Lagrange polynomials on axial (u_x) and torsion (θ_x) channels:

        L₁(ξ) = ½(1 − ξ),   dL₁/dξ = −½
        L₂(ξ) = ½(1 + ξ),   dL₂/dξ = +½

    Hermite cubic polynomials on bending channels (standard C¹ beam pair per plane):

        H₁(ξ) = ¼(1 − ξ)²(2 + ξ)        dH₁/dξ = −¾(1 − ξ²)            d²H₁/dξ² = (3/2)ξ
        H₂(ξ) = (L/8)(1 − ξ)²(1 + ξ)    dH₂/dξ = (L/8)(3ξ² − 2ξ − 1)   d²H₂/dξ² = (L/8)(6ξ − 2)
        H₃(ξ) = ¼(1 + ξ)²(2 − ξ)        dH₃/dξ = ¾(1 − ξ²)             d²H₃/dξ² = −(3/2)ξ
        H₄(ξ) = −(L/8)(1 + ξ)²(1 − ξ)   dH₄/dξ = −(L/8)(1 − 2ξ − 3ξ²)  d²H₄/dξ² = (L/8)(6ξ + 2)

    H₁, H₃ are displacement functions; H₂, H₄ are rotation functions (scaled by L
    to maintain consistent units and C¹ continuity across elements).

    **Sparsity structure of B (single Gauss point)**

    ```text
    ε = B U_e,   ε ∈ ℝ^6,   U_e ∈ ℝ^{12}
    DOF cols: [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x², u_y², u_z², θ_x², θ_y², θ_z²]
    B =
    [ b1,1  0     0     0     0     0    b1,7  0     0     0      0      0   ]  # ε_x
    [ 0     0    b2,3   0    b2,5   0     0    0    b2,9   0    b2,11    0   ]  # κ_y
    [ 0    b3,2   0     0     0    b3,6   0   b3,8   0     0      0    b3,12 ]  # κ_z
    [ 0     0     0     0     0     0     0    0     0     0      0      0   ]  # γ_xy = 0
    [ 0     0     0     0     0     0     0    0     0     0      0      0   ]  # γ_xz = 0
    [ 0     0     0    b6,4   0     0     0    0     0    b6,10   0      0   ]  # φ_x
    ```

    **Non-zero entries of B in physical coordinates**

    ```text
    ε_x row (row 0) — axial:
      B[0,0]  = (dL₁/dξ)(2/L) = −1/L          (u_x, node 1)
      B[0,6]  = (dL₂/dξ)(2/L) = +1/L          (u_x, node 2)

    κ_y row (row 1) — bending about y, Hermite second derivatives of XZ-plane functions:
      B[1,2]  = (d²H₁/dξ²)(4/L²) =  6ξ/L²          (u_z, node 1)
      B[1,4]  = (d²(−H₂)/dξ²)(4/L²) = (1 − 3ξ)/L   (θ_y, node 1; sign from EB rotation convention)
      B[1,8]  = (d²H₃/dξ²)(4/L²) = −6ξ/L²          (u_z, node 2)
      B[1,10] = (d²(−H₄)/dξ²)(4/L²) = −(1 + 3ξ)/L  (θ_y, node 2)

    κ_z row (row 2) — bending about z, Hermite second derivatives of XY-plane functions:
      B[2,1]  = (d²H₁/dξ²)(4/L²) =  6ξ/L²          (u_y, node 1)
      B[2,5]  = (d²H₂/dξ²)(4/L²)  = (3ξ − 1)/L     (θ_z, node 1)
      B[2,7]  = (d²H₃/dξ²)(4/L²) = −6ξ/L²          (u_y, node 2)
      B[2,11] = (d²H₄/dξ²)(4/L²)  = (1 + 3ξ)/L     (θ_z, node 2)

    γ_xy row (row 3): B[3,j] = 0 for all j   (EB shear-inextensibility)
    γ_xz row (row 4): B[4,j] = 0 for all j   (EB shear-inextensibility)

    φ_x row (row 5) — torsion:
      B[5,3]  = (dL₁/dξ)(2/L) = −1/L          (θ_x, node 1)
      B[5,9]  = (dL₂/dξ)(2/L) = +1/L          (θ_x, node 2)
    ```

    The natural-coordinate form B̃ omits the (2/L) and (4/L²) factors;
    `physical_coordinate_form` applies them. Shear rows remain zero in both forms.

    Weak-form assembly: `K_e += B.T @ D @ B * w_g * detJ` with `detJ = L/2`.
    Shear-deformable theories (Timoshenko, Levinson) populate rows 3 and 4 with
    non-zero kinematics and add GA or κGA to the constitutive diagonal.

    See Also
    --------
    linear_euler_bernoulli_3D.LinearEulerBernoulliBeamElement3D
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
                                d2N_dξ2: np.ndarray) -> np.ndarray:
        """
        Construct strain-displacement matrix `B̃` in natural coordinates (ξ-space).

        Parameters
        ----------
        dN_dξ : np.ndarray (n_gauss, 12, 6)
            First derivatives of shape functions
        d2N_dξ2 : np.ndarray (n_gauss, 12, 6)
            Second derivatives of shape functions

        Returns
        -------
        B : np.ndarray (n_gauss, 6, 12)
            Strain-displacement matrix in ξ-space
        """
        n_gauss = dN_dξ.shape[0]
        B = np.zeros((n_gauss, 6, 12))

        # Axial strain: ε_x = ∂u_x/∂ξ
        B[:, 0, [0, 6]] = dN_dξ[:, [0, 6], 0]       # u_x

        # Bending about y-axis: κ_y = ∂²u_z/∂ξ²
        B[:, 1, [2, 8]] = d2N_dξ2[:, [2, 8], 2]     # u_z
        B[:, 1, [4, 10]] = d2N_dξ2[:, [4, 10], 4]   # θ_y

        # Bending about z-axis: κ_z = ∂²u_y/∂ξ²
        B[:, 2, [1, 7]] = d2N_dξ2[:, [1, 7], 1]     # u_y
        B[:, 2, [5, 11]] = d2N_dξ2[:, [5, 11], 5]   # θ_z

        # Torsional strain: φ_x = ∂θ_x/∂ξ
        B[:, 5, [3, 9]] = dN_dξ[:, [3, 9], 3]       # θ_x

        return B

    def physical_coordinate_form(self,
                                 dN_dξ: np.ndarray,
                                 d2N_dξ2: np.ndarray) -> np.ndarray:
        """
        Construct strain-displacement matrix `B` in physical coordinates (x-space).

        Parameters
        ----------
        dN_dξ : np.ndarray (n_gauss, 12, 6)
            First derivatives of shape functions
        d2N_dξ2 : np.ndarray (n_gauss, 12, 6)
            Second derivatives of shape functions

        Returns
        -------
        B : np.ndarray (n_gauss, 6, 12)
            Physical strain-displacement matrix (ε = B @ u_e)
        """
        n_gauss = dN_dξ.shape[0]
        B = np.zeros((n_gauss, 6, 12))

        # Axial strain: ε_x = ∂u_x/∂x
        B[:, 0, [0, 6]] = dN_dξ[:, [0, 6], 0] * self.dξ_dx           # u_x

        # Bending about y-axis: κ_y = ∂²u_z/∂ξ²
        B[:, 1, [2, 8]] = d2N_dξ2[:, [2, 8], 2] * self.d2ξ_dx2       # u_z
        B[:, 1, [4, 10]] = d2N_dξ2[:, [4, 10], 4] * self.d2ξ_dx2     # θ_y
 
        # Bending about z-axis: κ_z = ∂²u_y/∂ξ²
        B[:, 2, [1, 7]] = d2N_dξ2[:, [1, 7], 1] * self.d2ξ_dx2       # u_y
        B[:, 2, [5, 11]] = d2N_dξ2[:, [5, 11], 5] * self.d2ξ_dx2     # θ_z

        # Torsional strain: φ_x = ∂θ_x/∂ξ
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