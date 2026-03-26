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
    Strain-displacement ``B`` (6, 12) per Gauss point for a 2-node 3-D Levinson beam.

    Voigt ε = B U_e with rows [ε_x, κ_z, κ_y, γ_xy, γ_xz, φ_x]
    (κ_z before κ_y, matching implementation):
    ε_x = ∂u_x/∂x; κ_z = ∂θ_z/∂x; κ_y = ∂θ_y/∂x;
    γ_xy = ∂u_y/∂x - θ_z + α ∂²θ_z/∂x²;
    γ_xz = ∂u_z/∂x - θ_y + α ∂²θ_y/∂x²;
    φ_x = ∂θ_x/∂x. Coefficient α (for example h²/12 on rectangular sections) is a section property.

    Map: ``x(xi)`` linear on chord, ``dx/dxi = L/2``, ``dxi_dx = 2/L``, ``d2xi_dx2 = 4/L**2``.

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

    Notes
    -----
    Canonical `B` block (single Gauss point, Levinson row order):

    ```text
    ε = B U_e
    ε = [ε_x, κ_z, κ_y, γ_xy, γ_xz, φ_x]^T

    B row meanings:
    row 0: d(u_x)/dx
    row 1: d(θ_z)/dx
    row 2: d(θ_y)/dx
    row 3: d(u_y)/dx - θ_z + α d2(θ_z)/dx2
    row 4: d(u_z)/dx - θ_y + α d2(θ_y)/dx2
    row 5: d(θ_x)/dx
    ```

    **B tensor (per Gauss point, shape (6, 12))**
    - row 0 ``eps_x``: ``d(u_x)/dx``.
    - row 1 ``kappa_z``: ``d(theta_z)/dx``.
    - row 2 ``kappa_y``: ``d(theta_y)/dx``.
    - row 3 ``gamma_xy``: ``d(u_y)/dx - theta_z + alpha*d2(theta_z)/dx2``.
    - row 4 ``gamma_xz``: ``d(u_z)/dx - theta_y + alpha*d2(theta_y)/dx2``.
    - row 5 ``phi_x``: ``d(theta_x)/dx``.

    **D linkage and zeros**
    - ``S = D @ eps`` with ``S = [N, M_z, M_y, V_y, V_z, T]`` for this row order.
    - Shear rows are active; unlike Timoshenko constitutive form, Levinson shear stiffness is ``G*A``
      (no ``kappa`` factor), while higher-order terms enter via ``B`` through ``alpha``.

    **N tensor linkage**
    - Shear rows require both derivatives and shape values; inputs ``N``, ``dN_dxi``, ``d2N_dxi2``
      use batch shape ``(n_gp, 12, 6)``.
    - Entries not referenced by the row definitions remain zero.

    Same Gauss weak form as the shear-deformable beam family; see module one-liner for Voigt row order.

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