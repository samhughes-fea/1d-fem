# pre_processing/element_library/linear/timoshenko/utilities/B_matrix.py
"""Strain-displacement ``B`` (6, 12) per Gauss point for 2-node 3-D Timoshenko beam.

``eps = B @ U_e``; Voigt ``eps`` = [eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x] with shear
``gamma_xy = d(u_y)/dx - theta_z``, ``gamma_xz = d(u_z)/dx - theta_y``. Parent sums ``K_e += B.T @ D @ B * w_g * detJ``.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

@dataclass(frozen=True)
class StrainDisplacementOperator:
    """
    Strain-displacement ``B`` (6, 12) per Gauss point for a 2-node 3-D Timoshenko beam.

    Voigt ``eps = B @ U_e`` with rows
    ``[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]``:
    ``eps_x = d(u_x)/dx``; ``kappa_y = d(theta_y)/dx``, ``kappa_z = d(theta_z)/dx`` (rotation-based bending);
    ``gamma_xy = d(u_y)/dx - theta_z``, ``gamma_xz = d(u_z)/dx - theta_y``;
    ``phi_x = d(theta_x)/dx``.

    Map: ``x(xi)`` linear on chord, ``dx/dxi = L/2``, ``dxi_dx = 2/L``, ``d2xi_dx2 = 4/L**2``.

    Parameters
    ----------
    element_length : float
        Length `L` of the beam element (must be > 0)

    Attributes
    ----------
    jacobian : float
        Jacobian of coordinate mapping (L/2)
    dξ_dx : float
        First derivative ∂ξ/∂x (2/L)
    d2ξ_dx2 : float
        Second derivative ∂²ξ/∂x² (4/L²)

    Notes
    -----
    **B tensor (per Gauss point, shape (6, 12))**
    - row 0 ``eps_x``: ``d(u_x)/dx``.
    - row 1 ``kappa_y``: ``d(theta_y)/dx``.
    - row 2 ``kappa_z``: ``d(theta_z)/dx``.
    - row 3 ``gamma_xy``: ``d(u_y)/dx - theta_z`` (non-zero).
    - row 4 ``gamma_xz``: ``d(u_z)/dx - theta_y`` (non-zero).
    - row 5 ``phi_x``: ``d(theta_x)/dx``.

    **D linkage and zeros**
    - ``S = D @ eps`` with ``S = [N, M_y, M_z, V_y, V_z, T]``.
    - Unlike EB, shear rows 3 and 4 are active in both ``B`` and ``D``.
    - Parent ``D`` uses ``kappa*G*A`` on shear rows.

    **N tensor linkage**
    - Shear rows use shape functions ``N`` (not only derivatives): input tensors
      ``N``, ``dN_dxi``, ``d2N_dxi2`` are batched ``(n_gp, 12, 6)``.
    - If ``N`` is omitted, shear rows cannot be assembled and remain zero in this utility call.

    Weak-form linkage: ``linear_timoshenko_3D`` uses ``physical_coordinate_form`` in the stiffness loop; ``detJ = L/2``.

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