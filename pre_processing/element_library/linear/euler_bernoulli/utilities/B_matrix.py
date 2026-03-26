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
    Builds the strain-displacement matrix ``B`` for a 2-node 3-D Euler-Bernoulli beam.

    The operator maps shape-function derivatives to beam strain rows in Voigt order:
    ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x].

    **B tensor (per Gauss point, shape (6, 12))**
    - row 0 ε_x: axial terms from ∂u_x/∂x at DOFs 0 and 6.
    - row 1 κ_y: bending terms from ∂²u_z/∂x² and ∂²θ_y/∂x².
    - row 2 κ_z: bending terms from ∂²u_y/∂x² and ∂²θ_z/∂x².
    - row 3 γ_xy: identically zero in EB kinematics.
    - row 4 γ_xz: identically zero in EB kinematics.
    - row 5 φ_x: torsion terms from ∂θ_x/∂x at DOFs 3 and 9.

    **D linkage and zeros**
    - Parent constitutive step is ``S = D @ eps`` with ``S = [N, M_y, M_z, V_y, V_z, T]``.
    - Because EB shear strain rows are zero and EB ``D`` shear rows are zero, constitutive
      shear resultants ``V_y`` and ``V_z`` are zero in this operator path.
    - If shear force is needed for reporting, use equilibrium relation ``V = dM/dx``.

    **N tensor linkage**
    - Shape functions come from ``shape_functions.natural_coordinate_form`` as
      ``N``, ``dN_dxi``, ``d2N_dxi2`` with batch shape ``(n_gp, 12, 6)``.
    - ``B`` uses the derivative tensors; entries not referenced by the row rules above remain zero.

    Coordinate mapping: ``x(xi)`` linear on chord, ``dx/dxi = L/2``,
    ``dxi_dx = 2/L``, ``d2xi_dx2 = 4/L**2``.

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
    Canonical `B` block (single Gauss point, representative sparse pattern):

    ```text
    ε = B U_e
    B =
    [ b1,1  0     0     0     0     0    b1,7  0     0     0      0      0   ]  # ε_x
    [ 0     0    b2,3   0    b2,5   0     0    0    b2,9   0    b2,11    0   ]  # κ_y
    [ 0    b3,2   0     0     0    b3,6   0   b3,8   0     0      0    b3,12 ]  # κ_z
    [ 0     0     0     0     0     0     0    0     0     0      0      0   ]  # γ_xy = 0
    [ 0     0     0     0     0     0     0    0     0     0      0      0   ]  # γ_xz = 0
    [ 0     0     0    b6,4   0     0     0    0     0    b6,10   0      0   ]  # φ_x
    ```

    **Contract:** same outer sizes as standard 12-DOF beam: ``B`` (6, 12), ``U_e`` (12,).
    **Diff vs shear-deformable theories:** shear strain rows of ``eps`` stay zero here; Timoshenko/Levinson
    populate those rows with non-zero kinematics and ``D`` adds ``G*A`` (or ``kappa*G*A``) stiffness.

    Weak-form linkage: the element sums ``K_e += B.T @ D @ B * w_g * detJ`` over Gauss points using
    ``physical_coordinate_form`` for ``B``. Natural-coordinate ``B_tilde`` uses the Jacobian chain
    (``dxi_dx``, ``d2xi_dx2``).

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