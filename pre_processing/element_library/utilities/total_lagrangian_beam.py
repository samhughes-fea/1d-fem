# pre_processing/element_library/utilities/total_lagrangian_beam.py
"""
Total Lagrangian formulation for 2-node 3D beams (Green–Lagrange strain, K_σ).
Operator classes: GreenLagrangeStrainOperator, StressResultantOperator, GeometricStiffnessOperator.
Reference: initial configuration; strain E = E_lin + E_nl; tangent K_T = K_0 + K_σ.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)


@dataclass(frozen=True)
class GreenLagrangeStrainOperator:
    """
    Green–Lagrange strain measures for a 2-node 3D beam in Total Lagrangian formulation.

    Reference configuration: initial (undeformed) geometry; all quantities referred to it.
    Strain vector E = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]ᵀ with nonlinear terms:

        ε_x  = ∂u_x/∂x + ½(∂u_x/∂x)² + ½(∂u_y/∂x)² + ½(∂u_z/∂x)²   (axial, Green–Lagrange)
        κ_y  = ∂²u_z/∂x² + nonlinear correction (bending about y)
        κ_z  = ∂²u_y/∂x² + nonlinear correction (bending about z)
        γ_xy = 0 (Euler–Bernoulli) or ∂u_y/∂x − θ_z + ... (Timoshenko)
        γ_xz = 0 (EB) or ∂u_z/∂x − θ_y + ... (Timoshenko)
        φ_x  = ∂θ_x/∂x (torsion, linear)

    Coordinate mapping (reference): x(ξ) = ((1−ξ)/2)x₁ + ((1+ξ)/2)x₂, dx/dξ = L/2, ∂ξ/∂x = 2/L.

    Parameters
    ----------
    element_length : float
        Reference length L of the beam element (must be > 0).
    include_shear : bool
        If True, include shear strain (Timoshenko); if False, shear = 0 (Euler–Bernoulli).
    """

    element_length: float
    include_shear: bool = False

    def __post_init__(self) -> None:
        if self.element_length <= 0:
            raise ValueError(f"element_length must be positive, got {self.element_length}")

    def strain_linear_part(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        u_e: np.ndarray,
    ) -> np.ndarray:
        """
        Linear part of Green–Lagrange strain at one point: E_lin = B_lin @ u_e.

        Parameters
        ----------
        dN_dx : np.ndarray, shape (12, 6)
            First derivatives of shape functions w.r.t. x (physical).
        d2N_dx2 : np.ndarray, shape (12, 6)
            Second derivatives w.r.t. x (for curvature).
        u_e : np.ndarray, shape (12,)
            Element displacement vector [u_x1, u_y1, u_z1, θ_x1, θ_y1, θ_z1, u_x2, ...].

        Returns
        -------
        E_lin : np.ndarray, shape (6,)
            Linear strain vector [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x].
        """
        B_lin = self.linearized_strain_displacement(dN_dx, d2N_dx2)
        return (B_lin @ u_e).ravel()

    def strain_nonlinear_part(
        self,
        dN_dx: np.ndarray,
        u_e: np.ndarray,
    ) -> np.ndarray:
        """
        Nonlinear part of Green–Lagrange strain (e.g. axial: ½(∂u/∂x)² terms).

        Parameters
        ----------
        dN_dx : np.ndarray, shape (12, 6)
            First derivatives of shape functions w.r.t. x.
        u_e : np.ndarray, shape (12,)
            Element displacement vector.

        Returns
        -------
        E_nl : np.ndarray, shape (6,)
            Nonlinear strain contribution (axial row filled; bending/shear as needed).
        """
        # Axial: ½ ( (du/dx)² + (dv/dx)² + (dw/dx)² ). DOF 0,6 -> u_x; 1,7 -> u_y; 2,8 -> u_z.
        du_dx = dN_dx[0, 0] * u_e[0] + dN_dx[6, 0] * u_e[6]
        dv_dx = dN_dx[1, 1] * u_e[1] + dN_dx[7, 1] * u_e[7]
        dw_dx = dN_dx[2, 2] * u_e[2] + dN_dx[8, 2] * u_e[8]
        e_nl = np.zeros(6, dtype=np.float64)
        e_nl[0] = 0.5 * (du_dx ** 2 + dv_dx ** 2 + dw_dx ** 2)
        return e_nl

    def linearized_strain_displacement(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        N: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Linearized strain-displacement ``B_lin`` (6, 12); element sums ``K_0 += sum_g B_lin.T @ D @ B_lin * w_g * detJ``.

        Same structure as the linear beam B matrix. For Timoshenko shear, pass N (12×6) at the Gauss point.

        Returns
        -------
        B_lin : np.ndarray, shape (6, 12)
        """
        B_lin = np.zeros((6, 12), dtype=np.float64)
        # Axial: ε_x = ∂u_x/∂x  -> rows 0,6
        B_lin[0, 0] = dN_dx[0, 0]
        B_lin[0, 6] = dN_dx[6, 0]
        # Bending κ_y = ∂²u_z/∂x² (and θ_y terms if any)
        B_lin[1, 2] = d2N_dx2[2, 2]
        B_lin[1, 8] = d2N_dx2[8, 2]
        B_lin[1, 4] = d2N_dx2[4, 4]
        B_lin[1, 10] = d2N_dx2[10, 4]
        # Bending κ_z = ∂²u_y/∂x²
        B_lin[2, 1] = d2N_dx2[1, 1]
        B_lin[2, 7] = d2N_dx2[7, 1]
        B_lin[2, 5] = d2N_dx2[5, 5]
        B_lin[2, 11] = d2N_dx2[11, 5]
        # Shear (zero for EB; caller can override for Timoshenko)
        # B_lin[3], B_lin[4]
        # Shear (Timoshenko): γ_xy = ∂u_y/∂x − θ_z, γ_xz = ∂u_z/∂x − θ_y
        if self.include_shear:
            B_lin[3, 1] = dN_dx[1, 1]
            B_lin[3, 7] = dN_dx[7, 1]
            B_lin[3, 5] = -N[5, 5] if N is not None else -1.0
            B_lin[3, 11] = -N[11, 5] if N is not None else -1.0
            B_lin[4, 2] = dN_dx[2, 2]
            B_lin[4, 8] = dN_dx[8, 2]
            B_lin[4, 4] = -N[4, 4] if N is not None else -1.0
            B_lin[4, 10] = -N[10, 4] if N is not None else -1.0
        # Torsion φ_x = ∂θ_x/∂x
        B_lin[5, 3] = dN_dx[3, 3]
        B_lin[5, 9] = dN_dx[9, 3]
        return B_lin


@dataclass(frozen=True)
class StressResultantOperator:
    """
    Section forces (stress resultants) from 2nd Piola–Kirchhoff stress in reference configuration.

    S = D @ E  at each Gauss point; integrate or project to get:
        N   : axial force
        M_y : bending moment about y
        M_z : bending moment about z
        (and optionally torsion T, shear V_y, V_z for Timoshenko).

    Parameters
    ----------
    None (stateless; D and E are passed to methods).
    """

    def section_forces_from_strain(
        self,
        E: np.ndarray,
        D: np.ndarray,
    ) -> Tuple[float, float, float]:
        """
        Compute section force resultants from strain and material stiffness at one point.

        For a single Gauss point or average: N = (D @ E)[0], M_y = (D @ E)[1], M_z = (D @ E)[2].
        Full integration over the element would sum over Gauss points with weights.

        Parameters
        ----------
        E : np.ndarray, shape (6,)
            Strain vector [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x].
        D : np.ndarray, shape (6, 6)
            Material stiffness (same as linear D matrix).

        Returns
        -------
        N : float
            Axial force.
        M_y : float
            Bending moment about y.
        M_z : float
            Bending moment about z.
        """
        S = D @ E
        N = float(S[0])
        M_y = float(S[1])
        M_z = float(S[2])
        return N, M_y, M_z
