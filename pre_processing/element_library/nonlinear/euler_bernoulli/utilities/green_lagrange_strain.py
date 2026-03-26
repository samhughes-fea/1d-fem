# pre_processing/element_library/nonlinear/euler_bernoulli/utilities/green_lagrange_strain.py
"""
Green-Lagrange strain utilities for 2-node 3D Euler-Bernoulli (Total Lagrangian).

Per Gauss point: ``E`` (6,), ``B_lin`` and ``B_nl`` each (6, 12); inputs ``dN_dx``, ``d2N_dx2`` (12, 6), ``u_e`` (12,).
Shear rows of ``E`` / ``B_lin`` / ``B_nl`` are zero. Parent element sums material tangent using ``B_lin + B_nl`` and ``D``.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GreenLagrangeStrainOperator:
    """
    Green–Lagrange strain measures for a 2-node 3D Euler–Bernoulli beam (Total Lagrangian).

    Reference configuration: initial (undeformed) geometry; all quantities referred to it.
    Strain vector E = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]ᵀ with nonlinear terms:

        ε_x  = ∂u_x/∂x + ½(∂u_x/∂x)² + ½(∂u_y/∂x)² + ½(∂u_z/∂x)²   (axial, Green–Lagrange)
        κ_y  = ∂²u_z/∂x² + (∂u_x/∂x)(∂²u_z/∂x²) + O(ε²)   (bending about y; axial–curvature coupling)
        κ_z  = ∂²u_y/∂x² + (∂u_x/∂x)(∂²u_y/∂x²) + O(ε²)   (bending about z; axial–curvature coupling)
        γ_xy = 0, γ_xz = 0   (Euler–Bernoulli: no shear deformation)
        φ_x  = ∂θ_x/∂x (torsion, linear)

    Full nonlinear curvature is included: κ_y and κ_z use the axial–curvature coupling
    terms; the tangent stiffness uses B_lin + B_nl (see nonlinear_strain_displacement_gradient).

    Coordinate mapping (reference): x(ξ) = ((1−ξ)/2)x₁ + ((1+ξ)/2)x₂, dx/dξ = L/2, ∂ξ/∂x = 2/L.

    Parameters
    ----------
    element_length : float
        Reference length L of the beam element (must be > 0).
    include_shear : bool
        For Euler-Bernoulli this is False (no shear strain); kept for API consistency with Timoshenko operator.

    Notes
    -----
    This module only evaluates strain and ``B`` operators at a station; the element loops Gauss points with ``w_g`` and ``detJ``.
    Chord map: ``dx/dxi = L/2``. Moderate-rotation TL: axial Green-Lagrange term and axial-bending curvature coupling as coded.

    See Also
    --------
    nonlinear_euler_bernoulli_3D.NonlinearEulerBernoulliBeamElement3D
    docs/element_library/total_lagrangian_beam_formulation.md
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
        d2N_dx2: np.ndarray,
        u_e: np.ndarray,
    ) -> np.ndarray:
        """
        Nonlinear part of Green–Lagrange strain: axial ½(∂u/∂x)² and curvature (∂u_x/∂x)(∂²u/∂x²).

        Parameters
        ----------
        dN_dx : np.ndarray, shape (12, 6)
            First derivatives of shape functions w.r.t. x.
        d2N_dx2 : np.ndarray, shape (12, 6)
            Second derivatives w.r.t. x (for curvature).
        u_e : np.ndarray, shape (12,)
            Element displacement vector.

        Returns
        -------
        E_nl : np.ndarray, shape (6,)
            Nonlinear strain [ε_nl, κ_y_nl, κ_z_nl, 0, 0, 0].
        """
        du_dx = dN_dx[0, 0] * u_e[0] + dN_dx[6, 0] * u_e[6]
        dv_dx = dN_dx[1, 1] * u_e[1] + dN_dx[7, 1] * u_e[7]
        dw_dx = dN_dx[2, 2] * u_e[2] + dN_dx[8, 2] * u_e[8]
        # ∂²u_z/∂x² (κ_y): dofs 2, 8 (u_z), 4, 10 (θ_y)
        d2u_z_dx2 = (
            d2N_dx2[2, 2] * u_e[2] + d2N_dx2[8, 2] * u_e[8]
            + d2N_dx2[4, 4] * u_e[4] + d2N_dx2[10, 4] * u_e[10]
        )
        # ∂²u_y/∂x² (κ_z): dofs 1, 7 (u_y), 5, 11 (θ_z)
        d2u_y_dx2 = (
            d2N_dx2[1, 1] * u_e[1] + d2N_dx2[7, 1] * u_e[7]
            + d2N_dx2[5, 5] * u_e[5] + d2N_dx2[11, 5] * u_e[11]
        )
        e_nl = np.zeros(6, dtype=np.float64)
        e_nl[0] = 0.5 * (du_dx ** 2 + dv_dx ** 2 + dw_dx ** 2)
        e_nl[1] = du_dx * d2u_z_dx2   # κ_y nonlinear
        e_nl[2] = du_dx * d2u_y_dx2   # κ_z nonlinear
        return e_nl

    def nonlinear_strain_displacement_gradient(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        u_e: np.ndarray,
    ) -> np.ndarray:
        """
        Gradient of E_nl w.r.t. u_e: B_nl such that dE_nl/du_e = B_nl (6×12).
        Used for consistent tangent: K_mat = ∫ (B_lin + B_nl)ᵀ D (B_lin + B_nl) dx.

        Parameters
        ----------
        dN_dx : np.ndarray, shape (12, 6)
            First derivatives of shape functions w.r.t. x.
        d2N_dx2 : np.ndarray, shape (12, 6)
            Second derivatives w.r.t. x.
        u_e : np.ndarray, shape (12,)
            Element displacement vector.

        Returns
        -------
        B_nl : np.ndarray, shape (6, 12)
        """
        du_dx = dN_dx[0, 0] * u_e[0] + dN_dx[6, 0] * u_e[6]
        dv_dx = dN_dx[1, 1] * u_e[1] + dN_dx[7, 1] * u_e[7]
        dw_dx = dN_dx[2, 2] * u_e[2] + dN_dx[8, 2] * u_e[8]
        d2u_z_dx2 = (
            d2N_dx2[2, 2] * u_e[2] + d2N_dx2[8, 2] * u_e[8]
            + d2N_dx2[4, 4] * u_e[4] + d2N_dx2[10, 4] * u_e[10]
        )
        d2u_y_dx2 = (
            d2N_dx2[1, 1] * u_e[1] + d2N_dx2[7, 1] * u_e[7]
            + d2N_dx2[5, 5] * u_e[5] + d2N_dx2[11, 5] * u_e[11]
        )
        B_nl = np.zeros((6, 12), dtype=np.float64)
        # Row 0: d(ε_nl)/du_e
        B_nl[0, 0] = du_dx * dN_dx[0, 0]
        B_nl[0, 6] = du_dx * dN_dx[6, 0]
        B_nl[0, 1] = dv_dx * dN_dx[1, 1]
        B_nl[0, 7] = dv_dx * dN_dx[7, 1]
        B_nl[0, 2] = dw_dx * dN_dx[2, 2]
        B_nl[0, 8] = dw_dx * dN_dx[8, 2]
        # Row 1: d(κ_y_nl)/du_e = (∂²u_z/∂x²) d(∂u_x/∂x)/du_e + (∂u_x/∂x) d(∂²u_z/∂x²)/du_e
        B_nl[1, 0] = d2u_z_dx2 * dN_dx[0, 0]
        B_nl[1, 6] = d2u_z_dx2 * dN_dx[6, 0]
        B_nl[1, 2] = du_dx * d2N_dx2[2, 2]
        B_nl[1, 8] = du_dx * d2N_dx2[8, 2]
        B_nl[1, 4] = du_dx * d2N_dx2[4, 4]
        B_nl[1, 10] = du_dx * d2N_dx2[10, 4]
        # Row 2: d(κ_z_nl)/du_e
        B_nl[2, 0] = d2u_y_dx2 * dN_dx[0, 0]
        B_nl[2, 6] = d2u_y_dx2 * dN_dx[6, 0]
        B_nl[2, 1] = du_dx * d2N_dx2[1, 1]
        B_nl[2, 7] = du_dx * d2N_dx2[7, 1]
        B_nl[2, 5] = du_dx * d2N_dx2[5, 5]
        B_nl[2, 11] = du_dx * d2N_dx2[11, 5]
        return B_nl

    def linearized_strain_displacement(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        N: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Linearized strain-displacement ``B_lin`` (6, 12); element uses ``B_lin`` in ``F_int`` and in linear stiffness part.

        Same structure as the linear Euler–Bernoulli B matrix. Shear rows (3, 4) are zero.

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
        # Shear: γ_xy = γ_xz = 0 for Euler–Bernoulli (rows 3, 4 left zero)
        # Torsion φ_x = ∂θ_x/∂x
        B_lin[5, 3] = dN_dx[3, 3]
        B_lin[5, 9] = dN_dx[9, 3]
        return B_lin
