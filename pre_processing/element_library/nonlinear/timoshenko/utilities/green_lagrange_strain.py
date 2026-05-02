# pre_processing/element_library/nonlinear/timoshenko/utilities/green_lagrange_strain.py
"""
Green-Lagrange strain for 2-node 3D Timoshenko (Total Lagrangian).

``E`` (6,), ``B_lin`` (6, 12) matches linear Timoshenko ``B`` at each station when ``N`` is supplied for shear.
With ``include_shear=True``, nonlinear centroid shear rows ``gamma_xy``, ``gamma_xz`` use ``N`` at the Gauss point.
Inputs ``dN_dx``, ``d2N_dx2`` (12, 6), ``u_e`` (12,).

Parent element uses ``B`` from this operator with ``S = D @ E`` and Gauss sums for ``F_int`` and tangents.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GreenLagrangeStrainOperator:
    """
    Green–Lagrange strain measures for a 2-node 3D beam in Total Lagrangian formulation.

    Reference configuration: initial (undeformed) geometry; all quantities referred to it.
    Strain vector E = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]ᵀ with nonlinear terms:

        ε_x  = ∂u_x/∂x + ½(∂u_x/∂x)² + ½(∂u_y/∂x)² + ½(∂u_z/∂x)²   (axial, Green–Lagrange)
        κ_y  = ∂θ_y/∂x + (∂u_x/∂x)(∂²θ_y/∂x²) + …   (linear Timoshenko row plus axial–curvature coupling on θ_y)
        κ_z  = ∂θ_z/∂x + (∂u_x/∂x)(∂²θ_z/∂x²) + …
        γ_xy = ∂u_y/∂x − θ_z + γ_nl_xy   (Timoshenko; γ_nl_xy = −θ_z·u_x' + θ_x·u_z' at centroid when shear NL on)
        γ_xz = ∂u_z/∂x − θ_y + γ_nl_xz   (γ_nl_xz = θ_y·u_x' − θ_x·u_y')
        φ_x  = ∂θ_x/∂x (torsion, linear)

    When ``include_shear`` is False, γ_xy and γ_xz rows have no nonlinear supplement and shear rows of ``B_lin`` are zero.

    Coordinate mapping (reference): x(ξ) = ((1−ξ)/2)x₁ + ((1+ξ)/2)x₂, dx/dξ = L/2, ∂ξ/∂x = 2/L.

    Parameters
    ----------
    element_length : float
        Reference length L of the beam element (must be > 0).
    include_shear : bool
        If True, Timoshenko shear strain rows and centroid nonlinear shear; if False, shear rows zero (EB-like).

    Notes
    -----
    Per-Gauss-point only; element supplies ``w_g``, ``detJ``. Nonlinear ``B_nl`` is defined in
    ``nonlinear_strain_displacement_gradient``.     Parent TL Timoshenko uses ``B_tot = B_lin + B_nl`` for ``F_int`` and ``K_delta`` vs ``K_0`` (see ``nonlinear_timoshenko_3D``).

    See Also
    --------
    nonlinear_timoshenko_3D.NonlinearTimoshenkoBeamElement3D
    nonlinear.euler_bernoulli.utilities.green_lagrange_strain (EB variant)
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
            Second derivatives w.r.t. x (for curvature rows in this operator).
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
        N: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Nonlinear part of Green–Lagrange strain: axial quadratic terms, axial–curvature coupling on κ_y/κ_z,
        and centroid shear remainders when ``include_shear`` and ``N`` are provided.

        Parameters
        ----------
        dN_dx : np.ndarray, shape (12, 6)
            First derivatives of shape functions w.r.t. x.
        d2N_dx2 : np.ndarray, shape (12, 6)
            Second derivatives w.r.t. x (for κ_y, κ_z nonlinear terms).
        u_e : np.ndarray, shape (12,)
            Element displacement vector.
        N : np.ndarray, optional, shape (12, 6)
            Shape functions at the Gauss point (required for nonlinear shear when ``include_shear`` is True).

        Returns
        -------
        E_nl : np.ndarray, shape (6,)
            Nonlinear strain contribution.
        """
        du_dx = dN_dx[0, 0] * u_e[0] + dN_dx[6, 0] * u_e[6]
        dv_dx = dN_dx[1, 1] * u_e[1] + dN_dx[7, 1] * u_e[7]
        dw_dx = dN_dx[2, 2] * u_e[2] + dN_dx[8, 2] * u_e[8]
        # ∂²θ_y/∂x², ∂²θ_z/∂x² — align κ_nl with linear κ rows (rotation-based bending)
        d2theta_y_dx2 = d2N_dx2[4, 4] * u_e[4] + d2N_dx2[10, 4] * u_e[10]
        d2theta_z_dx2 = d2N_dx2[5, 5] * u_e[5] + d2N_dx2[11, 5] * u_e[11]
        e_nl = np.zeros(6, dtype=np.float64)
        e_nl[0] = 0.5 * (du_dx**2 + dv_dx**2 + dw_dx**2)
        e_nl[1] = du_dx * d2theta_y_dx2
        e_nl[2] = du_dx * d2theta_z_dx2

        if self.include_shear and N is not None:
            theta_x = N[3, 3] * u_e[3] + N[9, 3] * u_e[9]
            theta_y = N[4, 4] * u_e[4] + N[10, 4] * u_e[10]
            theta_z = N[5, 5] * u_e[5] + N[11, 5] * u_e[11]
            # Centroid: 2E_xy − γ_lin_xy and 2E_xz − γ_lin_xz with γ_lin = v' − θ_z, w' − θ_y
            e_nl[3] = -theta_z * du_dx + theta_x * dw_dx
            e_nl[4] = theta_y * du_dx - theta_x * dv_dx

        return e_nl

    def nonlinear_strain_displacement_gradient(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        u_e: np.ndarray,
        N: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Gradient of E_nl w.r.t. u_e: B_nl such that dE_nl/du_e ≈ B_nl (6×12).

        Row 0 matches the EB axial NL gradient; rows 1–2 use rotation-based κ_nl (θ_y, θ_z); rows 3–4 are γ_nl when shear NL is active.

        Parameters
        ----------
        dN_dx : np.ndarray, shape (12, 6)
            First derivatives of shape functions w.r.t. x.
        d2N_dx2 : np.ndarray, shape (12, 6)
            Second derivatives w.r.t. x.
        u_e : np.ndarray, shape (12,)
            Element displacement vector.
        N : np.ndarray, optional, shape (12, 6)
            Shape functions at the Gauss point (required when ``include_shear`` and using shear rows).

        Returns
        -------
        B_nl : np.ndarray, shape (6, 12)
        """
        du_dx = dN_dx[0, 0] * u_e[0] + dN_dx[6, 0] * u_e[6]
        dv_dx = dN_dx[1, 1] * u_e[1] + dN_dx[7, 1] * u_e[7]
        dw_dx = dN_dx[2, 2] * u_e[2] + dN_dx[8, 2] * u_e[8]
        d2theta_y_dx2 = d2N_dx2[4, 4] * u_e[4] + d2N_dx2[10, 4] * u_e[10]
        d2theta_z_dx2 = d2N_dx2[5, 5] * u_e[5] + d2N_dx2[11, 5] * u_e[11]
        B_nl = np.zeros((6, 12), dtype=np.float64)
        # Row 0: d(ε_nl)/du_e
        B_nl[0, 0] = du_dx * dN_dx[0, 0]
        B_nl[0, 6] = du_dx * dN_dx[6, 0]
        B_nl[0, 1] = dv_dx * dN_dx[1, 1]
        B_nl[0, 7] = dv_dx * dN_dx[7, 1]
        B_nl[0, 2] = dw_dx * dN_dx[2, 2]
        B_nl[0, 8] = dw_dx * dN_dx[8, 2]
        # Row 1: d(κ_y_nl)/du_e with κ_y_nl = u_x' · ∂²θ_y/∂x²
        B_nl[1, 0] = d2theta_y_dx2 * dN_dx[0, 0]
        B_nl[1, 6] = d2theta_y_dx2 * dN_dx[6, 0]
        B_nl[1, 4] = du_dx * d2N_dx2[4, 4]
        B_nl[1, 10] = du_dx * d2N_dx2[10, 4]
        # Row 2: d(κ_z_nl)/du_e with κ_z_nl = u_x' · ∂²θ_z/∂x²
        B_nl[2, 0] = d2theta_z_dx2 * dN_dx[0, 0]
        B_nl[2, 6] = d2theta_z_dx2 * dN_dx[6, 0]
        B_nl[2, 5] = du_dx * d2N_dx2[5, 5]
        B_nl[2, 11] = du_dx * d2N_dx2[11, 5]

        if self.include_shear and N is not None:
            theta_x = N[3, 3] * u_e[3] + N[9, 3] * u_e[9]
            theta_y = N[4, 4] * u_e[4] + N[10, 4] * u_e[10]
            theta_z = N[5, 5] * u_e[5] + N[11, 5] * u_e[11]
            # γ_nl_xy = -theta_z * du_dx + theta_x * dw_dx
            # ∂/∂u_a: -∂theta_z/∂u_a * du_dx - theta_z * ∂(du_dx)/∂u_a + ∂theta_x/∂u_a * dw_dx + theta_x * ∂(dw_dx)/∂u_a
            B_nl[3, 0] = -theta_z * dN_dx[0, 0]
            B_nl[3, 6] = -theta_z * dN_dx[6, 0]
            B_nl[3, 2] = theta_x * dN_dx[2, 2]
            B_nl[3, 8] = theta_x * dN_dx[8, 2]
            B_nl[3, 3] = dw_dx * N[3, 3]
            B_nl[3, 9] = dw_dx * N[9, 3]
            B_nl[3, 5] = -du_dx * N[5, 5]
            B_nl[3, 11] = -du_dx * N[11, 5]
            # γ_nl_xz = theta_y * du_dx - theta_x * dv_dx
            B_nl[4, 0] = theta_y * dN_dx[0, 0]
            B_nl[4, 6] = theta_y * dN_dx[6, 0]
            B_nl[4, 1] = -theta_x * dN_dx[1, 1]
            B_nl[4, 7] = -theta_x * dN_dx[7, 1]
            B_nl[4, 3] = -dv_dx * N[3, 3]
            B_nl[4, 9] = -dv_dx * N[9, 3]
            B_nl[4, 4] = du_dx * N[4, 4]
            B_nl[4, 10] = du_dx * N[10, 4]

        return B_nl

    def linearized_strain_displacement(
        self,
        dN_dx: np.ndarray,
        d2N_dx2: np.ndarray,
        N: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Linearized strain-displacement ``B_lin`` (6, 12); parent caches ``K_0 += sum_g B_lin.T @ D @ B_lin * w_g * detJ``.

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
