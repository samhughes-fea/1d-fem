# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/B_matrix.py
"""Strain–displacement B for Vlasov warping EB: batch (n_gp, 7, 14) in physical x.

Composes the linear EB ``B`` (n_gp, 6, 12) on ``U_e`` indices 0–11 with a seventh strain row coupling ``θ_x``
and warping ``χ`` DOFs. At each Gauss point, ``ε = B @ U_e`` with ``ε ∈ ℝ^7``, ``U_e ∈ ℝ^{14}``.
Rows 0–5 follow Voigt order in ``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md``; row 6 is the warping
extension (bimoment-type strain rate). Parent element uses ``K_e += B.T @ D @ B * w_g * detJ`` with ``detJ = L/2``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator

from .constants import N_DOF, N_STRAIN


@dataclass(frozen=True)
class WarpingStrainDisplacementOperator:
    """
    Strain–displacement operator assembling ``B ∈ ℝ^{7×14}`` (per Gauss point) for warping EB.

    **Baseline embedding:** ``B[:6, :12]`` equals the linear EB ``StrainDisplacementOperator`` in physical
    coordinates. Row 6 has nonzeros only on ``θ_x`` DOFs (indices 3, 9 within the first 12 columns) and on
    ``χ`` DOFs (indices 12, 13), with coefficients ``±1/L`` from linear interpolation along the element.

    **Strain vector** ``ε ∈ ℝ^7`` (Voigt order):

        Row 0: ``ε_x``; 1: ``κ_y``; 2: ``κ_z``; 3: ``γ_xy``; 4: ``γ_xz``; 5: ``φ_x`` (St. Venant twist rate
        ``∂θ_x/∂x`` from the standard six rows); row 6: warping strain ``φ_x′ = ∂θ_x/∂x + ∂χ/∂x`` (non-uniform
        torsion / bimoment-type strain rate). Rows 3–4 remain zero for EB.

    **Displacement vector** ``U_e ∈ ℝ^{14}`` (node-major):

        ``[u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, χ¹, u_x², u_y², u_z², θ_x², θ_y², θ_z², χ²]``.

    Row 6 of ``B`` is **constant** in ξ (linear 2-node interpolation of ``θ_x`` and ``χ``): each Gauss point
    uses the same coefficients ``±1/L`` on DOFs 3, 9, 12, 13.

    Parameters
    ----------
    element_length : float
        Chord length L [m].
    base_strain_operator : StrainDisplacementOperator
        Linear EB operator producing ``B`` (n_gp, 6, 12) from ``(dN/dξ, d²N/dξ²)`` in physical x.

    See Also
    --------
    StrainDisplacementOperator
        Linear EB (6, 12) operator in ``euler_bernoulli/utilities/B_matrix.py``.
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
        Voigt order for rows 0--5; extensions for ``(14,) U_e`` and ``(7, 14) B``.
    """

    element_length: float
    base_strain_operator: StrainDisplacementOperator

    def warping_row(self) -> np.ndarray:
        """
        Row 6 of ``B``: ``φ_x′ = ∂θ_x/∂x + ∂χ/∂x`` with linear 2-node interpolation (``±1/L`` on ``θ_x`` and ``χ`` DOFs).

        Returns
        -------
        np.ndarray
            Shape (14,).
        """
        L = self.element_length
        row = np.zeros(N_DOF, dtype=np.float64)
        row[3] = -1.0 / L
        row[9] = 1.0 / L
        row[12] = -1.0 / L
        row[13] = 1.0 / L
        return row

    def physical_coordinate_form(
        self,
        dN_dξ: np.ndarray,
        d2N_dξ2: np.ndarray,
    ) -> np.ndarray:
        """
        Build ``B`` in physical coordinates for all Gauss points.

        Parameters
        ----------
        dN_dξ : np.ndarray
            First derivatives ``∂N/∂ξ``, shape (n_gauss, 12, 6).
        d2N_dξ2 : np.ndarray
            Second derivatives ``∂²N/∂ξ²``, shape (n_gauss, 12, 6).

        Returns
        -------
        np.ndarray
            Strain–displacement batch, shape (n_gauss, 7, 14).
        """
        n_gauss = dN_dξ.shape[0]
        B_6x12 = self.base_strain_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)
        B = np.zeros((n_gauss, N_STRAIN, N_DOF), dtype=np.float64)
        B[:, :6, :12] = B_6x12
        wr = self.warping_row()
        for g in range(n_gauss):
            B[g, 6, :] = wr
        return B
