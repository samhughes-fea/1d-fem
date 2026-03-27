# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/B_matrix.py
"""Strain–displacement B (7, 14) per Gauss point: linear EB block plus warping row."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator

from .constants import N_DOF, N_STRAIN


@dataclass(frozen=True)
class WarpingStrainDisplacementOperator:
    """
    Strain–displacement B ∈ ℝ^{7×14} for warping EB: rows 0–5 from linear EB on DOFs 0–11;
    row 6 is φ_x′ = ∂θ_x/∂x + ∂χ/∂x (bimoment strain).

    Parameters
    ----------
    element_length : float
        Chord length L.
    base_strain_operator : StrainDisplacementOperator
        Linear EB (6, 12) operator on standard DOFs.
    """

    element_length: float
    base_strain_operator: StrainDisplacementOperator

    def warping_row(self) -> np.ndarray:
        """Row 6 of B: φ_x′ = dθ_x/dx + dχ/dx. Linear: ±1/L for θ_x and warping DOFs."""
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
        Build B (n_gauss, 7, 14) from EB B (n_gauss, 6, 12) plus constant warping row per Gauss point.

        Parameters
        ----------
        dN_dξ : np.ndarray
            Shape (n_gauss, 12, 6).
        d2N_dξ2 : np.ndarray
            Shape (n_gauss, 12, 6).

        Returns
        -------
        np.ndarray
            Shape (n_gauss, 7, 14).
        """
        n_gauss = dN_dξ.shape[0]
        B_6x12 = self.base_strain_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)
        B = np.zeros((n_gauss, N_STRAIN, N_DOF), dtype=np.float64)
        B[:, :6, :12] = B_6x12
        wr = self.warping_row()
        for g in range(n_gauss):
            B[g, 6, :] = wr
        return B
