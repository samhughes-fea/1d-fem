# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/D_matrix.py
"""Material stiffness D (7, 7): linear EB upper block and E·Γ on warping strain."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator

from .constants import N_STRAIN


@dataclass(frozen=True)
class WarpingMaterialStiffnessOperator:
    """
    Material stiffness D ∈ ℝ^{7×7}: upper 6×6 from linear EB; ``D[6, 6] = E * Gamma``.

    Parameters
    ----------
    base_material_operator : MaterialStiffnessOperator
        Linear EB (6, 6) assembly.
    youngs_modulus : float
        E (Pa).
    warping_gamma : float
        Γ (m⁶) warping constant; bimoment stiffness is E·Γ.
    """

    base_material_operator: MaterialStiffnessOperator
    youngs_modulus: float
    warping_gamma: float

    def assembly_form(self) -> np.ndarray:
        """Return D (7, 7)."""
        D = np.zeros((N_STRAIN, N_STRAIN), dtype=np.float64)
        D[:6, :6] = self.base_material_operator.assembly_form()
        D[6, 6] = self.youngs_modulus * self.warping_gamma
        return D
