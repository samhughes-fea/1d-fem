# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/D_matrix.py
"""Material stiffness D (7, 7) for Vlasov warping EB: block-embeds linear EB (6, 6) and ``D[6,6] = E·Γ``.

Constitutive law ``S = D @ ε`` with ``ε``, ``S ∈ ℝ^7``. Rows 0–5 pair with the standard beam resultants
(``N``, ``M_y``, ``M_z``, ``V_y``, ``V_z``, ``T``); row 6 pairs with the warping strain and the corresponding
generalized stress for bimoment-type stiffness ``E·Γ`` in this implementation. Used in
``K_e += B.T @ D @ B * w_g * detJ``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator

from .constants import N_STRAIN


@dataclass(frozen=True)
class WarpingMaterialStiffnessOperator:
    """
    Constitutive tensor ``D ∈ ℝ^{7×7}`` for 3-D Euler–Bernoulli beam with Vlasov warping.

    **Strain and resultant vectors** (Voigt extension of the 6×6 table):

        ``ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x, ε_w]^T``   (7 strains; ``ε_w`` is the seventh row, warping strain)
        ``S = [N, M_y, M_z, V_y, V_z, T, S_w]^T``   (7 work-conjugate resultants; ``S_w`` pairs with ``ε_w``)

    In this implementation ``D[:6, :6]`` is identical to ``MaterialStiffnessOperator.assembly_form()`` for linear EB
    (shear rows 3–4 zero). The only coupling to the seventh DOF is ``D[6, 6] = E·Γ``; all other entries in row/column
    6 are zero.

    Notes
    -----
    **Sparsity structure of D (warping EB)**

    ```text
    D[:6, :6] = D_EB   (6, 6) from linear EB
    D[6, 6]   = E·Γ    (warping / bimoment stiffness; Γ [m⁶] from section)
    D[i, 6] = D[6, i] = 0  for i = 0..5
    ```

    **Weak-form linkage**

    ``B ∈ ℝ^{7×14}`` from ``WarpingStrainDisplacementOperator``; ``detJ = L/2``.

    Parameters
    ----------
    base_material_operator : MaterialStiffnessOperator
        Linear EB ``D`` (6, 6) for rows/columns 0–5.
    youngs_modulus : float
        Young's modulus E [Pa].
    warping_gamma : float
        Warping constant Γ [m⁶]; ``D[6, 6] = E * warping_gamma``.

    See Also
    --------
    MaterialStiffnessOperator
        Linear EB ``D`` in ``euler_bernoulli/utilities/D_matrix.py``.
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
        Baseline ``ε`` / ``S`` rows 0--5.
    """

    base_material_operator: MaterialStiffnessOperator
    youngs_modulus: float
    warping_gamma: float

    def assembly_form(self) -> np.ndarray:
        """
        Return the symmetric material matrix for stiffness assembly.

        Returns
        -------
        np.ndarray
            Shape (7, 7), symmetric.
        """
        D = np.zeros((N_STRAIN, N_STRAIN), dtype=np.float64)
        D[:6, :6] = self.base_material_operator.assembly_form()
        D[6, 6] = self.youngs_modulus * self.warping_gamma
        return D
