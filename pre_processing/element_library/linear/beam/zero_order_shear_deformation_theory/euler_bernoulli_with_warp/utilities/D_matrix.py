# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/D_matrix.py
"""Material stiffness D (7, 7) for Euler–Bernoulli beam with Vlasov warping.

``S = D ε`` embeds linear EB ``D`` on rows/columns 0–5, zero shear rows, St. Venant torsion on ``D[5,5]``,
and warping stiffness ``E·Γ`` on ``D[6,6]``. Used in ``K_e += B.T @ D @ B * w_g * detJ`` in the parent element
(see ``linear_warping_euler_bernoulli_3D.py``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator

from .constants import N_STRAIN


@dataclass(frozen=True)
class WarpingMaterialStiffnessOperator:
    """
    Constitutive tensor ``D ∈ ℝ^{7×7}`` for 3-D Euler–Bernoulli beam elements with Vlasov warping.

    ``D`` is a rank-2 symmetric tensor relating the generalised strain vector ``ε ∈ ℝ^7`` to the
    section resultant vector ``S ∈ ℝ^7`` via ``S = D ε``:

        ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x, ε_w]^T   (Voigt strains; ``ε_w`` is the seventh row)
        S = [N,   M_y, M_z, V_y, V_z, T,   S_w]^T   (work-conjugate resultants; ``S_w`` pairs with ``ε_w``)

    ``D[:6, :6]`` is identical to ``MaterialStiffnessOperator.assembly_form()`` for linear EB (shear rows
    3–4 zero). The only entry involving the warping strain is ``D[6, 6] = E·Γ``; all other entries in
    row or column 6 are zero. Shear resultants ``V_y``, ``V_z`` remain zero from the constitutive path;
    shear forces are recovered from equilibrium where applicable, as for pure EB.

    Parameters
    ----------
    base_material_operator : MaterialStiffnessOperator
        Linear EB ``D`` (6, 6) for rows/columns 0–5 (``EA``, ``EI_y``, ``EI_z``, zero shear, ``G·J_t``).
    youngs_modulus : float
        Young's modulus ``E`` [Pa]; scales ``D[6, 6]`` together with ``warping_gamma``.
    warping_gamma : float
        Warping constant ``Γ`` [m⁶] from the section; ``D[6, 6] = E * warping_gamma``.

    Notes
    -----
    **Sparsity structure of D (warping EB, Voigt order)**

    ```text
    S = D ε
    ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x, ε_w]^T
    S = [N,   M_y, M_z, V_y, V_z, T,   S_w]^T

    D =
    [ EA     0      0      0    0    0     0    ]
    [ 0     EI_y    0      0    0    0     0    ]
    [ 0      0     EI_z    0    0    0     0    ]
    [ 0      0      0      0    0    0     0    ]
    [ 0      0      0      0    0    0     0    ]
    [ 0      0      0      0    0   G·J_t  0    ]
    [ 0      0      0      0    0    0    E·Γ   ]
    ```

    **Component definitions — D ∈ ℝ^{7×7}, block-diagonal extension of EB**

    ```text
    D[0,0] = E·A        (axial stiffness)
    D[1,1] = E·I_y      (bending stiffness about y)
    D[2,2] = E·I_z      (bending stiffness about z)
    D[3,3] = 0          (no constitutive shear; γ_xy = 0 by EB kinematic constraint)
    D[4,4] = 0          (no constitutive shear; γ_xz = 0 by EB kinematic constraint)
    D[5,5] = G·J_t      (St. Venant torsional stiffness)
    D[6,6] = E·Γ        (warping / bimoment stiffness; Γ from section)
    D[i,j] = 0  for all i ≠ j   (uncoupled across the seven modes in this implementation)
    ```

    **Weak-form assembly linkage**

    The element stiffness is accumulated as ``K_e += B.T @ D @ B * w_g * detJ`` with ``ξ ∈ [−1, 1]`` and
    ``detJ = L/2``. ``B ∈ ℝ^{7×14}`` comes from ``WarpingStrainDisplacementOperator``; the shape-function
    tensors ``N``, ``∂N/∂ξ``, ``∂²N/∂ξ²`` of batch shape ``(n_gp, 12, 6)`` come from the registry as for
    linear EB.

    See Also
    --------
    MaterialStiffnessOperator
        Linear EB ``D`` (6, 6) in ``euler_bernoulli/utilities/D_matrix.py``.
    WarpingStrainDisplacementOperator
        ``B`` (7, 14) in ``euler_bernoulli_with_warp/utilities/B_matrix.py``.
    linear_warping_euler_bernoulli_3D.LinearWarpingEulerBernoulliBeamElement3D
        Parent element.
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
        Baseline ``ε`` / ``S`` rows 0–5; extensions for seventh row/column.
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
            Shape ``(7, 7)``, symmetric.
        """
        D = np.zeros((N_STRAIN, N_STRAIN), dtype=np.float64)
        D[:6, :6] = self.base_material_operator.assembly_form()
        D[6, 6] = self.youngs_modulus * self.warping_gamma
        return D
