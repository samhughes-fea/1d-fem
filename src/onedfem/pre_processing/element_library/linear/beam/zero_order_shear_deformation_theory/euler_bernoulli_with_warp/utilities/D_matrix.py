# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/D_matrix.py
"""Material stiffness D (7, 7) for Euler–Bernoulli beam with Vlasov warping.

``S = D ε`` embeds linear EB on rows/columns 0–5, St. Venant torsion ``G·J_t`` on ``D[5,5]`` (shear-dominated
uniform twist), and warping / bimoment-type stiffness ``E·Γ`` on ``D[6,6]`` (see class docstring for **G** vs **E**).
Used in ``K_e += B.T @ D @ B * w_g * detJ`` for unified linear EB with warping (see ``linear_euler_bernoulli_3D.py``).
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

    For **open** thin-walled sections, a **warping DOF** (or warping amplitude) together with a section
    **warping constant** ``Γ`` (m⁶; sectorial / warping inertia) is standard in formulations that extend
    **St. Venant** uniform torsion with **non-uniform** torsion and **bimoment**-type effects. In this 1D
    model, row 5 pairs twist rate ``φ_x`` with torque ``T`` via **``G·J_t``**: uniform twist is
    **shear-dominated**, so **G** is the appropriate modulus. Row 6 pairs the warping strain with ``S_w``
    via **``E·Γ``**: restrained warping drives **longitudinal normal** stresses in the beam direction
    (fibers stretch as the section warps out of plane and is constrained)—an **axial-stretch-type**
    mechanism in the warping sense—so the usual Vlasov-type 1D stiffness uses **Young’s modulus E**
    times ``Γ``, not **G·Γ**.

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
    D[5,5] = G·J_t      (St. Venant torsional stiffness; shear-dominated uniform twist)
    D[6,6] = E·Γ        (warping / bimoment stiffness; Γ [m⁶]; normal-stress-type warping mode)
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
    linear_euler_bernoulli_3D.LinearEulerBernoulliBeamElement3D
        Unified linear EB; use with `[warping]` / warping mesh for 14 local DOFs.
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
