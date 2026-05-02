# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/B_matrix.py
"""Strain–displacement B (7, 14) per Gauss point for 2-node 3-D Euler–Bernoulli beam with Vlasov warping.

At each Gauss point, ``ε = B U_e`` with ``ε ∈ ℝ^7`` and ``U_e ∈ ℝ^{14}``; rows 0–5 match
``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`` / linear EB; row 6 is the Vlasov-type warping strain
rate (``∂θ_x/∂x + ∂χ/∂x``), work-conjugate to ``S_w`` with ``E·Γ`` on ``D`` (see class docstring).
Parent element uses ``K_e += B.T @ D @ B * w_g * detJ`` with ``detJ = L/2``.
Voigt order follows ``FORMULATION_DOCSTRING_STANDARDS.md`` (extensions for ``(14,) U_e`` and ``(7, 14) B``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator

from .constants import N_DOF, N_STRAIN


@dataclass(frozen=True)
class WarpingStrainDisplacementOperator:
    """
    Strain–displacement tensor ``B ∈ ℝ^{7×14}`` for a 2-node 3-D Euler–Bernoulli beam with Vlasov warping.

    ``B`` is a rank-2 tensor defined at each Gauss point such that ``ε = B U_e``,
    where ``ε ∈ ℝ^7`` is the generalised strain vector and ``U_e ∈ ℝ^{14}`` is the
    element displacement vector. **Column layout (stiffness assembly):** indices 0–11 are the
    standard 12-DOF EB packing; indices 12–13 are nodal warping intensities ``χ¹``, ``χ²``:

        U_e = [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x², u_y², u_z², θ_x², θ_y², θ_z², χ¹, χ²]^T

    **Kinematic equations (Voigt rows 0–6)**

    Rows 0–5 are identical in meaning to linear EB (see ``StrainDisplacementOperator`` and
    ``FORMULATION_DOCSTRING_STANDARDS.md``); columns 12–13 do not enter these rows:

        ε_x   = ∂u_x/∂x                          (axial extension)
        κ_y   = ∂²u_z/∂x²                       (curvature about y)
        κ_z   = ∂²u_y/∂x²                       (curvature about z)
        γ_xy  = 0                                (EB shear-inextensibility)
        γ_xz  = 0                                (EB shear-inextensibility)
        φ_x   = ∂θ_x/∂x                          (St. Venant twist rate; row 5)

    Row 6 is the warping / non-uniform torsion strain rate (work-conjugate to ``S_w`` via ``D``;
    see ``WarpingMaterialStiffnessOperator``):

        φ_x′  = ∂θ_x/∂x + ∂χ/∂x                  (row 6; bimoment-type strain rate)

    **Theory context (Vlasov-type 1D beam)**

    For open sections, codes that add **non-uniform** torsion and **bimoment**-type physics to **St. Venant**
    torsion introduce a warping DOF (or amplitude) and section data **Γ** (m⁶). In that family, the kinematic
    extension is **twist rate** plus **warping gradient**; this implementation’s seventh strain row
    ``φ_x′ = ∂θ_x/∂x + ∂χ/∂x`` is paired with ``S_w`` via ``D`` with **``E·Γ``** on ``D[6,6]``, while row 5
    ``φ_x = ∂θ_x/∂x`` pairs with **``G·J_t``** for uniform (shear-dominated) torsion—see
    ``WarpingMaterialStiffnessOperator`` for **G** vs **E**.

    **Constitutive hook**

    ``S = D ε`` with ``S ∈ ℝ^7`` and ``D ∈ ℝ^{7×7}`` from
    ``euler_bernoulli_with_warp/utilities/D_matrix.py``. Rows 0–5 pair with the usual section
    resultants; row 6 pairs with the warping resultant ``S_w`` (notation as in that module).

    Parameters
    ----------
    element_length : float
        Chord length ``L`` of the beam element [m]; must match ``base_strain_operator.element_length``
        (same physical element).
    base_strain_operator : StrainDisplacementOperator
        Linear EB operator that maps ``(dN/dξ, d²N/dξ²)`` to ``B`` of shape ``(n_gp, 6, 12)`` in
        physical ``x``. Its ``physical_coordinate_form`` supplies rows 0–5 and columns 0–11;
        ``dξ_dx`` and ``d²ξ/dx²`` live on that object (``2/L``, ``4/L²``).

    Notes
    -----
    **Isoparametric mapping**

    The same map ``x(ξ) = L(1+ξ)/2`` applies as for linear EB; ``∂ξ/∂x = 2/L``, ``∂²ξ/∂x² = 4/L²``.
    Rows 0–5 are assembled by ``base_strain_operator.physical_coordinate_form``, which applies
    the chain rule to the registry Hermite/linear shape-function derivatives. Row 6 is assembled
    **from that same** ``B`` batch ``B_6x12`` (no separate hardcoded ``1/L``): columns ``0:12``
    copy **row 5** (``∂θ_x/∂x``); columns **12–13** take **row 0** entries ``B[0,0]`` and ``B[0,6]``
    (``∂L₁/∂x``, ``∂L₂/∂x`` for axial ``u_x``), which are the same linear Lagrange slopes as
    ``∂χ/∂x`` when χ uses the same ``L₁, L₂`` as ``u_x``.

    **Shape function basis**

    Rows 0–5 use the same ``N``, ``∂N/∂ξ``, ``∂²N/∂ξ²`` pipeline and column coupling as
    ``StrainDisplacementOperator`` (see ``euler_bernoulli/utilities/B_matrix.py``). Row 6 reuses
    rows 0 and 5 of ``B_6x12`` as above. For **consistent mass**, χ interpolation in the extended
    ``N`` tensor is implemented in ``shape_functions.extend_natural_shape_to_warping``.

    **Sparsity structure of B (single Gauss point)**

    ```text
    ε = B U_e,   ε ∈ ℝ^7,   U_e ∈ ℝ^{14}
    Cols 0–6:   [u_x¹, u_y¹, u_z¹, θ_x¹, θ_y¹, θ_z¹, u_x²]
    Cols 7–13:  [u_y², u_z², θ_x², θ_y², θ_z², χ¹, χ²]

    Rows 0–5 on cols 0–11: identical symbolic pattern to linear EB B (6×12); cols 12–13 are 0.

    B =
    [ b1,1  0     0     0     0     0    b1,7  0     0     0      0      0      0   0   ]  # ε_x
    [ 0     0    b2,3   0    b2,5   0     0    0    b2,9   0    b2,11    0      0   0   ]  # κ_y
    [ 0    b3,2   0     0     0    b3,6   0   b3,8   0     0      0    b3,12    0   0   ]  # κ_z
    [ 0     0     0     0     0     0     0    0     0     0      0      0      0   0   ]  # γ_xy = 0
    [ 0     0     0     0     0     0     0    0     0     0      0      0      0   0   ]  # γ_xz = 0
    [ 0     0     0    b6,4   0     0     0    0     0    b6,10   0      0      0   0   ]  # φ_x = ∂θ_x/∂x
    [ same as row 5 on cols 0–11;  b1,1  0  …  b1,7  on cols 12–13 ]  # row 6: copy row 5; χ cols = row-0 axial slopes
    ```

    **Non-zero entries of B in physical coordinates**

    Rows 0–5: identical formulas to ``euler_bernoulli/utilities/B_matrix.py`` (section *Non-zero
    entries of B in physical coordinates*); columns 12–13 are zero.

    Row 6 (warping strain rate ``φ_x′ = ∂θ_x/∂x + ∂χ/∂x``), from ``B_6x12``:

    ```text
      B[6, j]  = B[5, j]   for j = 0..11   (∂θ_x/∂x)
      B[6, 12] = B[0, 0]    (∂χ/∂x on χ¹; same ∂L₁/∂x as ε_x row)
      B[6, 13] = B[0, 6]    (∂χ/∂x on χ²; same ∂L₂/∂x as ε_x row)
    ```

    For 2-node linear ``L₁, L₂``, these match ``±1/L`` everywhere.

    **Batch shape**

    ``physical_coordinate_form`` returns ``(n_gp, 7, 14)``. Rows 0–5 vary with ``ξ_g``; row 6 is
    independent of ``ξ_g`` for linear ``L₁, L₂`` (same entries at every Gauss point).

    Weak-form assembly: ``K_e += B.T @ D @ B * w_g * detJ`` with ``detJ = L/2``.

    See Also
    --------
    StrainDisplacementOperator
        Linear EB ``B`` (6, 12) in ``euler_bernoulli/utilities/B_matrix.py``.
    WarpingMaterialStiffnessOperator
        ``D`` (7, 7) in ``euler_bernoulli_with_warp/utilities/D_matrix.py``.
    linear_euler_bernoulli_3D.LinearEulerBernoulliBeamElement3D
        Unified linear EB; use with `[warping]` for operators extended to 14 DOFs.
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
        Voigt order for rows 0–5; extensions for ``(14,) U_e`` and ``(7, 14) B``.
    """

    element_length: float
    base_strain_operator: StrainDisplacementOperator

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
            First derivatives ``∂N/∂ξ``, shape ``(n_gauss, 12, 6)``.
        d2N_dξ2 : np.ndarray
            Second derivatives ``∂²N/∂ξ²``, shape ``(n_gauss, 12, 6)``.

        Returns
        -------
        np.ndarray
            Strain–displacement batch, shape ``(n_gauss, 7, 14)``.
        """
        B_6x12 = self.base_strain_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)
        n_gauss = B_6x12.shape[0]
        B = np.zeros((n_gauss, N_STRAIN, N_DOF), dtype=np.float64)
        B[:, :6, :12] = B_6x12
        # Row 6: ∂θ_x/∂x + ∂χ/∂x from B_6x12 — row 5 + χ slopes = same L₁/L₂ as axial row 0, cols 0 and 6.
        B[:, 6, :12] = B_6x12[:, 5, :]
        B[:, 6, 12] = B_6x12[:, 0, 0]
        B[:, 6, 13] = B_6x12[:, 0, 6]
        return B
