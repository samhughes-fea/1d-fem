# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/B_matrix.py
"""Strain–displacement B (7, 14) per Gauss point for 2-node 3-D Euler–Bernoulli beam with Vlasov warping.

At each Gauss point, ``ε = B U_e`` with ``ε ∈ ℝ^7`` and ``U_e ∈ ℝ^{14}``; rows 0–5 match
``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`` / linear EB; row 6 is the warping-extension strain.
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

    **Constitutive hook**

    ``S = D ε`` with ``S ∈ ℝ^7`` and ``D ∈ ℝ^{7×7}`` from
    ``euler_bernoulli_with_warp/utilities/D_matrix.py``. Rows 0–5 pair with the usual section
    resultants; row 6 pairs with the warping resultant ``S_w`` (notation as in that module).

    Parameters
    ----------
    element_length : float
        Chord length ``L`` of the beam element [m]; used for row 6 coefficients ``±1/L`` and
        passed consistently with the base EB operator (same element as ``base_strain_operator``).
    base_strain_operator : StrainDisplacementOperator
        Linear EB operator that maps ``(dN/dξ, d²N/dξ²)`` to ``B`` of shape ``(n_gp, 6, 12)`` in
        physical ``x``. Its ``physical_coordinate_form`` supplies rows 0–5 and columns 0–11;
        ``dξ_dx`` and ``d²ξ/dx²`` live on that object (``2/L``, ``4/L²``).

    Notes
    -----
    **Isoparametric mapping**

    The same map ``x(ξ) = L(1+ξ)/2`` applies as for linear EB; ``∂ξ/∂x = 2/L``, ``∂²ξ/∂x² = 4/L²``.
    Rows 0–5 are assembled by ``base_strain_operator.physical_coordinate_form``, which applies
    the chain rule to the registry Hermite/linear shape-function derivatives. Row 6 is **not**
    built from ``dN/dξ`` at ``ξ_g``: it is filled once per row in physical ``x`` with constant
    entries ``±1/L`` on ``θ_x`` and ``χ`` DOFs (2-node linear Lagrange derivatives in ``x``).

    **Shape function basis**

    Rows 0–5 use the same ``N``, ``∂N/∂ξ``, ``∂²N/∂ξ²`` pipeline and column coupling as
    ``StrainDisplacementOperator`` (see ``euler_bernoulli/utilities/B_matrix.py``). Row 6 uses
    linear Lagrange ``L₁, L₂`` on ``ξ`` for both ``θ_x`` (columns 3 and 9) and ``χ`` (columns 12
    and 13), i.e. the same 2-node pattern as axial/torsion in EB. For **consistent mass**, χ
    interpolation in the extended ``N`` tensor is implemented in
    ``shape_functions.extend_natural_shape_to_warping``.

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
    [ 0     0     0    −1/L   0     0     0    0     0    +1/L    0      0    −1/L +1/L ]  # φ_x′ (row 6)
    ```

    **Non-zero entries of B in physical coordinates**

    Rows 0–5: identical formulas to ``euler_bernoulli/utilities/B_matrix.py`` (section *Non-zero
    entries of B in physical coordinates*); columns 12–13 are zero.

    Row 6 (warping strain rate ``φ_x′ = ∂θ_x/∂x + ∂χ/∂x``), constant in ``ξ``:

    ```text
      B[6,3]  = −1/L    (θ_x, node 1)
      B[6,9]  = +1/L    (θ_x, node 2)
      B[6,12] = −1/L    (χ¹)
      B[6,13] = +1/L    (χ²)
    ```

    **Batch shape**

    ``physical_coordinate_form`` returns ``(n_gp, 7, 14)``. Rows 0–5 vary with ``ξ_g``; row 6 is
    the same for every Gauss point.

    Weak-form assembly: ``K_e += B.T @ D @ B * w_g * detJ`` with ``detJ = L/2``.

    See Also
    --------
    StrainDisplacementOperator
        Linear EB ``B`` (6, 12) in ``euler_bernoulli/utilities/B_matrix.py``.
    WarpingMaterialStiffnessOperator
        ``D`` (7, 7) in ``euler_bernoulli_with_warp/utilities/D_matrix.py``.
    linear_warping_euler_bernoulli_3D.LinearWarpingEulerBernoulliBeamElement3D
        Parent element using this operator.
    docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md
        Voigt order for rows 0–5; extensions for ``(14,) U_e`` and ``(7, 14) B``.
    """

    element_length: float
    base_strain_operator: StrainDisplacementOperator

    def warping_row(self) -> np.ndarray:
        """
        Row 6 of ``B``: ``φ_x′ = ∂θ_x/∂x + ∂χ/∂x`` with 2-node linear interpolation in ``x``
        (coefficients ``±1/L`` on ``θ_x`` and ``χ`` DOFs).

        Returns
        -------
        np.ndarray
            Shape ``(14,)``.
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
            First derivatives ``∂N/∂ξ``, shape ``(n_gauss, 12, 6)``.
        d2N_dξ2 : np.ndarray
            Second derivatives ``∂²N/∂ξ²``, shape ``(n_gauss, 12, 6)``.

        Returns
        -------
        np.ndarray
            Strain–displacement batch, shape ``(n_gauss, 7, 14)``.
        """
        n_gauss = dN_dξ.shape[0]
        B_6x12 = self.base_strain_operator.physical_coordinate_form(dN_dξ, d2N_dξ2)
        B = np.zeros((n_gauss, N_STRAIN, N_DOF), dtype=np.float64)
        B[:, :6, :12] = B_6x12
        wr = self.warping_row()
        for g in range(n_gauss):
            B[g, 6, :] = wr
        return B
