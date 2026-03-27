# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/shape_functions.py
"""Extend registry EB shape tensor (12, 6) to (14, 6) for warping χ DOFs (mass and kinematics).

The shape-function registry returns ``N``, ``∂N/∂ξ``, ``∂²N/∂ξ²`` with batch shape (n_gp, 12, 6) for
``LinearEulerBernoulliBeamElement3D`` — one row per global DOF (0–11), one column per displacement component
(0–5). This module defines the **extended** ``N`` slice (14, 6) used for **consistent mass** on 14 DOFs: rows 0–11
copy the EB block; rows 12–13 interpolate χ in the **θ_x component column** (index 3) with linear Lagrange in ξ.

Stiffness ``B`` is **not** built from this extended ``N``; it uses ``StrainDisplacementOperator`` + warping row
(see ``B_matrix.py``). See ``docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`` for column semantics.
"""

import numpy as np

from .constants import N_DOF, N_STANDARD_DOF


def extend_natural_shape_to_warping(N12: np.ndarray, xi: float) -> np.ndarray:
    """
    Extend standard (12, 6) ``N`` at one Gauss point to (14, 6).

    Rows ``0:N_STANDARD_DOF`` copy ``N12``. Rows 12 and 13 set χ¹ and χ² shape functions in the ``θ_x`` column
    (component index 3): ``L₁(ξ) = ½(1−ξ)``, ``L₂(ξ) = ½(1+ξ)``, matching the same axial/torsion linear Lagrange
    convention as the underlying EB operator.

    Parameters
    ----------
    N12 : np.ndarray
        Shape (12, 6), natural-coordinate shape functions from ``ShapeFunctionOperator.natural_coordinate_form``
        for one point.
    xi : float
        Natural coordinate in [-1, 1].

    Returns
    -------
    np.ndarray
        Extended ``N``, shape (14, 6).

    Notes
    -----
    Used in ``element_mass_matrix`` with per-DOF weights including ``ρ·Γ`` on χ DOFs. Not used for
    ``K_e`` assembly (stiffness uses ``B`` from ``WarpingStrainDisplacementOperator``).

    **Sparsity (χ rows)**

    ```text
    N_ext[12, 3] = L₁(ξ),   N_ext[13, 3] = L₂(ξ);   other entries in rows 12–13 are zero.
    ```

    See Also
    --------
    linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.shape_functions.ShapeFunctionOperator
        Base (12, 6) / (n_gp, 12, 6) tensors.
    """
    Nf = np.zeros((N_DOF, 6), dtype=np.float64)
    Nf[:N_STANDARD_DOF, :] = N12
    xi_f = float(xi)
    Nf[12, 3] = 0.5 * (1.0 - xi_f)
    Nf[13, 3] = 0.5 * (1.0 + xi_f)
    return Nf
