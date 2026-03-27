# pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli_with_warp/utilities/shape_functions.py
"""Extend registry EB shape tensor (12, 6) to (14, 6) for warping χ DOFs."""

import numpy as np

from .constants import N_DOF, N_STANDARD_DOF


def extend_natural_shape_to_warping(N12: np.ndarray, xi: float) -> np.ndarray:
    """
    Extend standard (12, 6) N at one Gauss point to (14, 6).

    Rows 0–11 copy ``N12``; rows 12–13 use linear Lagrange on ξ in the θ_x component
    column (index 3), matching the twist/warping kinematics used for mass and load maps.

    Parameters
    ----------
    N12 : np.ndarray
        Shape (12, 6), natural-coordinate shape functions from EB ``ShapeFunctionOperator``.
    xi : float
        Natural coordinate in [-1, 1].

    Returns
    -------
    np.ndarray
        Shape (14, 6).
    """
    Nf = np.zeros((N_DOF, 6), dtype=np.float64)
    Nf[:N_STANDARD_DOF, :] = N12
    xi_f = float(xi)
    Nf[12, 3] = 0.5 * (1.0 - xi_f)
    Nf[13, 3] = 0.5 * (1.0 + xi_f)
    return Nf
