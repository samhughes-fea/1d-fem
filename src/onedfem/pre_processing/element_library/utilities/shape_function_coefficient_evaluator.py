# pre_processing/element_library/utilities/shape_function_coefficient_evaluator.py

"""
Generic evaluator for shape functions from polynomial coefficients (B2 format).

Coefficient convention: N_coefficients[dof, comp, k] = coefficient of ξ^k,
so N_dof,comp(ξ) = sum_k N_coefficients[dof, comp, k] * ξ**k.
Shape (12, 6, 4) for 1D beam with max degree 3.
"""

from typing import Tuple
import numpy as np


def evaluate_shape_functions_from_coefficients(
    N_coeffs: np.ndarray,
    dN_dxi_coeffs: np.ndarray,
    d2N_dxi2_coeffs: np.ndarray,
    xi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluate N(ξ), dN/dξ, d²N/dξ² from monomial coefficients.

    Parameters
    ----------
    N_coeffs : np.ndarray
        Shape (12, 6, 4). N_coeffs[dof, comp, k] = coefficient of ξ^k for N.
    dN_dxi_coeffs : np.ndarray
        Shape (12, 6, 4). Coefficients for dN/dξ.
    d2N_dxi2_coeffs : np.ndarray
        Shape (12, 6, 4). Coefficients for d²N/dξ².
    xi : np.ndarray
        Natural coordinates, shape (n_points,) or broadcastable.

    Returns
    -------
    N : np.ndarray
        Shape (n_points, 12, 6). Shape function values.
    dN_dxi : np.ndarray
        Shape (n_points, 12, 6). First derivatives.
    d2N_dxi2 : np.ndarray
        Shape (n_points, 12, 6). Second derivatives.
    """
    xi = np.asarray(xi, dtype=np.float64)
    if xi.ndim == 0:
        xi = xi.reshape(1)
    n_points = xi.size
    xi_flat = xi.ravel()

    N = np.zeros((n_points, 12, 6), dtype=np.float64)
    dN_dxi = np.zeros((n_points, 12, 6), dtype=np.float64)
    d2N_dxi2 = np.zeros((n_points, 12, 6), dtype=np.float64)

    # np.polyval expects coeffs from highest degree to constant; we store [c0, c1, c2, c3]
    for dof in range(12):
        for comp in range(6):
            p_N = N_coeffs[dof, comp, :][::-1]
            p_dN = dN_dxi_coeffs[dof, comp, :][::-1]
            p_d2N = d2N_dxi2_coeffs[dof, comp, :][::-1]
            N[:, dof, comp] = np.polyval(p_N, xi_flat)
            dN_dxi[:, dof, comp] = np.polyval(p_dN, xi_flat)
            d2N_dxi2[:, dof, comp] = np.polyval(p_d2N, xi_flat)

    return N, dN_dxi, d2N_dxi2
