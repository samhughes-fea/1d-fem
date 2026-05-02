"""Regression: TL geometric stiffness is a Gauss sum equivalent to legacy beam-column form for EB Hermite."""

import numpy as np
from numpy.polynomial.legendre import leggauss

from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)


def _plane_k4(L: float) -> np.ndarray:
    return np.array(
        [
            [36.0, 3.0 * L, -36.0, 3.0 * L],
            [3.0 * L, 4.0 * L * L, -3.0 * L, -L * L],
            [-36.0, -3.0 * L, 36.0, -3.0 * L],
            [3.0 * L, -L * L, -3.0 * L, 4.0 * L * L],
        ],
        dtype=np.float64,
    )


def _dN_dxi_eb_xy(L: float, xi: float) -> np.ndarray:
    """Match linear EB shape_functions (xy bending rows 1,5,7,11)."""
    dN = np.zeros((12, 6))
    dN[[0, 6], 0] = 0.5 * np.array([-1.0, 1.0])
    dN1_dxi = -0.75 + 0.75 * xi**2
    dN3_dxi = 0.75 - 0.75 * xi**2
    dN[[1, 7], 1] = [dN1_dxi, dN3_dxi]
    dN2_dxi = (L / 8) * (-1 - 2 * xi + 3 * xi**2)
    dN4_dxi = -(L / 8) * (1 - 2 * xi - 3 * xi**2)
    dN[[5, 11], 5] = [dN2_dxi, dN4_dxi]
    dN[[2, 8], 2] = dN[[1, 7], 1]
    dN[[4, 10], 4] = -dN[[5, 11], 5]
    dN[[3, 9], 3] = dN[[0, 6], 0]
    return dN


def test_geometric_stiffness_gauss_matches_legacy_constant_forces():
    L = 2.7
    detJ = L / 2.0
    dxi_dx = 2.0 / L
    xi, w = leggauss(8)
    dN_dx_stack = np.stack([_dN_dxi_eb_xy(L, float(xg)) * dxi_dx for xg in xi], axis=0)
    n_g = len(xi)
    Ntest, Mytest, Mztest = 1.25e3, 4.5e2, 6.0e2
    N_gp = np.full(n_g, Ntest)
    M_y_gp = np.full(n_g, Mytest)
    M_z_gp = np.full(n_g, Mztest)

    op = GeometricStiffnessOperator(element_length=L)
    K_new = op.assemble_K_sigma(N_gp, M_y_gp, M_z_gp, w, dN_dx_stack, detJ)

    K_old = np.zeros((12, 12), dtype=np.float64)
    cN = Ntest / (30.0 * L)
    cMy = Mytest / (30.0 * L * L)
    cMz = Mztest / (30.0 * L * L)
    K4 = _plane_k4(L)

    def _embed(K: np.ndarray, idx, K4b: np.ndarray) -> None:
        for a, ia in enumerate(idx):
            for b, ib in enumerate(idx):
                K[ia, ib] += K4b[a, b]

    for k, wk in enumerate(w):
        dN = dN_dx_stack[k]
        for i in (0, 6):
            for j in (0, 6):
                K_old[i, j] += Ntest * dN[i, 0] * dN[j, 0] * wk * detJ
    _embed(K_old, [1, 5, 7, 11], (cN + cMz) * K4)
    _embed(K_old, [2, 4, 8, 10], (cN + cMy) * K4)
    K_old = 0.5 * (K_old + K_old.T)

    np.testing.assert_allclose(K_new, K_old, rtol=1e-10, atol=1e-9)


def test_geometric_stiffness_symmetric():
    L = 1.0
    detJ = L / 2.0
    dxi_dx = 2.0 / L
    xi, w = leggauss(3)
    dN_dx_stack = np.stack([_dN_dxi_eb_xy(L, float(xg)) * dxi_dx for xg in xi], axis=0)
    rng = np.random.default_rng(0)
    n_g = len(xi)
    Ks = GeometricStiffnessOperator(element_length=L).assemble_K_sigma(
        rng.standard_normal(n_g),
        rng.standard_normal(n_g),
        rng.standard_normal(n_g),
        w,
        dN_dx_stack,
        detJ,
    )
    np.testing.assert_allclose(Ks, Ks.T, rtol=0, atol=1e-14)
