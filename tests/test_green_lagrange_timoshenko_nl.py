"""Tests for Timoshenko Green–Lagrange nonlinear strain operator (TL)."""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture
def L() -> float:
    return 1.0


def test_timoshenko_E_nl_axial_matches_euler_bernoulli(L: float) -> None:
    """Row 0 (axial) of E_nl matches EB utility; κ rows differ (rotation-based Timoshenko)."""
    from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.green_lagrange_strain import (
        GreenLagrangeStrainOperator as EBOperator,
    )
    from pre_processing.element_library.nonlinear.timoshenko.utilities.green_lagrange_strain import (
        GreenLagrangeStrainOperator as TimoOperator,
    )

    rng = np.random.default_rng(42)
    dN_dx = rng.standard_normal((12, 6))
    d2N_dx2 = rng.standard_normal((12, 6))
    u_e = rng.standard_normal(12)

    eb = EBOperator(element_length=L, include_shear=False)
    timo = TimoOperator(element_length=L, include_shear=False)

    e_eb = eb.strain_nonlinear_part(dN_dx, d2N_dx2, u_e)
    e_timo = timo.strain_nonlinear_part(dN_dx, d2N_dx2, u_e, N=None)

    np.testing.assert_allclose(e_timo[0], e_eb[0], rtol=1e-15, atol=1e-15)


def test_shear_nl_algebra_and_B_nl_fd(L: float) -> None:
    """Centroid shear NL matches closed form; B_nl columns ~ finite difference of E_nl."""
    from pre_processing.element_library.nonlinear.timoshenko.utilities.green_lagrange_strain import (
        GreenLagrangeStrainOperator,
    )

    op = GreenLagrangeStrainOperator(element_length=L, include_shear=True)
    dN_dx = np.zeros((12, 6), dtype=np.float64)
    dN_dx[0, 0] = 2.0
    dN_dx[6, 0] = 0.5
    dN_dx[1, 1] = -1.0
    dN_dx[7, 1] = 1.0
    dN_dx[2, 2] = 0.25
    dN_dx[8, 2] = 0.75

    d2N_dx2 = np.zeros((12, 6), dtype=np.float64)
    N = np.zeros((12, 6), dtype=np.float64)
    N[3, 3] = 1.0
    N[9, 3] = 0.0
    N[4, 4] = 0.5
    N[10, 4] = 0.5
    N[5, 5] = 1.0
    N[11, 5] = 0.0

    u_e = np.zeros(12, dtype=np.float64)
    u_e[0] = 1.0
    u_e[3] = 0.1
    u_e[4] = 0.2
    u_e[5] = 0.3
    u_e[10] = 0.2

    du_dx = dN_dx[0, 0] * u_e[0] + dN_dx[6, 0] * u_e[6]
    dv_dx = dN_dx[1, 1] * u_e[1] + dN_dx[7, 1] * u_e[7]
    dw_dx = dN_dx[2, 2] * u_e[2] + dN_dx[8, 2] * u_e[8]
    theta_x = N[3, 3] * u_e[3] + N[9, 3] * u_e[9]
    theta_y = N[4, 4] * u_e[4] + N[10, 4] * u_e[10]
    theta_z = N[5, 5] * u_e[5] + N[11, 5] * u_e[11]

    expected_xy = -theta_z * du_dx + theta_x * dw_dx
    expected_xz = theta_y * du_dx - theta_x * dv_dx

    e_nl = op.strain_nonlinear_part(dN_dx, d2N_dx2, u_e, N)
    assert e_nl[3] == pytest.approx(expected_xy)
    assert e_nl[4] == pytest.approx(expected_xz)

    B_nl = op.nonlinear_strain_displacement_gradient(dN_dx, d2N_dx2, u_e, N)
    h = 1e-8
    for j in range(12):
        up = u_e.copy()
        up[j] += h
        e_p = op.strain_nonlinear_part(dN_dx, d2N_dx2, up, N)
        col_fd = (e_p - e_nl) / h
        np.testing.assert_allclose(B_nl[:, j], col_fd, rtol=5e-5, atol=5e-5)


def test_include_shear_false_skips_shear_nl_even_with_N(L: float) -> None:
    from pre_processing.element_library.nonlinear.timoshenko.utilities.green_lagrange_strain import (
        GreenLagrangeStrainOperator,
    )

    op = GreenLagrangeStrainOperator(element_length=L, include_shear=False)
    dN_dx = np.random.default_rng(0).standard_normal((12, 6))
    d2N_dx2 = np.random.default_rng(1).standard_normal((12, 6))
    u_e = np.random.default_rng(2).standard_normal(12)
    N = np.ones((12, 6)) * 0.1

    e_nl = op.strain_nonlinear_part(dN_dx, d2N_dx2, u_e, N)
    assert e_nl[3] == 0.0
    assert e_nl[4] == 0.0
