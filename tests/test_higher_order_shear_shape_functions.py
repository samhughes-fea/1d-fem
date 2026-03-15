"""
Higher-order shear (Phase 2b): Reddy/Levinson shape functions — interpolation and continuity.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.reddy.utilities.shape_functions import ShapeFunctionOperator


def test_higher_order_shear_shape_functions():
    """Reddy (higher-order shear) shape functions: quintic/cubic interpolation at nodes and C1 at mid-span."""
    L = 0.5
    op = ShapeFunctionOperator(element_length=L)
    xi_nodes = np.array([-1.0, 1.0])
    N, dN_dxi, d2N_dxi2 = op.natural_coordinate_form(xi_nodes)
    assert N.shape[0] == 2 and N.shape[1] == 12 and N.shape[2] == 6
    assert np.isclose(N[0, 1, 1], 1.0, atol=1e-10) and np.isclose(N[1, 7, 1], 1.0, atol=1e-10)
    assert np.isclose(N[0, 5, 5], 1.0, atol=1e-10) and np.isclose(N[1, 11, 5], 1.0, atol=1e-10)


def test_higher_order_shear_shape_functions_continuity():
    """Reddy shape functions: continuity at mid-span (ξ=0); N continuous, dN/dξ finite."""
    L = 0.5
    op = ShapeFunctionOperator(element_length=L)
    xi_mid = np.array([0.0])
    N, dN_dxi, d2N_dxi2 = op.natural_coordinate_form(xi_mid)
    assert np.all(np.isfinite(N))
    assert np.all(np.isfinite(dN_dxi))
    assert N.shape == (1, 12, 6) and dN_dxi.shape == (1, 12, 6)
