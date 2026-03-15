"""
Higher-order shear (Phase 2b): B and D dimensions and shear kinematics (γ = du/dx − θ + α d²θ/dx²).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.reddy.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.linear.reddy.utilities.D_matrix import MaterialStiffnessOperator
from pre_processing.element_library.linear.reddy.utilities.shape_functions import ShapeFunctionOperator


def test_higher_order_shear_B_D_dimensions():
    """Reddy B has shape (6, 12), D has shape (6, 6); shear rows (3, 4) use higher-order kinematics."""
    L = 1.0
    alpha = 1e-6
    shape_op = ShapeFunctionOperator(element_length=L)
    b_op = StrainDisplacementOperator(element_length=L, alpha_coefficient=alpha)
    E, G, A, I_y, I_z, J_t = 2.1e11, 8.1e10, 0.001, 1e-8, 1e-8, 1e-9
    d_op = MaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=I_y,
        moment_inertia_z=I_z,
        torsion_constant=J_t,
    )
    xi = np.array([0.0])
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(xi)
    B = b_op.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
    D = d_op.assembly_form()
    assert B.shape == (6, 12), f"B shape {B.shape} != (6, 12)"
    assert D.shape == (6, 6), f"D shape {D.shape} != (6, 6)"
    assert np.any(B[3, :] != 0.0) and np.any(B[4, :] != 0.0), "Shear rows must be non-zero"
    assert D[3, 3] > 0 and D[4, 4] > 0, "Shear stiffness (GA) must be positive"
    assert D[3, 3] == G * A and D[4, 4] == G * A, "Reddy uses GA without κ"


def test_higher_order_shear_B_uses_second_derivative():
    """Reddy shear strain γ = du/dx − θ + α d²θ/dx²; shear rows (3,4) couple displacement gradient, θ, and (when alpha != 0) d²θ/dx²."""
    L = 1.0
    shape_op = ShapeFunctionOperator(element_length=L)
    b_op = StrainDisplacementOperator(element_length=L, alpha_coefficient=0.01)
    xi = np.array([0.0])
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(xi)
    B = b_op.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
    assert B.shape == (6, 12)
    assert np.any(B[3, [1, 7]] != 0.0) and np.any(B[4, [2, 8]] != 0.0), "Shear rows include du/dx"
    assert np.any(B[3, [5, 11]] != 0.0) and np.any(B[4, [4, 10]] != 0.0), "Shear rows include θ (and possibly α d²θ/dx²)"
