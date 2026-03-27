"""
Verify curved beam: at κ0=0, B-matrix matches straight Timoshenko B (Phase 3b regression).
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def verify_curved_beam_straight_limit():
    """At κ0=0, curved strain-displacement B equals straight Timoshenko B at a Gauss point."""
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.curved_beam.utilities.B_matrix import CurvedStrainDisplacementOperator
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.B_matrix import StrainDisplacementOperator
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.shape_functions import ShapeFunctionOperator

    L = 1.0
    xi = np.array([0.0])
    shape_op = ShapeFunctionOperator(element_length=L)
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(xi)
    B_straight = StrainDisplacementOperator(element_length=L).physical_coordinate_form(
        dN_dxi, d2N_dxi2, N
    )[0]
    B_curved = CurvedStrainDisplacementOperator(element_length=L, curvature=0.0).physical_coordinate_form(
        dN_dxi, d2N_dxi2, N
    )[0]
    if not np.allclose(B_curved, B_straight, atol=1e-12):
        print("FAIL: kappa0=0 curved B should equal straight Timoshenko B")
        return 1
    print("PASS: Curved beam straight limit (kappa0=0) B matches Timoshenko B")
    return 0


if __name__ == "__main__":
    sys.exit(verify_curved_beam_straight_limit())
