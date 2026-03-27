"""
Verify shape function formulation against standard Hermite cubic polynomials.

Standard Euler-Bernoulli beam element shape functions should satisfy:
- N1(-1) = 1, N1(1) = 0, N1'(-1) = 0, N1'(1) = 0
- N2(-1) = 0, N2(1) = 0, N2'(-1) = 1, N2'(1) = 0
- N3(-1) = 0, N3(1) = 1, N3'(-1) = 0, N3'(1) = 0
- N4(-1) = 0, N4(1) = 0, N4'(-1) = 0, N4'(1) = 1

Where N1, N3 are displacement shape functions and N2, N4 are rotation shape functions.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.shape_functions import ShapeFunctionOperator

def verify_shape_functions():
    """Verify shape functions match standard Hermite cubic formulation."""
    
    L = 0.2
    shape_op = ShapeFunctionOperator(element_length=L)
    
    # Evaluate at nodes: xi = -1 (node 1) and xi = 1 (node 2)
    xi_nodes = np.array([-1.0, 1.0])
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(xi_nodes)
    
    print("=" * 70)
    print("SHAPE FUNCTION VERIFICATION")
    print("=" * 70)
    print(f"Element length L = {L:.3f} m")
    print()
    
    # Check displacement shape functions (u_y)
    # N1 = shape function for u_y1 (DOF 1)
    # N3 = shape function for u_y2 (DOF 7)
    print("=== Displacement Shape Functions (u_y) ===")
    print("At xi = -1 (node 1):")
    print(f"  N[1,1] (u_y1) = {N[0,1,1]:.6f} (should be 1.0)")
    print(f"  N[7,1] (u_y2) = {N[0,7,1]:.6f} (should be 0.0)")
    print(f"  dN_dxi[1,1] = {dN_dxi[0,1,1]:.6f} (should be 0.0)")
    print(f"  dN_dxi[7,1] = {dN_dxi[0,7,1]:.6f}")
    print()
    
    print("At xi = 1 (node 2):")
    print(f"  N[1,1] (u_y1) = {N[1,1,1]:.6f} (should be 0.0)")
    print(f"  N[7,1] (u_y2) = {N[1,7,1]:.6f} (should be 1.0)")
    print(f"  dN_dxi[1,1] = {dN_dxi[1,1,1]:.6f}")
    print(f"  dN_dxi[7,1] = {dN_dxi[1,7,1]:.6f} (should be 0.0)")
    print()
    
    # Check rotation shape functions (theta_z)
    # N2 = shape function for theta_z1 (DOF 5)
    # N4 = shape function for theta_z2 (DOF 11)
    print("=== Rotation Shape Functions (theta_z) ===")
    print("At xi = -1 (node 1):")
    print(f"  N[5,5] (theta_z1) = {N[0,5,5]:.6f} (should be 0.0)")
    print(f"  N[11,5] (theta_z2) = {N[0,11,5]:.6f} (should be 0.0)")
    print(f"  dN_dxi[5,5] = {dN_dxi[0,5,5]:.6f} (should be 1.0 for standard, or L for some formulations)")
    print(f"  dN_dxi[11,5] = {dN_dxi[0,11,5]:.6f} (should be 0.0)")
    print()
    
    print("At xi = 1 (node 2):")
    print(f"  N[5,5] (theta_z1) = {N[1,5,5]:.6f} (should be 0.0)")
    print(f"  N[11,5] (theta_z2) = {N[1,11,5]:.6f} (should be 0.0)")
    print(f"  dN_dxi[5,5] = {dN_dxi[1,5,5]:.6f} (should be 0.0)")
    print(f"  dN_dxi[11,5] = {dN_dxi[1,11,5]:.6f} (should be 1.0 for standard, or L for some formulations)")
    print()
    
    # Check polynomial expressions
    print("=== Polynomial Expressions ===")
    xi = 0.0
    N_test, dN_dxi_test, d2N_dxi2_test = shape_op.natural_coordinate_form(np.array([xi]))
    N_test = N_test[0]
    dN_dxi_test = dN_dxi_test[0]
    d2N_dxi2_test = d2N_dxi2_test[0]
    
    print(f"At xi = 0.0:")
    print(f"  N[1,1] (u_y1) = {N_test[1,1]:.6f}")
    print(f"  N[7,1] (u_y2) = {N_test[7,1]:.6f}")
    print(f"  N[5,5] (theta_z1) = {N_test[5,5]:.6f}")
    print(f"  N[11,5] (theta_z2) = {N_test[11,5]:.6f}")
    print()
    
    # Standard Hermite cubics (for comparison)
    # N1 = (1/4)(1-xi)^2(2+xi) = 0.5 - 0.75*xi + 0.25*xi^3
    # N2 = (L/8)(1-xi)^2(1+xi) = (L/8)(1 - xi - xi^2 + xi^3)
    # N3 = (1/4)(1+xi)^2(2-xi) = 0.5 + 0.75*xi - 0.25*xi^3
    # N4 = -(L/8)(1+xi)^2(1-xi) = -(L/8)(1 + xi - xi^2 - xi^3)
    
    N1_std = 0.5 - 0.75*xi + 0.25*xi**3
    N2_std = (L/8) * (1 - xi - xi**2 + xi**3)
    N3_std = 0.5 + 0.75*xi - 0.25*xi**3
    N4_std = -(L/8) * (1 + xi - xi**2 - xi**3)
    
    print("Standard Hermite cubics (for comparison):")
    print(f"  N1_std = {N1_std:.6f}")
    print(f"  N2_std = {N2_std:.6f}")
    print(f"  N3_std = {N3_std:.6f}")
    print(f"  N4_std = {N4_std:.6f}")
    print()
    
    # Code uses different formulation
    N1_code = 1 - 3*xi**2 + 2*xi**3
    N2_code = xi - 2*xi**2 + xi**3
    N3_code = 3*xi**2 - 2*xi**3
    N4_code = -xi**2 + xi**3
    
    print("Code formulation:")
    print(f"  N1_code = 1 - 3*xi^2 + 2*xi^3 = {N1_code:.6f}")
    print(f"  N2_code = xi - 2*xi^2 + xi^3 = {N2_code:.6f}")
    print(f"  N3_code = 3*xi^2 - 2*xi^3 = {N3_code:.6f}")
    print(f"  N4_code = -xi^2 + xi^3 = {N4_code:.6f}")
    print()
    
    print("These are DIFFERENT formulations! The code uses a different normalization.")
    print("This may be the source of the scaling error.")

if __name__ == "__main__":
    verify_shape_functions()

