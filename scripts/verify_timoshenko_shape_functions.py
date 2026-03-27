"""
Verify Timoshenko shape functions satisfy standard Hermite cubic properties.

Timoshenko uses the same shape functions as Euler-Bernoulli for bending,
which we already fixed. This script verifies they are correct.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.shape_functions import ShapeFunctionOperator

def verify_timoshenko_shape_functions():
    """Verify Timoshenko shape functions match standard Hermite cubic formulation."""
    
    L = 0.2
    shape_op = ShapeFunctionOperator(element_length=L)
    
    # Evaluate at nodes: xi = -1 (node 1) and xi = 1 (node 2)
    xi_nodes = np.array([-1.0, 1.0])
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(xi_nodes)
    
    print("=" * 70)
    print("TIMOSHENKO SHAPE FUNCTION VERIFICATION")
    print("=" * 70)
    print(f"Element length L = {L:.3f} m")
    print()
    
    # Check displacement shape functions (u_y)
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
    print("=== Rotation Shape Functions (theta_z) ===")
    print("At xi = -1 (node 1):")
    print(f"  N[5,5] (theta_z1) = {N[0,5,5]:.6f} (should be 0.0)")
    print(f"  N[11,5] (theta_z2) = {N[0,11,5]:.6f} (should be 0.0)")
    print(f"  dN_dxi[5,5] = {dN_dxi[0,5,5]:.6f} (should be L/2 = {L/2:.6f})")
    print(f"  dN_dxi[11,5] = {dN_dxi[0,11,5]:.6f} (should be 0.0)")
    print()
    
    print("At xi = 1 (node 2):")
    print(f"  N[5,5] (theta_z1) = {N[1,5,5]:.6f} (should be 0.0)")
    print(f"  N[11,5] (theta_z2) = {N[1,11,5]:.6f} (should be 0.0)")
    print(f"  dN_dxi[5,5] = {dN_dxi[1,5,5]:.6f} (should be 0.0)")
    print(f"  dN_dxi[11,5] = {dN_dxi[1,11,5]:.6f} (should be L/2 = {L/2:.6f})")
    print()
    
    # Verify properties
    print("=== VERIFICATION ===")
    checks = []
    
    # Displacement at node 1
    checks.append(("N1(-1) = 1", abs(N[0,1,1] - 1.0) < 1e-10))
    checks.append(("N3(-1) = 0", abs(N[0,7,1]) < 1e-10))
    checks.append(("N1'(-1) = 0", abs(dN_dxi[0,1,1]) < 1e-10))
    
    # Displacement at node 2
    checks.append(("N1(1) = 0", abs(N[1,1,1]) < 1e-10))
    checks.append(("N3(1) = 1", abs(N[1,7,1] - 1.0) < 1e-10))
    checks.append(("N3'(1) = 0", abs(dN_dxi[1,7,1]) < 1e-10))
    
    # Rotation at node 1
    checks.append(("N2(-1) = 0", abs(N[0,5,5]) < 1e-10))
    checks.append(("N2'(-1) = L/2", abs(dN_dxi[0,5,5] - L/2) < 1e-10))
    
    # Rotation at node 2
    checks.append(("N4(1) = 0", abs(N[1,11,5]) < 1e-10))
    checks.append(("N4'(1) = L/2", abs(dN_dxi[1,11,5] - L/2) < 1e-10))
    
    all_passed = True
    for name, passed in checks:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False
    
    print()
    if all_passed:
        print("SUCCESS: All shape function properties verified!")
    else:
        print("FAILURE: Some shape function properties are incorrect!")
    
    return all_passed

if __name__ == "__main__":
    verify_timoshenko_shape_functions()

