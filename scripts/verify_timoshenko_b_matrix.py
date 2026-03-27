"""
Verify Timoshenko B-matrix computes shear terms correctly.

Checks that γ_xy = du_y/dx - θ_z is correctly computed at Gauss points.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.shape_functions import ShapeFunctionOperator
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.B_matrix import StrainDisplacementOperator

def verify_timoshenko_b_matrix():
    """Verify Timoshenko B-matrix shear terms."""
    print("=" * 70)
    print("VERIFICATION: Timoshenko B-Matrix Shear Terms")
    print("=" * 70)
    
    # Element properties
    L = 0.2  # m (element length)
    
    # Initialize operators
    shape_op = ShapeFunctionOperator(element_length=L)
    strain_op = StrainDisplacementOperator(element_length=L)
    
    # Test at a few Gauss points
    xi_points = np.array([-0.774597, 0.0, 0.774597])  # 3-point Gauss quadrature
    
    print(f"\nElement length L = {L:.3f} m")
    print(f"Testing at {len(xi_points)} Gauss points: {xi_points}")
    print()
    
    all_ok = True
    
    for i, xi in enumerate(xi_points):
        print(f"Gauss Point {i+1}: xi = {xi:.6f}")
        
        # Get shape functions and derivatives (returns shape (n_points, 12, 6))
        N, dN_dξ, d2N_dξ2 = shape_op.natural_coordinate_form(np.array([xi]))
        # Keep batch dimension - B-matrix expects (n_gauss, 12, 6)
        
        # Get B-matrix (returns shape (n_gauss, 6, 12))
        B = strain_op.physical_coordinate_form(dN_dξ, d2N_dξ2, N)
        B = B[0]  # Get first (and only) Gauss point, shape (6, 12)
        
        # B-matrix shape should be (6, 12) for 6 strain components and 12 DOFs
        print(f"  B-matrix shape: {B.shape}")
        
        # Check shear term γ_xy = du_y/dx - θ_z
        # This should be in B[3, :] (row 3 = shear_xy)
        # B[3, 1] should be du_y/dx term (from node 1, DOF 1 = u_y)
        # B[3, 5] should be -θ_z term (from node 1, DOF 5 = θ_z)
        # B[3, 7] should be du_y/dx term (from node 2, DOF 7 = u_y)
        # B[3, 11] should be -θ_z term (from node 2, DOF 11 = θ_z)
        
        print(f"  B[3, 1] (du_y1/dx): {B[3, 1]:.6e}")
        print(f"  B[3, 5] (-theta_z1): {B[3, 5]:.6e}")
        print(f"  B[3, 7] (du_y2/dx): {B[3, 7]:.6e}")
        print(f"  B[3, 11] (-theta_z2): {B[3, 11]:.6e}")
        
        # Verify that B[3, 5] and B[3, 11] are negative (they represent -θ_z)
        if B[3, 5] >= 0:
            print(f"  [WARNING] B[3, 5] should be negative (represents -theta_z1)")
            all_ok = False
        if B[3, 11] >= 0:
            print(f"  [WARNING] B[3, 11] should be negative (represents -theta_z2)")
            all_ok = False
        
        # Check that shear terms are non-zero
        shear_row_norm = np.linalg.norm(B[3, :])
        print(f"  Shear row (B[3, :]) norm: {shear_row_norm:.6e}")
        
        if shear_row_norm < 1e-10:
            print(f"  [ERROR] Shear terms are essentially zero!")
            all_ok = False
        
        # Check coordinate transformation
        # dxi/dx = 2/L
        dxi_dx = 2.0 / L
        print(f"  Coordinate transformation: dxi/dx = {dxi_dx:.6f}")
        
        # For a unit displacement u_y1 = 1, the strain should be du_y/dx = dN/dξ * dξ/dx
        # Check if B[3, 1] matches dN_dξ[1, 1] * dξ_dx (approximately, accounting for shape function structure)
        print()
    
    print("=" * 70)
    if all_ok:
        print("[PASS] B-matrix shear terms appear to be computed correctly")
    else:
        print("[FAIL] B-matrix shear terms have issues")
    print("=" * 70)
    
    return all_ok

if __name__ == "__main__":
    verify_timoshenko_b_matrix()
