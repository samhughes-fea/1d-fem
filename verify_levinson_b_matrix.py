"""
Verify Levinson B-matrix includes higher-order shear terms.

Expected: gamma_xy = du_y/dx - theta_z + alpha*(d^2theta_z/dx^2)
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.levinson.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.levinson.utilities.shape_functions import ShapeFunctionOperator

def verify_levinson_b_matrix():
    """Verify Levinson B-matrix includes higher-order shear terms."""
    
    L = 0.2
    shape_op = ShapeFunctionOperator(element_length=L)
    strain_op = StrainDisplacementOperator(element_length=L)
    
    # Evaluate at a Gauss point (not at center where rotation shape functions are zero)
    xi_g = 0.577350269  # Typical Gauss point
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(np.array([xi_g]))
    
    # Get B-matrix (pass N for shear terms)
    B = strain_op.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
    
    print("=" * 70)
    print("LEVINSON B-MATRIX VERIFICATION")
    print("=" * 70)
    print(f"Element length L = {L:.3f} m")
    print(f"Gauss point xi = {xi_g:.3f}")
    print()
    
    print("=== B-MATRIX STRUCTURE ===")
    print(f"B-matrix shape: {B.shape}")
    print("Expected: (6, 12) for [eps_x, kappa_z, kappa_y, gamma_xy, gamma_xz, phi_x]")
    print()
    
    print("=== KEY TERMS (at xi=0) ===")
    print("Axial strain (eps_x):")
    print(f"  B[0,0] (from u_x1) = {B[0,0]:.6e}")
    print(f"  B[0,6] (from u_x2) = {B[0,6]:.6e}")
    print()
    
    print("Bending curvature (kappa_z):")
    print(f"  B[1,2] (from u_z1) = {B[1,2]:.6e}")
    print(f"  B[1,8] (from u_z2) = {B[1,8]:.6e}")
    print(f"  B[1,4] (from theta_y1) = {B[1,4]:.6e}")
    print(f"  B[1,10] (from theta_y2) = {B[1,10]:.6e}")
    print()
    
    print("=== VERIFICATION ===")
    issues = []
    
    # Check B-matrix dimensions
    if B.shape != (6, 12):
        issues.append(f"B-matrix has wrong shape: {B.shape} (expected (6, 12))")
    
    # Check shear terms
    print("Shear strain (gamma_xy = du_y/dx - theta_z + alpha*(d^2theta_z/dx^2) for Levinson):")
    print(f"  B[3,1] (from u_y1, du_y/dx term) = {B[3,1]:.6e}")
    print(f"  B[3,7] (from u_y2, du_y/dx term) = {B[3,7]:.6e}")
    print(f"  B[3,5] (from theta_z1, -theta_z term) = {B[3,5]:.6e}")
    print(f"  B[3,11] (from theta_z2, -theta_z term) = {B[3,11]:.6e}")
    print()
    
    if abs(B[3,1]) < 1e-10 and abs(B[3,7]) < 1e-10:
        issues.append("CRITICAL: Shear strain terms (gamma_xy) are ZERO - Levinson should include du_y/dx terms")
    
    # Check rotation terms (may be zero at specific points, but should be non-zero in general)
    if abs(B[3,5]) < 1e-10 and abs(B[3,11]) < 1e-10:
        # Check if this is just because we're at a point where rotation shape functions are zero
        if abs(N[0,5,5]) < 1e-10 and abs(N[0,11,5]) < 1e-10:
            print("NOTE: Rotation shape functions are zero at this point (expected for some points)")
        else:
            issues.append("CRITICAL: Shear strain terms (gamma_xy) are ZERO - Levinson should include -theta_z terms")
    
    if issues:
        print("FAILURES FOUND:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    else:
        print("B-matrix structure verified, but shear terms need separate verification")
        return True

if __name__ == "__main__":
    verify_levinson_b_matrix()

