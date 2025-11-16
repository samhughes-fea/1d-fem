"""
Verify Timoshenko B-matrix includes shear terms correctly.

CRITICAL ISSUE FOUND: Timoshenko B-matrix is currently identical to Euler-Bernoulli,
which means it does NOT include shear strain terms. This is incorrect for Timoshenko theory.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.timoshenko.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.timoshenko.utilities.shape_functions import ShapeFunctionOperator

def verify_timoshenko_b_matrix():
    """Verify Timoshenko B-matrix includes shear terms."""
    
    L = 0.2
    shape_op = ShapeFunctionOperator(element_length=L)
    strain_op = StrainDisplacementOperator(element_length=L)
    
    # Evaluate at a Gauss point
    xi_g = 0.0
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(np.array([xi_g]))
    
    # Get B-matrix (pass N for shear terms)
    B = strain_op.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
    
    print("=" * 70)
    print("TIMOSHENKO B-MATRIX VERIFICATION")
    print("=" * 70)
    print(f"Element length L = {L:.3f} m")
    print(f"Gauss point xi = {xi_g:.3f}")
    print()
    
    print("=== B-MATRIX STRUCTURE ===")
    print("Strain vector: [eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]")
    print("DOF vector: [u_x1, u_y1, u_z1, theta_x1, theta_y1, theta_z1, u_x2, u_y2, u_z2, theta_x2, theta_y2, theta_z2]")
    print()
    
    print("=== KEY TERMS (at xi=0) ===")
    print("Axial strain (eps_x):")
    print(f"  B[0,0] (from u_x1) = {B[0,0]:.6e}")
    print(f"  B[0,6] (from u_x2) = {B[0,6]:.6e}")
    print()
    
    print("Bending curvature (kappa_z = dtheta_z/dx for Timoshenko):")
    print(f"  B[2,5] (from theta_z1) = {B[2,5]:.6e}")
    print(f"  B[2,11] (from theta_z2) = {B[2,11]:.6e}")
    print()
    
    print("CRITICAL: Shear strain (gamma_xy = du_y/dx - theta_z for Timoshenko):")
    print(f"  B[3,1] (from u_y1, du_y/dx term) = {B[3,1]:.6e}")
    print(f"  B[3,7] (from u_y2, du_y/dx term) = {B[3,7]:.6e}")
    print(f"  B[3,5] (from theta_z1, -theta_z term) = {B[3,5]:.6e}")
    print(f"  B[3,11] (from theta_z2, -theta_z term) = {B[3,11]:.6e}")
    print()
    
    print("=== VERIFICATION ===")
    issues = []
    
    # Check if shear terms are non-zero
    if abs(B[3,1]) < 1e-10 and abs(B[3,7]) < 1e-10:
        issues.append("CRITICAL: Shear strain terms (gamma_xy) are ZERO - Timoshenko should include du_y/dx terms")
    
    if abs(B[3,5]) < 1e-10 and abs(B[3,11]) < 1e-10:
        issues.append("CRITICAL: Shear strain terms (gamma_xy) are ZERO - Timoshenko should include -theta_z terms")
    
    # Check if bending uses dtheta/dx (Timoshenko) or d^2u/dx^2 (Euler-Bernoulli)
    # For Timoshenko: kappa_z = dtheta_z/dx (first derivative of rotation)
    # For Euler-Bernoulli: kappa_z = d^2u_y/dx^2 (second derivative of displacement)
    # Current implementation uses d^2u_y/dx^2, which is Euler-Bernoulli, not Timoshenko
    
    if abs(B[2,1]) > 1e-10 or abs(B[2,7]) > 1e-10:
        issues.append("CRITICAL: Bending uses d^2u_y/dx^2 (Euler-Bernoulli) instead of dtheta_z/dx (Timoshenko)")
    
    if issues:
        print("FAILURES FOUND:")
        for issue in issues:
            print(f"  - {issue}")
        print()
        print("RECOMMENDATION: Timoshenko B-matrix needs to be corrected to include:")
        print("  1. Shear strain: gamma_xy = du_y/dx - theta_z")
        print("  2. Bending curvature: kappa_z = dtheta_z/dx (not d^2u_y/dx^2)")
        return False
    else:
        print("SUCCESS: B-matrix includes correct Timoshenko terms")
        return True

if __name__ == "__main__":
    verify_timoshenko_b_matrix()

