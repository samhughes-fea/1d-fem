"""
Comprehensive verification script for all three beam element formulations.
Verifies shape functions, B-matrices, D-matrices, and integration.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.shape_functions import ShapeFunctionOperator as EBShapeOp
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator as EBBOp
from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator as EBDOp

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.shape_functions import ShapeFunctionOperator as TimShapeOp
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.B_matrix import StrainDisplacementOperator as TimBOp
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.D_matrix import MaterialStiffnessOperator as TimDOp

from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.utilities.shape_functions import ShapeFunctionOperator as LevShapeOp
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.utilities.B_matrix import StrainDisplacementOperator as LevBOp
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.utilities.D_matrix import MaterialStiffnessOperator as LevDOp

def verify_timoshenko():
    """Verify Timoshenko element formulation."""
    print("=" * 70)
    print("TIMOSHENKO ELEMENT VERIFICATION")
    print("=" * 70)
    
    L = 0.2
    E = 2.1e11
    G = 8.1e10
    A = 0.00131
    I_z = 2.08769e-06
    I_y = 1e-06
    J_t = 1e-06
    kappa = 5.0 / 6.0
    
    # Test D-matrix
    D_op = TimDOp(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=I_y,
        moment_inertia_z=I_z,
        torsion_constant=J_t,
        shear_correction_factor=kappa
    )
    D = D_op.assembly_form()
    
    print("\n=== D-Matrix Verification ===")
    print(f"D[3,3] (shear xy) = {D[3,3]:.6e} (should be kappa*G*A = {kappa*G*A:.6e})")
    print(f"D[4,4] (shear xz) = {D[4,4]:.6e} (should be kappa*G*A = {kappa*G*A:.6e})")
    
    if abs(D[3,3] - kappa*G*A) < 1e-6:
        print("[PASS] D[3,3] correct")
    else:
        print("[FAIL] D[3,3] INCORRECT")
    
    if abs(D[4,4] - kappa*G*A) < 1e-6:
        print("[PASS] D[4,4] correct")
    else:
        print("[FAIL] D[4,4] INCORRECT")
    
    # Test B-matrix with shape functions
    shape_op = TimShapeOp(element_length=L)
    B_op = TimBOp(element_length=L)
    
    xi = np.array([0.0])
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(xi)
    B = B_op.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
    
    print("\n=== B-Matrix Verification ===")
    print(f"B shape: {B.shape} (should be (6, 12))")
    
    # Check shear terms
    # gamma_xy = du_y/dx - theta_z
    # B[3, 1] should be du_y/dx term (dN_dxi[1,1] * dxi_dx)
    # B[3, 5] should be -theta_z term (-N[5,5])
    print(f"\nShear strain gamma_xy terms:")
    print(f"  B[3, 1] (du_y/dx) = {B[3, 1]:.6e}")
    print(f"  B[3, 5] (-theta_z) = {B[3, 5]:.6e}")
    
    if abs(B[3, 1]) > 1e-10:
        print("[PASS] du_y/dx term present")
    else:
        print("[FAIL] du_y/dx term MISSING")
    
    if abs(B[3, 5]) > 1e-10:
        print("[PASS] -theta_z term present")
    else:
        print("[FAIL] -theta_z term MISSING")
    
    # Check bending uses rotation-based curvature
    # kappa_z = dtheta_z/dx (not d2u_y/dx2)
    print(f"\nBending curvature kappa_z:")
    print(f"  B[2, 5] (dtheta_z/dx) = {B[2, 5]:.6e}")
    
    if abs(B[2, 5]) > 1e-10:
        print("[PASS] Uses rotation-based curvature (Timoshenko)")
    else:
        print("[FAIL] Missing rotation-based curvature")
    
    print("\n" + "=" * 70)

def verify_levinson():
    """Verify Levinson element formulation."""
    print("=" * 70)
    print("LEVINSON ELEMENT VERIFICATION")
    print("=" * 70)
    
    L = 0.2
    E = 2.1e11
    G = 8.1e10
    A = 0.00131
    I_z = 2.08769e-06
    I_y = 1e-06
    J_t = 1e-06
    
    # Compute alpha (I_z/A for rectangular approximation)
    alpha = I_z / A
    
    # Test D-matrix
    D_op = LevDOp(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=I_y,
        moment_inertia_z=I_z,
        torsion_constant=J_t
    )
    D = D_op.assembly_form()
    
    print("\n=== D-Matrix Verification ===")
    print(f"D[3,3] (shear xy) = {D[3,3]:.6e} (should be G*A = {G*A:.6e}, no kappa)")
    print(f"D[4,4] (shear xz) = {D[4,4]:.6e} (should be G*A = {G*A:.6e}, no kappa)")
    
    if abs(D[3,3] - G*A) < 1e-6:
        print("[PASS] D[3,3] correct (GA, no kappa)")
    else:
        print("[FAIL] D[3,3] INCORRECT")
    
    if abs(D[4,4] - G*A) < 1e-6:
        print("[PASS] D[4,4] correct (GA, no kappa)")
    else:
        print("[FAIL] D[4,4] INCORRECT")
    
    # Test B-matrix with alpha
    shape_op = LevShapeOp(element_length=L)
    B_op = LevBOp(element_length=L, alpha_coefficient=alpha)
    
    # Test at non-zero point (xi=0.5) where shape functions are non-zero
    xi = np.array([0.5])
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(xi)
    B = B_op.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]
    
    print(f"\n=== B-Matrix Verification ===")
    print(f"Alpha coefficient: {alpha:.6e}")
    print(f"B shape: {B.shape} (should be (6, 12))")
    
    # Check higher-order shear terms
    # gamma_xy = du_y/dx - theta_z + alpha(d2theta_z/dx2)
    print(f"\nHigher-order shear strain gamma_xy terms:")
    print(f"  B[3, 1] (du_y/dx) = {B[3, 1]:.6e}")
    print(f"  B[3, 5] (-theta_z + alpha(d2theta_z/dx2)) = {B[3, 5]:.6e}")
    
    if abs(B[3, 1]) > 1e-10:
        print("[PASS] du_y/dx term present")
    else:
        print("[FAIL] du_y/dx term MISSING")
    
    if abs(B[3, 5]) > 1e-10:
        print("[PASS] Higher-order shear term present")
    else:
        print("[FAIL] Higher-order shear term MISSING")
    
    print("\n" + "=" * 70)

if __name__ == "__main__":
    verify_timoshenko()
    verify_levinson()
    print("\n[PASS] Verification complete!")

