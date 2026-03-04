"""
Verify Timoshenko D-matrix includes shear stiffness correctly.

CRITICAL ISSUE FOUND: Timoshenko D-matrix is currently identical to Euler-Bernoulli,
which means it does NOT include shear stiffness (κGA). This is incorrect for Timoshenko theory.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.timoshenko.utilities.D_matrix import MaterialStiffnessOperator

def verify_timoshenko_d_matrix():
    """Verify Timoshenko D-matrix includes shear stiffness."""
    
    # Job 0000 parameters
    E = 2.1e11  # Pa
    G = 8.1e10  # Pa
    A = 0.00131  # m^2
    I_y = 3.23400e-07  # m^4
    I_z = 2.08769e-06  # m^4
    J_t = 2.60673e-08  # m^4
    kappa = 5.0 / 6.0  # Shear correction factor for rectangular section
    
    material_op = MaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=I_y,
        moment_inertia_z=I_z,
        torsion_constant=J_t
    )
    
    D = material_op.assembly_form()
    
    print("=" * 70)
    print("TIMOSHENKO D-MATRIX VERIFICATION")
    print("=" * 70)
    print(f"E = {E:.2e} Pa")
    print(f"G = {G:.2e} Pa")
    print(f"A = {A:.6f} m^2")
    print(f"I_z = {I_z:.2e} m^4")
    print(f"kappa = {kappa:.6f} (shear correction factor)")
    print()
    
    print("=== D-MATRIX STRUCTURE ===")
    print("Strain: [eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]")
    print("Stress: [N, M_y, M_z, V_xy, V_xz, M_x]")
    print()
    
    print("=== KEY TERMS ===")
    print(f"D[0,0] (EA, axial) = {D[0,0]:.6e} N")
    print(f"D[1,1] (EI_y, bending about y) = {D[1,1]:.6e} N*m^2")
    print(f"D[2,2] (EI_z, bending about z) = {D[2,2]:.6e} N*m^2")
    print(f"D[3,3] (shear xy) = {D[3,3]:.6e} N")
    print(f"D[4,4] (shear xz) = {D[4,4]:.6e} N")
    print(f"D[5,5] (GJ_t, torsion) = {D[5,5]:.6e} N*m^2")
    print()
    
    print("=== VERIFICATION ===")
    issues = []
    
    # Expected values
    EA_expected = E * A
    EI_z_expected = E * I_z
    kappa_GA_expected = kappa * G * A
    GJ_t_expected = G * J_t
    
    # Check axial
    if abs(D[0,0] - EA_expected) > 1e-6:
        issues.append(f"EA incorrect: {D[0,0]:.6e} vs {EA_expected:.6e}")
    else:
        print(f"PASS: EA = {D[0,0]:.6e}")
    
    # Check bending
    if abs(D[2,2] - EI_z_expected) > 1e-6:
        issues.append(f"EI_z incorrect: {D[2,2]:.6e} vs {EI_z_expected:.6e}")
    else:
        print(f"PASS: EI_z = {D[2,2]:.6e}")
    
    # CRITICAL: Check shear stiffness
    if abs(D[3,3]) < 1e-10:
        issues.append("CRITICAL: D[3,3] (shear xy) is ZERO - Timoshenko should have kappa*G*A")
        print(f"FAIL: D[3,3] = {D[3,3]:.6e} (should be {kappa_GA_expected:.6e})")
    elif abs(D[3,3] - kappa_GA_expected) > 1e-6:
        issues.append(f"D[3,3] incorrect: {D[3,3]:.6e} vs {kappa_GA_expected:.6e}")
        print(f"FAIL: D[3,3] = {D[3,3]:.6e} (should be {kappa_GA_expected:.6e})")
    else:
        print(f"PASS: D[3,3] (kappa*G*A) = {D[3,3]:.6e}")
    
    if abs(D[4,4]) < 1e-10:
        issues.append("CRITICAL: D[4,4] (shear xz) is ZERO - Timoshenko should have kappa*G*A")
        print(f"FAIL: D[4,4] = {D[4,4]:.6e} (should be {kappa_GA_expected:.6e})")
    elif abs(D[4,4] - kappa_GA_expected) > 1e-6:
        issues.append(f"D[4,4] incorrect: {D[4,4]:.6e} vs {kappa_GA_expected:.6e}")
        print(f"FAIL: D[4,4] = {D[4,4]:.6e} (should be {kappa_GA_expected:.6e})")
    else:
        print(f"PASS: D[4,4] (kappa*G*A) = {D[4,4]:.6e}")
    
    # Check torsion
    if abs(D[5,5] - GJ_t_expected) > 1e-6:
        issues.append(f"GJ_t incorrect: {D[5,5]:.6e} vs {GJ_t_expected:.6e}")
    else:
        print(f"PASS: GJ_t = {D[5,5]:.6e}")
    
    print()
    if issues:
        print("FAILURES FOUND:")
        for issue in issues:
            print(f"  - {issue}")
        print()
        print("RECOMMENDATION: Timoshenko D-matrix needs to be corrected to include:")
        print(f"  D[3,3] = kappa*G*A = {kappa_GA_expected:.6e} N")
        print(f"  D[4,4] = kappa*G*A = {kappa_GA_expected:.6e} N")
        return False
    else:
        print("SUCCESS: D-matrix includes correct Timoshenko shear stiffness")
        return True

if __name__ == "__main__":
    verify_timoshenko_d_matrix()

