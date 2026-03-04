"""
Verify Timoshenko element stiffness matrix integration.

Checks that K^e = ∫ B^T D B |J| dξ is computed correctly,
focusing on shear terms.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.timoshenko.utilities.shape_functions import ShapeFunctionOperator
from pre_processing.element_library.linear.timoshenko.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.linear.timoshenko.utilities.D_matrix import MaterialStiffnessOperator

def verify_timoshenko_stiffness_integration():
    """Verify Timoshenko stiffness matrix integration."""
    print("=" * 70)
    print("VERIFICATION: Timoshenko Stiffness Matrix Integration")
    print("=" * 70)
    
    # Element properties (from job_0001)
    L = 0.2  # m
    E = 2.1e11  # Pa
    G = 8.1e10  # Pa
    I_z = 2.08769e-06  # m^4
    A = 0.00131  # m^2
    kappa = 5.0 / 6.0
    
    print(f"\nElement Properties:")
    print(f"  L = {L:.3f} m")
    print(f"  E = {E:.2e} Pa")
    print(f"  G = {G:.2e} Pa")
    print(f"  I_z = {I_z:.2e} m^4")
    print(f"  A = {A:.6f} m^2")
    print(f"  kappa = {kappa:.4f}")
    
    # Initialize operators
    shape_op = ShapeFunctionOperator(element_length=L)
    strain_op = StrainDisplacementOperator(element_length=L)
    material_op = MaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=3.23400e-07,
        moment_inertia_z=I_z,
        torsion_constant=2.60673e-08,
        shear_correction_factor=kappa
    )
    
    # Get D-matrix
    D = material_op.assembly_form()
    print(f"\nD-Matrix:")
    print(f"  D[2,2] (EI_z, bending): {D[2,2]:.6e} N*m^2")
    print(f"  D[3,3] (kappa*G*A, shear): {D[3,3]:.6e} N")
    print(f"  D[4,4] (kappa*G*A, shear): {D[4,4]:.6e} N")
    
    # Expected values
    expected_EI_z = E * I_z
    expected_kappa_GA = kappa * G * A
    print(f"\nExpected D-Matrix values:")
    print(f"  D[2,2] should be: {expected_EI_z:.6e} N*m^2")
    print(f"  D[3,3] should be: {expected_kappa_GA:.6e} N")
    
    # Check D-matrix
    d_ok = True
    if abs(D[2,2] - expected_EI_z) > 1e-6 * expected_EI_z:
        print(f"  [ERROR] D[2,2] mismatch!")
        d_ok = False
    if abs(D[3,3] - expected_kappa_GA) > 1e-6 * expected_kappa_GA:
        print(f"  [ERROR] D[3,3] mismatch!")
        d_ok = False
    if D[3,3] < 1e-6:
        print(f"  [ERROR] D[3,3] is essentially zero - shear stiffness missing!")
        d_ok = False
    
    # Compute stiffness matrix using 3-point Gauss quadrature
    quadrature_order = 3
    xi_gauss, weights = np.polynomial.legendre.leggauss(quadrature_order)
    detJ = L / 2.0
    
    print(f"\nGauss Quadrature (order {quadrature_order}):")
    print(f"  Points: {xi_gauss}")
    print(f"  Weights: {weights}")
    print(f"  detJ = {detJ:.6f} m")
    
    Ke = np.zeros((12, 12))
    
    print(f"\nIntegrating K^e = integral(B^T D B |J| dxi):")
    for g, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
        # Get shape functions and B-matrix
        N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(np.array([xi_g]))
        B = strain_op.physical_coordinate_form(dN_dxi, d2N_dxi2, N)[0]  # Shape (6, 12)
        
        # Contribution: B^T D B * w * detJ
        Ke_contrib = B.T @ D @ B * w_g * detJ
        Ke += Ke_contrib
        
        # Check shear terms in contribution
        print(f"\n  Gauss Point {g+1} (xi = {xi_g:.6f}):")
        print(f"    B[3, :] norm (shear row): {np.linalg.norm(B[3, :]):.6e}")
        print(f"    Ke_contrib[1, 1] (u_y1-u_y1): {Ke_contrib[1, 1]:.6e}")
        print(f"    Ke_contrib[1, 5] (u_y1-theta_z1): {Ke_contrib[1, 5]:.6e}")
        print(f"    Ke_contrib[5, 5] (theta_z1-theta_z1): {Ke_contrib[5, 5]:.6e}")
        print(f"    Ke_contrib[7, 7] (u_y2-u_y2): {Ke_contrib[7, 7]:.6e}")
        print(f"    Ke_contrib[11, 11] (theta_z2-theta_z2): {Ke_contrib[11, 11]:.6e}")
    
    print(f"\nFinal Stiffness Matrix (selected terms):")
    print(f"  Ke[1, 1] (u_y1-u_y1): {Ke[1, 1]:.6e}")
    print(f"  Ke[1, 5] (u_y1-theta_z1): {Ke[1, 5]:.6e}")
    print(f"  Ke[5, 5] (theta_z1-theta_z1): {Ke[5, 5]:.6e}")
    print(f"  Ke[7, 7] (u_y2-u_y2): {Ke[7, 7]:.6e}")
    print(f"  Ke[11, 11] (theta_z2-theta_z2): {Ke[11, 11]:.6e}")
    
    # Check if stiffness is reasonable
    # For a cantilever, the tip deflection should be δ = PL^3/(3EI) + PL/(κGA)
    # The stiffness should relate to this
    
    print(f"\n" + "=" * 70)
    if d_ok:
        print("[PASS] D-matrix values are correct")
    else:
        print("[FAIL] D-matrix values are incorrect")
    print("=" * 70)
    
    return d_ok, Ke

if __name__ == "__main__":
    verify_timoshenko_stiffness_integration()

