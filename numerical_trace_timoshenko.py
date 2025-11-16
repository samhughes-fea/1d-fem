"""
Numerical trace of Timoshenko stiffness matrix computation.

This will show the issues with missing shear terms in B-matrix and D-matrix.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.timoshenko.utilities.shape_functions import ShapeFunctionOperator
from pre_processing.element_library.timoshenko.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.timoshenko.utilities.D_matrix import MaterialStiffnessOperator

def trace_timoshenko_stiffness():
    """Trace Timoshenko stiffness computation to identify issues."""
    
    # Element properties
    L = 0.2  # m
    E = 2.1e11  # Pa
    I_z = 2.08769e-06  # m^4
    A = 0.00131  # m^2
    G = 8.1e10  # Pa
    kappa = 5.0 / 6.0  # Shear correction factor
    
    # Quadrature
    quadrature_order = 3
    xi_gauss, weights = np.polynomial.legendre.leggauss(quadrature_order)
    
    print("=" * 70)
    print("NUMERICAL TRACE: TIMOSHENKO STIFFNESS MATRIX")
    print("=" * 70)
    print(f"Element length L = {L:.3f} m")
    print(f"E = {E:.2e} Pa")
    print(f"G = {G:.2e} Pa")
    print(f"I_z = {I_z:.2e} m^4")
    print(f"A = {A:.6f} m^2")
    print(f"kappa = {kappa:.6f}")
    print()
    
    # Initialize operators
    shape_op = ShapeFunctionOperator(element_length=L)
    strain_op = StrainDisplacementOperator(element_length=L)
    material_op = MaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=3.23400e-07,
        moment_inertia_z=I_z,
        torsion_constant=2.60673e-08
    )
    
    # Get D-matrix
    D = material_op.assembly_form()
    
    print("=== D-MATRIX (Material Stiffness) ===")
    print(f"D[2,2] (EI_z, bending) = {D[2,2]:.6e} N*m^2")
    print(f"D[3,3] (shear xy) = {D[3,3]:.6e} N (should be kappa*G*A = {kappa*G*A:.6e})")
    print(f"D[4,4] (shear xz) = {D[4,4]:.6e} N (should be kappa*G*A = {kappa*G*A:.6e})")
    print()
    
    if abs(D[3,3]) < 1e-10:
        print("CRITICAL: D[3,3] is ZERO - Timoshenko needs shear stiffness!")
    if abs(D[4,4]) < 1e-10:
        print("CRITICAL: D[4,4] is ZERO - Timoshenko needs shear stiffness!")
    print()
    
    # Initialize stiffness matrix
    K_e = np.zeros((12, 12))
    
    # Trace for one Gauss point
    xi_g = 0.0
    w_g = weights[1]  # Middle point
    
    print("=" * 70)
    print(f"GAUSS POINT: xi = {xi_g:.6f}, weight = {w_g:.6f}")
    print("=" * 70)
    
    # Shape functions
    N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(np.array([xi_g]))
    N = N[0]
    dN_dxi = dN_dxi[0]
    d2N_dxi2 = d2N_dxi2[0]
    
    # B-matrix
    B = strain_op.physical_coordinate_form(dN_dxi[np.newaxis, :, :], d2N_dxi2[np.newaxis, :, :])[0]
    
    print("=== B-MATRIX (Strain-Displacement) ===")
    print("Shear strain terms (gamma_xy = du_y/dx - theta_z):")
    print(f"  B[3,1] (from u_y1, du_y/dx) = {B[3,1]:.6e}")
    print(f"  B[3,7] (from u_y2, du_y/dx) = {B[3,7]:.6e}")
    print(f"  B[3,5] (from theta_z1, -theta_z) = {B[3,5]:.6e}")
    print(f"  B[3,11] (from theta_z2, -theta_z) = {B[3,11]:.6e}")
    print()
    
    if abs(B[3,1]) < 1e-10 and abs(B[3,7]) < 1e-10:
        print("CRITICAL: B[3,1] and B[3,7] are ZERO - Timoshenko needs du_y/dx terms!")
    if abs(B[3,5]) < 1e-10 and abs(B[3,11]) < 1e-10:
        print("CRITICAL: B[3,5] and B[3,11] are ZERO - Timoshenko needs -theta_z terms!")
    print()
    
    # B^T @ D @ B
    BT_D = B.T @ D
    BT_D_B = BT_D @ B
    
    print("=== B^T @ D @ B (before integration) ===")
    print(f"(B^T @ D @ B)[1,1] (u_y1-u_y1) = {BT_D_B[1,1]:.6e}")
    print(f"(B^T @ D @ B)[1,5] (u_y1-theta_z1) = {BT_D_B[1,5]:.6e}")
    print()
    
    # Integration
    detJ = L / 2
    K_contrib = BT_D_B * w_g * detJ
    
    print("=== CONTRIBUTION TO K_e ===")
    print(f"K_contrib[1,1] = {K_contrib[1,1]:.6e}")
    print()
    
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("Issues identified:")
    print("  1. D-matrix missing shear stiffness (D[3,3] and D[4,4] are zero)")
    print("  2. B-matrix missing shear strain terms (B[3,:] terms are zero)")
    print("  3. Bending uses d^2u_y/dx^2 instead of dtheta_z/dx")
    print()
    print("These issues prevent Timoshenko elements from modeling shear deformation.")

if __name__ == "__main__":
    trace_timoshenko_stiffness()

