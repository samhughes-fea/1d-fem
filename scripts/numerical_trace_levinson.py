"""
Numerical trace of Levinson stiffness matrix computation.

This will show the issues with shape functions and potential missing shear terms.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.utilities.shape_functions import ShapeFunctionOperator
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.levinson.utilities.D_matrix import MaterialStiffnessOperator

def trace_levinson_stiffness():
    """Trace Levinson stiffness computation to identify issues."""
    
    # Element properties
    L = 0.2  # m
    E = 2.1e11  # Pa
    I_z = 2.08769e-06  # m^4
    A = 0.00131  # m^2
    G = 8.1e10  # Pa
    
    # Quadrature (may need higher order for quintic)
    quadrature_order = 4  # Higher order for quintic shape functions
    xi_gauss, weights = np.polynomial.legendre.leggauss(quadrature_order)
    
    print("=" * 70)
    print("NUMERICAL TRACE: LEVINSON STIFFNESS MATRIX")
    print("=" * 70)
    print(f"Element length L = {L:.3f} m")
    print(f"E = {E:.2e} Pa")
    print(f"G = {G:.2e} Pa")
    print(f"I_z = {I_z:.2e} m^4")
    print(f"A = {A:.6f} m^2")
    print(f"Quadrature order: {quadrature_order} (higher order for quintic)")
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
    print(f"D-matrix shape: {D.shape}")
    print(f"D[0,0] (EA, axial) = {D[0,0]:.6e} N")
    print(f"D[1,1] (EI_z, bending) = {D[1,1]:.6e} N*m^2")
    print()
    print("NOTE: D-matrix is 4x4 (axial, bending, torsion)")
    print("      Levinson should include GA for shear (no kappa factor)")
    print("      Need to verify if shear terms are included")
    print()
    
    # Check shape functions at nodes
    print("=== SHAPE FUNCTION VERIFICATION ===")
    xi_nodes = np.array([-1.0, 1.0])
    N_nodes, _, _ = shape_op.natural_coordinate_form(xi_nodes)
    
    print("At xi = -1 (node 1):")
    print(f"  N[1,1] (u_y1) = {N_nodes[0,1,1]:.6f} (should be 1.0)")
    print(f"  N[7,1] (u_y2) = {N_nodes[0,7,1]:.6f} (should be 0.0)")
    print()
    
    if abs(N_nodes[0,1,1] - 1.0) > 0.1:
        print("CRITICAL: Shape function N1_v(-1) is NOT 1.0!")
        print("          This will cause incorrect interpolation and stiffness.")
    print()
    
    # Trace for one Gauss point
    xi_g = 0.0
    w_g = weights[2]  # Middle point (for order 4)
    
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
    print(f"B-matrix shape: {B.shape}")
    print("Expected: (4, 12) for [eps_x, kappa_z, kappa_y, phi_x]")
    print()
    print("NOTE: B-matrix is 4x12 (axial, bending, torsion only)")
    print("      Levinson should include higher-order shear terms:")
    print("      gamma_xy = du_y/dx - theta_z + alpha*(d^2theta_z/dx^2)")
    print("      Need to verify if shear terms are included")
    print()
    
    # B^T @ D @ B
    BT_D = B.T @ D
    BT_D_B = BT_D @ B
    
    print("=== B^T @ D @ B (before integration) ===")
    print(f"(B^T @ D @ B)[1,1] (u_y1-u_y1) = {BT_D_B[1,1]:.6e}")
    print()
    
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("Issues identified:")
    print("  1. Shape functions do not satisfy interpolation properties")
    print("  2. B-matrix structure suggests shear terms may be missing")
    print("  3. D-matrix structure needs verification for shear terms")
    print()
    print("These issues need to be addressed for correct Levinson formulation.")

if __name__ == "__main__":
    trace_levinson_stiffness()

