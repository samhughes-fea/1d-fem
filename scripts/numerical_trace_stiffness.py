"""
Numerical trace of stiffness matrix computation for one element.

Traces step-by-step the computation of K_e to identify where the scaling error occurs.
"""

import numpy as np
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.euler_bernoulli.utilities.shape_functions import ShapeFunctionOperator
from pre_processing.element_library.linear.euler_bernoulli.utilities.B_matrix import StrainDisplacementOperator
from pre_processing.element_library.linear.euler_bernoulli.utilities.D_matrix import MaterialStiffnessOperator

def trace_stiffness_computation():
    """Trace stiffness computation for job_0000 element 0."""
    
    # Element properties
    L = 0.2  # m
    E = 2.1e11  # Pa
    I_z = 2.08769e-06  # m^4
    A = 0.00131  # m^2
    I_y = 3.23400e-07  # m^4
    G = 8.1e10  # Pa
    J_t = 2.60673e-08  # m^4
    
    # Quadrature order
    quadrature_order = 3
    xi_gauss, weights = np.polynomial.legendre.leggauss(quadrature_order)
    
    print("=" * 70)
    print("NUMERICAL TRACE: Stiffness Matrix Computation")
    print("=" * 70)
    print(f"Element length L = {L:.3f} m")
    print(f"E = {E:.2e} Pa")
    print(f"I_z = {I_z:.2e} m^4")
    print(f"EI_z = {E * I_z:.2e} N*m^2")
    print(f"Quadrature order: {quadrature_order}")
    print(f"Gauss points: {xi_gauss}")
    print(f"Gauss weights: {weights}")
    print()
    
    # Initialize operators
    shape_op = ShapeFunctionOperator(element_length=L)
    strain_op = StrainDisplacementOperator(element_length=L)
    material_op = MaterialStiffnessOperator(
        youngs_modulus=E,
        shear_modulus=G,
        cross_section_area=A,
        moment_inertia_y=I_y,
        moment_inertia_z=I_z,
        torsion_constant=J_t
    )
    
    # Get D-matrix
    D = material_op.assembly_form()
    print("=== D-Matrix (Material Stiffness) ===")
    print(f"D[2,2] (EI_z for bending about z) = {D[2,2]:.6e} N*m^2")
    print()
    
    # Coordinate transformation factors
    detJ = L / 2
    dxi_dx = 2 / L
    d2xi_dx2 = 4 / (L**2)
    
    print("=== Coordinate Transformation Factors ===")
    print(f"detJ = L/2 = {detJ:.6f} m")
    print(f"dxi_dx = 2/L = {dxi_dx:.6f}")
    print(f"d2xi_dx2 = 4/L^2 = {d2xi_dx2:.6e} (1/m^2)")
    print(f"Note: d2xi_dx2 is actually (dxi/dx)^2, not d^2xi/dx^2")
    print()
    
    # Initialize stiffness matrix
    K_e = np.zeros((12, 12))
    
    # Trace computation for each Gauss point
    for g, (xi_g, w_g) in enumerate(zip(xi_gauss, weights)):
        print(f"{'=' * 70}")
        print(f"GAUSS POINT {g+1}: xi = {xi_g:.6f}, weight = {w_g:.6f}")
        print(f"{'=' * 70}")
        
        # Step 1: Shape functions
        N, dN_dxi, d2N_dxi2 = shape_op.natural_coordinate_form(np.array([xi_g]))
        N = N[0]  # Remove first dimension
        dN_dxi = dN_dxi[0]
        d2N_dxi2 = d2N_dxi2[0]
        
        print(f"\nStep 1: Shape Function Derivatives (at xi={xi_g:.6f})")
        print(f"  d2N_dxi2[1,1] (u_y1, component 1) = {d2N_dxi2[1,1]:.6f}")
        print(f"  d2N_dxi2[7,1] (u_y2, component 1) = {d2N_dxi2[7,1]:.6f}")
        print(f"  d2N_dxi2[5,5] (theta_z1, component 5) = {d2N_dxi2[5,5]:.6f}")
        print(f"  d2N_dxi2[11,5] (theta_z2, component 5) = {d2N_dxi2[11,5]:.6f}")
        
        # Step 2: B-matrix
        B = strain_op.physical_coordinate_form(dN_dxi[np.newaxis, :, :], d2N_dxi2[np.newaxis, :, :])[0]
        
        print(f"\nStep 2: B-Matrix (physical coordinates)")
        print(f"  B[2,1] (kappa_z from u_y1) = d2N_dxi2[1,1] * d2xi_dx2 = {d2N_dxi2[1,1]:.6f} * {d2xi_dx2:.6e} = {B[2,1]:.6e}")
        print(f"  B[2,7] (kappa_z from u_y2) = d2N_dxi2[7,1] * d2xi_dx2 = {d2N_dxi2[7,1]:.6f} * {d2xi_dx2:.6e} = {B[2,7]:.6e}")
        print(f"  B[2,5] (kappa_z from theta_z1) = d2N_dxi2[5,5] * d2xi_dx2 = {d2N_dxi2[5,5]:.6f} * {d2xi_dx2:.6e} = {B[2,5]:.6e}")
        print(f"  B[2,11] (kappa_z from theta_z2) = d2N_dxi2[11,5] * d2xi_dx2 = {d2N_dxi2[11,5]:.6f} * {d2xi_dx2:.6e} = {B[2,11]:.6e}")
        
        # Step 3: B^T @ D @ B
        BT_D = B.T @ D  # [12, 6]
        BT_D_B = BT_D @ B  # [12, 12]
        
        print(f"\nStep 3: B^T @ D @ B (before integration factors)")
        print(f"  (B^T @ D @ B)[1,1] = {BT_D_B[1,1]:.6e}")
        print(f"  (B^T @ D @ B)[1,5] = {BT_D_B[1,5]:.6e}")
        print(f"  (B^T @ D @ B)[5,5] = {BT_D_B[5,5]:.6e}")
        
        # Step 4: Multiply by weight and detJ
        K_contrib = BT_D_B * w_g * detJ
        
        print(f"\nStep 4: Contribution to K_e (multiply by weight and detJ)")
        print(f"  weight = {w_g:.6f}")
        print(f"  detJ = {detJ:.6f} m")
        print(f"  K_contrib[1,1] = {K_contrib[1,1]:.6e}")
        print(f"  K_contrib[1,5] = {K_contrib[1,5]:.6e}")
        print(f"  K_contrib[5,5] = {K_contrib[5,5]:.6e}")
        
        K_e += K_contrib
        print()
    
    print("=" * 70)
    print("FINAL STIFFNESS MATRIX (sum over all Gauss points)")
    print("=" * 70)
    print(f"K_e[1,1] (u_y1-u_y1) = {K_e[1,1]:.6e}")
    print(f"K_e[1,5] (u_y1-theta_z1) = {K_e[1,5]:.6e}")
    print(f"K_e[5,5] (theta_z1-theta_z1) = {K_e[5,5]:.6e}")
    print(f"K_e[7,7] (u_y2-u_y2) = {K_e[7,7]:.6e}")
    print(f"K_e[11,11] (theta_z2-theta_z2) = {K_e[11,11]:.6e}")
    print()
    
    # Compare with analytical
    EI = E * I_z
    K11_analytical = 12 * EI / (L**3)
    K15_analytical = 6 * EI / (L**2)
    K55_analytical = 4 * EI / L
    
    print("=" * 70)
    print("COMPARISON WITH ANALYTICAL")
    print("=" * 70)
    print(f"Analytical K[1,1] = 12*EI/L^3 = {K11_analytical:.6e}")
    print(f"Computed   K[1,1] = {K_e[1,1]:.6e}")
    print(f"Ratio: {K_e[1,1] / K11_analytical:.2f}x")
    print()
    print(f"Analytical K[1,5] = 6*EI/L^2 = {K15_analytical:.6e}")
    print(f"Computed   K[1,5] = {K_e[1,5]:.6e}")
    print(f"Ratio: {K_e[1,5] / K15_analytical:.2f}x")
    print()
    print(f"Analytical K[5,5] = 4*EI/L = {K55_analytical:.6e}")
    print(f"Computed   K[5,5] = {K_e[5,5]:.6e}")
    print(f"Ratio: {K_e[5,5] / K55_analytical:.2f}x")
    print()
    
    # Check scaling factors
    print("=" * 70)
    print("SCALING ANALYSIS")
    print("=" * 70)
    print("B-matrix scaling: d2xi_dx2 = 4/L^2")
    print("B^T D B scaling: (4/L^2)^2 = 16/L^4")
    print("Integration: detJ = L/2")
    print("Total scaling: (16/L^4) * (L/2) = 8/L^3")
    print("Expected scaling: 1/L^3")
    print("Extra factor: 8x")
    print()
    print("But we see 112x error, not 8x. Investigating further...")

if __name__ == "__main__":
    trace_stiffness_computation()

