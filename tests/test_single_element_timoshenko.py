"""
Single-element test for Timoshenko beam element.

Tests isolated element to determine if issue is in element formulation
or global assembly/solver.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
    LinearTimoshenkoBeamElement3D,
)
from scripts.analytical_timoshenko_benchmark import analytical_timoshenko_tip_deflection

def test_single_element_timoshenko():
    """Test single Timoshenko element with direct solve."""
    print("=" * 70)
    print("SINGLE ELEMENT TEST: Timoshenko")
    print("=" * 70)
    
    # Element properties (same as job_0001)
    L = 2.0  # m (full beam length in single element)
    E = 2.1e11  # Pa
    G = 8.1e10  # Pa
    I_z = 2.08769e-06  # m^4
    A = 0.00131  # m^2
    kappa = 5.0 / 6.0
    
    # Load
    P = -500.0  # N (point load at tip)
    
    print(f"\nElement Properties:")
    print(f"  L = {L:.3f} m")
    print(f"  E = {E:.2e} Pa")
    print(f"  G = {G:.2e} Pa")
    print(f"  I_z = {I_z:.2e} m^4")
    print(f"  A = {A:.6f} m^2")
    print(f"  kappa = {kappa:.4f}")
    print(f"  P = {P:.1f} N")
    
    # Create element dictionaries (minimal setup)
    element_dictionary = {
        "ids": np.array([0]),
        "connectivity": np.array([[0, 1]]),
        "types": np.array(["TimoshenkoBeamElement3D"]),
        "integration_orders": {
            "axial": np.array([3]),
            "bending_y": np.array([3]),
            "bending_z": np.array([3]),
            "shear_y": np.array([2]),
            "shear_z": np.array([2]),
            "torsion": np.array([3]),
            "load": np.array([2]),
        }
    }
    
    grid_dictionary = {
        "coordinates": np.array([[0.0, 0.0, 0.0], [L, 0.0, 0.0]])
    }
    
    material_dictionary = {
        "E": np.array([E]),
        "G": np.array([G]),
        "nu": np.array([0.3]),
        "rho": np.array([7850.0])
    }
    
    section_dictionary = {
        "A": np.array([A]),
        "I_x": np.array([0.0]),
        "I_y": np.array([3.23400e-07]),
        "I_z": np.array([I_z]),
        "J_t": np.array([2.60673e-08])
    }
    
    # Point load at tip (node 1, x = L)
    point_load_array = np.array([[L, 0.0, 0.0, 0.0, P, 0.0, 0.0, 0.0, 0.0]])
    distributed_load_array = np.empty((0, 9))
    
    # Create temporary results directory with required subdirectories
    import tempfile
    import os
    temp_dir = tempfile.mkdtemp()
    job_results_dir = os.path.join(temp_dir, "test_results")
    os.makedirs(job_results_dir, exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_stiffness_matrices"), exist_ok=True)
    os.makedirs(os.path.join(job_results_dir, "element_force_vectors"), exist_ok=True)
    
    try:
        # Create element
        element = LinearTimoshenkoBeamElement3D(
            element_id=0,
            element_dictionary=element_dictionary,
            grid_dictionary=grid_dictionary,
            material_dictionary=material_dictionary,
            section_dictionary=section_dictionary,
            point_load_array=point_load_array,
            distributed_load_array=distributed_load_array,
            job_results_dir=job_results_dir
        )
        
        # Compute stiffness and force
        element_obj = element.element_stiffness_matrix()
        Ke = element_obj.K_e
        
        force_obj = element.element_force_vector()
        Fe = force_obj.F_e
        
        print(f"\nElement Stiffness Matrix (selected terms):")
        print(f"  Ke[1, 1] (u_y1-u_y1): {Ke[1, 1]:.6e}")
        print(f"  Ke[1, 5] (u_y1-theta_z1): {Ke[1, 5]:.6e}")
        print(f"  Ke[5, 5] (theta_z1-theta_z1): {Ke[5, 5]:.6e}")
        print(f"  Ke[7, 7] (u_y2-u_y2): {Ke[7, 7]:.6e}")
        print(f"  Ke[11, 11] (theta_z2-theta_z2): {Ke[11, 11]:.6e}")
        
        print(f"\nElement Force Vector:")
        print(f"  Fe[1] (node 1, u_y): {Fe[1]:.6f} N")
        print(f"  Fe[7] (node 2, u_y): {Fe[7]:.6f} N")
        
        # Apply boundary conditions: node 0 fixed (all 6 DOFs = 0)
        # DOFs: [0:u_x, 1:u_y, 2:u_z, 3:theta_x, 4:theta_y, 5:theta_z] for each node
        fixed_dofs = [0, 1, 2, 3, 4, 5]  # Node 0 DOFs
        free_dofs = [6, 7, 8, 9, 10, 11]  # Node 1 DOFs
        
        # Extract free-free submatrix
        Ke_ff = Ke[np.ix_(free_dofs, free_dofs)]
        Fe_f = Fe[free_dofs]
        
        # Solve: u_f = Ke_ff^(-1) * Fe_f
        u_f = np.linalg.solve(Ke_ff, Fe_f)
        
        print(f"\nFree DOF Displacements:")
        print(f"  u[6] (node 1, u_x): {u_f[0]:.6e} m")
        print(f"  u[7] (node 1, u_y): {u_f[1]:.6e} m = {u_f[1]*1000:.6f} mm")
        print(f"  u[8] (node 1, u_z): {u_f[2]:.6e} m")
        print(f"  u[9] (node 1, theta_x): {u_f[3]:.6e} rad")
        print(f"  u[10] (node 1, theta_y): {u_f[4]:.6e} rad")
        print(f"  u[11] (node 1, theta_z): {u_f[5]:.6e} rad = {u_f[5]*180/np.pi:.6f} deg")
        
        # Analytical solution
        analytical_u_y, delta_bending, delta_shear = analytical_timoshenko_tip_deflection(
            P, E, I_z, L, G, A, kappa
        )
        analytical_theta_z = (P * L**2) / (2 * E * I_z)  # Same as EB for rotation
        
        print(f"\nAnalytical Solution:")
        print(f"  u_y (total): {analytical_u_y*1000:.6f} mm")
        print(f"    Bending: {delta_bending*1000:.6f} mm")
        print(f"    Shear: {delta_shear*1000:.6f} mm")
        print(f"  theta_z: {analytical_theta_z*180/np.pi:.6f} deg")
        
        # Compare
        u_y_error = abs(u_f[1] - analytical_u_y) * 1000  # mm
        theta_z_error = abs(u_f[5] - analytical_theta_z) * 180 / np.pi  # deg
        
        print(f"\nComparison:")
        print(f"  u_y error: {u_y_error:.6f} mm")
        print(f"  theta_z error: {theta_z_error:.6f} deg")
        
        tol = 0.01  # 0.01 mm tolerance
        u_y_match = u_y_error < tol
        theta_z_match = theta_z_error < 0.01  # 0.01 deg tolerance
        
        print(f"\n" + "=" * 70)
        if u_y_match and theta_z_match:
            print("[PASS] Single element test passed")
            print("  Issue is likely in global assembly or solver")
        else:
            print("[FAIL] Single element test failed")
            print("  Issue is in element formulation")
            if not u_y_match:
                print(f"  u_y mismatch: {u_y_error:.6f} mm > {tol:.3f} mm")
            if not theta_z_match:
                print(f"  theta_z mismatch: {theta_z_error:.6f} deg")
        print("=" * 70)
        
        return u_y_match and theta_z_match
        
    finally:
        # Cleanup
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    test_single_element_timoshenko()

