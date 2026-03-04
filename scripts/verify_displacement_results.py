"""
Verify displacement results against analytical solutions.

Compares FEM computed displacements with analytical beam theory solutions
for cantilever beam with point load.
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

def analytical_cantilever_tip_deflection(P, E, I_z, L):
    """
    Analytical tip deflection for cantilever beam with point load at tip.
    
    δ = PL³/(3EI)
    
    Parameters
    ----------
    P : float
        Point load at tip (N, negative for downward)
    E : float
        Young's modulus (Pa)
    I_z : float
        Moment of inertia about z-axis (m⁴)
    L : float
        Beam length (m)
    
    Returns
    -------
    delta : float
        Tip deflection (m, negative for downward)
    """
    EI = E * I_z
    delta = (P * L**3) / (3 * EI)
    return delta

def analytical_cantilever_rotation(P, E, I_z, L):
    """
    Analytical tip rotation for cantilever beam with point load at tip.
    
    θ = PL²/(2EI)
    
    Parameters
    ----------
    P : float
        Point load at tip (N)
    E : float
        Young's modulus (Pa)
    I_z : float
        Moment of inertia about z-axis (m⁴)
    L : float
        Beam length (m)
    
    Returns
    -------
    theta : float
        Tip rotation (rad)
    """
    EI = E * I_z
    theta = (P * L**2) / (2 * EI)
    return theta

def verify_job_results(job_name, results_dir):
    """Verify displacement results for a single job."""
    print("=" * 70)
    print(f"VERIFICATION: {job_name}")
    print("=" * 70)
    
    # Read input files
    job_dir = PROJECT_ROOT / "jobs" / job_name
    
    # Material properties
    material_file = job_dir / "material.txt"
    with open(material_file, 'r') as f:
        lines = f.readlines()
    # Extract E (first material, second column after header)
    E = float(lines[4].split()[1])  # Line 4 (0-indexed) has first material data
    
    # Section properties
    section_file = job_dir / "section.txt"
    with open(section_file, 'r') as f:
        lines = f.readlines()
    # Extract I_z (first section, 4th column after header)
    I_z = float(lines[4].split()[4])  # I_z is 5th column (index 4)
    
    # Grid (to get beam length)
    grid_file = job_dir / "grid.txt"
    # Read grid file - format: [node_id] [x] [y] [z]
    grid_data = pd.read_csv(grid_file, comment='#', skipinitialspace=True, sep=r'\s+', 
                            names=['node_id', 'x', 'y', 'z'], skiprows=3)
    x_coords = grid_data['x'].values
    L = x_coords[-1] - x_coords[0]  # Total beam length
    
    # Point load
    point_load_file = job_dir / "point_load.txt"
    with open(point_load_file, 'r') as f:
        lines = f.readlines()
    # Extract load (line 5, F_y is 5th column)
    load_line = lines[4].split()
    P = float(load_line[4])  # F_y is 5th column (index 4)
    
    print(f"\nInput Parameters:")
    print(f"  E = {E:.2e} Pa")
    print(f"  I_z = {I_z:.2e} m^4")
    print(f"  EI_z = {E * I_z:.2e} N*m^2")
    print(f"  L = {L:.3f} m")
    print(f"  P = {P:.1f} N (point load at tip)")
    
    # Read computed results
    U_file = Path(results_dir) / "primary_results" / "global" / "U_global.csv"
    if not U_file.exists():
        print(f"\n[ERROR] Results file not found: {U_file}")
        return False
    
    U_df = pd.read_csv(U_file)
    
    # Find tip node (last node, node 10 for 11 nodes)
    # DOF mapping: node i has DOFs [6*i, 6*i+1, 6*i+2, 6*i+3, 6*i+4, 6*i+5]
    # For node 10: DOFs 60-65
    # u_y is DOF 61 (6*10 + 1)
    # theta_z is DOF 65 (6*10 + 5)
    
    n_nodes = len(x_coords)
    tip_node = n_nodes - 1
    tip_u_y_dof = 6 * tip_node + 1
    tip_theta_z_dof = 6 * tip_node + 5
    
    # Get computed values
    computed_u_y = U_df[U_df['Global DOF'] == tip_u_y_dof]['Value'].values[0]
    computed_theta_z = U_df[U_df['Global DOF'] == tip_theta_z_dof]['Value'].values[0]
    
    # Analytical solutions
    analytical_u_y = analytical_cantilever_tip_deflection(P, E, I_z, L)
    analytical_theta_z = analytical_cantilever_rotation(P, E, I_z, L)
    
    print(f"\nTip Node (Node {tip_node}, x = {x_coords[tip_node]:.3f} m):")
    print(f"  DOF {tip_u_y_dof} (u_y):")
    print(f"    Computed:   {computed_u_y*1000:.6f} mm")
    print(f"    Analytical: {analytical_u_y*1000:.6f} mm")
    print(f"    Error:      {abs(computed_u_y - analytical_u_y)*1000:.6f} mm")
    print(f"    Relative:   {abs((computed_u_y - analytical_u_y) / analytical_u_y) * 100:.4f}%")
    
    print(f"\n  DOF {tip_theta_z_dof} (theta_z):")
    print(f"    Computed:   {computed_theta_z*180/np.pi:.6f} deg")
    print(f"    Analytical: {analytical_theta_z*180/np.pi:.6f} deg")
    print(f"    Error:      {abs(computed_theta_z - analytical_theta_z)*180/np.pi:.6f} deg")
    print(f"    Relative:   {abs((computed_theta_z - analytical_theta_z) / analytical_theta_z) * 100:.4f}%")
    
    # Check boundary conditions at fixed end (node 0)
    fixed_u_y_dof = 1  # Node 0, u_y
    fixed_theta_z_dof = 5  # Node 0, theta_z
    fixed_u_y = U_df[U_df['Global DOF'] == fixed_u_y_dof]['Value'].values[0]
    fixed_theta_z = U_df[U_df['Global DOF'] == fixed_theta_z_dof]['Value'].values[0]
    
    print(f"\nBoundary Conditions (Node 0):")
    print(f"  DOF {fixed_u_y_dof} (u_y): {fixed_u_y:.2e} m (should be 0.0)")
    print(f"  DOF {fixed_theta_z_dof} (theta_z): {fixed_theta_z:.2e} rad (should be 0.0)")
    
    # Verification
    tol = 1e-6  # 1 micron tolerance
    u_y_ok = abs(computed_u_y - analytical_u_y) < tol
    theta_z_ok = abs(computed_theta_z - analytical_theta_z) < tol
    bc_ok = abs(fixed_u_y) < 1e-12 and abs(fixed_theta_z) < 1e-12
    
    print(f"\n{'='*70}")
    if u_y_ok and theta_z_ok and bc_ok:
        print("[PASS] All verifications passed!")
        print(f"  Tip deflection matches analytical within {tol*1000:.3f} mm")
        print(f"  Tip rotation matches analytical")
        print(f"  Boundary conditions correctly applied")
        return True
    else:
        print("[FAIL] Some verifications failed:")
        if not u_y_ok:
            print(f"  Tip deflection error: {abs(computed_u_y - analytical_u_y)*1000:.6f} mm > {tol*1000:.3f} mm")
        if not theta_z_ok:
            print(f"  Tip rotation error: {abs(computed_theta_z - analytical_theta_z)*180/np.pi:.6f} deg")
        if not bc_ok:
            print(f"  Boundary conditions not properly fixed")
        return False

if __name__ == "__main__":
    # Find latest results directory for job_0000
    results_base = PROJECT_ROOT / "post_processing" / "results"
    job_0000_dirs = list(results_base.glob("job_0000_*"))
    
    if not job_0000_dirs:
        print("[ERROR] No results found for job_0000")
        sys.exit(1)
    
    # Use most recent
    latest_dir = max(job_0000_dirs, key=lambda p: p.stat().st_mtime)
    
    print(f"Using results directory: {latest_dir.name}")
    print()
    
    success = verify_job_results("job_0000", latest_dir)
    
    if success:
        print("\n" + "="*70)
        print("VERIFICATION SUMMARY: Results are CORRECT")
        print("="*70)
        sys.exit(0)
    else:
        print("\n" + "="*70)
        print("VERIFICATION SUMMARY: Results have ERRORS")
        print("="*70)
        sys.exit(1)

