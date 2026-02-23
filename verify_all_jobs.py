"""
Verify displacement results for all three jobs against analytical solutions.

Compares FEM computed displacements with analytical beam theory solutions
for all three beam element formulations: Euler-Bernoulli, Timoshenko, and Levinson.
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from analytical_benchmark import analytical_euler_bernoulli_stiffness
from analytical_timoshenko_benchmark import analytical_timoshenko_tip_deflection
from analytical_levinson_benchmark import analytical_levinson_tip_deflection

def analytical_cantilever_tip_deflection_eb(P, E, I_z, L):
    """Euler-Bernoulli tip deflection: δ = PL³/(3EI)"""
    EI = E * I_z
    delta = (P * L**3) / (3 * EI)
    return delta

def analytical_cantilever_rotation(P, E, I_z, L):
    """Tip rotation: θ = PL²/(2EI)"""
    EI = E * I_z
    theta = (P * L**2) / (2 * EI)
    return theta

def verify_job_results(job_name, results_dir, element_type):
    """Verify displacement results for a single job."""
    print("=" * 70)
    print(f"VERIFICATION: {job_name} ({element_type})")
    print("=" * 70)
    
    # Read input files
    job_dir = PROJECT_ROOT / "jobs" / job_name
    
    # Material properties
    material_file = job_dir / "material.txt"
    with open(material_file, 'r') as f:
        lines = f.readlines()
    E = float(lines[4].split()[1])
    G = float(lines[4].split()[2])
    
    # Section properties
    section_file = job_dir / "section.txt"
    with open(section_file, 'r') as f:
        lines = f.readlines()
    A = float(lines[4].split()[1])
    I_z = float(lines[4].split()[4])
    
    # Grid (to get beam length)
    grid_file = job_dir / "grid.txt"
    grid_data = pd.read_csv(grid_file, comment='#', skipinitialspace=True, sep=r'\s+', 
                            names=['node_id', 'x', 'y', 'z'], skiprows=3)
    x_coords = grid_data['x'].values
    L = x_coords[-1] - x_coords[0]
    n_nodes = len(x_coords)
    
    # Point load
    point_load_file = job_dir / "point_load.txt"
    with open(point_load_file, 'r') as f:
        lines = f.readlines()
    load_line = lines[4].split()
    P = float(load_line[4])
    
    print(f"\nInput Parameters:")
    print(f"  E = {E:.2e} Pa")
    print(f"  G = {G:.2e} Pa")
    print(f"  I_z = {I_z:.2e} m^4")
    print(f"  A = {A:.6f} m^2")
    print(f"  EI_z = {E * I_z:.2e} N*m^2")
    print(f"  L = {L:.3f} m")
    print(f"  P = {P:.1f} N (point load at tip)")
    
    # Read computed results
    U_file = Path(results_dir) / "primary_results" / "global" / "U_global.csv"
    if not U_file.exists():
        print(f"\n[ERROR] Results file not found: {U_file}")
        return False, None, None
    
    U_df = pd.read_csv(U_file)
    
    # Find tip node
    tip_node = n_nodes - 1
    tip_u_y_dof = 6 * tip_node + 1
    tip_theta_z_dof = 6 * tip_node + 5
    
    # Get computed values
    computed_u_y = U_df[U_df['Global DOF'] == tip_u_y_dof]['Value'].values[0]
    computed_theta_z = U_df[U_df['Global DOF'] == tip_theta_z_dof]['Value'].values[0]
    
    # Analytical solutions based on element type
    if element_type == "EulerBernoulliBeamElement3D":
        analytical_u_y = analytical_cantilever_tip_deflection_eb(P, E, I_z, L)
        analytical_theta_z = analytical_cantilever_rotation(P, E, I_z, L)
        shear_contribution = 0.0
    elif element_type == "TimoshenkoBeamElement3D":
        kappa = 5.0 / 6.0  # Shear correction factor
        analytical_u_y, delta_bending, delta_shear = analytical_timoshenko_tip_deflection(
            P, E, I_z, L, G, A, kappa
        )
        analytical_theta_z = analytical_cantilever_rotation(P, E, I_z, L)
        shear_contribution = delta_shear
    elif element_type == "LevinsonBeamElement3D":
        analytical_u_y, delta_bending, delta_shear = analytical_levinson_tip_deflection(
            P, E, I_z, L, G, A
        )
        analytical_theta_z = analytical_cantilever_rotation(P, E, I_z, L)
        shear_contribution = delta_shear
    elif element_type in ("BarElement3D", "TrussElement3D"):
        # Sanity check only: no analytical deflection comparison
        print(f"\nTip Node (Node {tip_node}, x = {x_coords[tip_node]:.3f} m):")
        print(f"  DOF {tip_u_y_dof} (u_y): {computed_u_y*1000:.6f} mm (sanity check, no analytical)")
        print(f"  DOF {tip_theta_z_dof} (theta_z): {computed_theta_z*180/np.pi:.6f} deg")
        fixed_u_y_dof = 1
        fixed_theta_z_dof = 5
        fixed_u_y = U_df[U_df['Global DOF'] == fixed_u_y_dof]['Value'].values[0]
        fixed_theta_z = U_df[U_df['Global DOF'] == fixed_theta_z_dof]['Value'].values[0]
        print(f"\nBoundary Conditions (Node 0): u_y={fixed_u_y:.2e}, theta_z={fixed_theta_z:.2e}")
        bc_ok = abs(fixed_u_y) < 1e-12 and abs(fixed_theta_z) < 1e-12
        if bc_ok:
            print("\n[PASS] OK (sanity): primary results loaded, BCs zero at support.")
        else:
            print("\n[WARNING] OK (sanity): primary results loaded; support BCs non-zero.")
        return True, computed_u_y, None
    else:
        print(f"\n[ERROR] Unknown element type: {element_type}")
        return False, None, None

    print(f"\nTip Node (Node {tip_node}, x = {x_coords[tip_node]:.3f} m):")
    print(f"  DOF {tip_u_y_dof} (u_y):")
    print(f"    Computed:   {computed_u_y*1000:.6f} mm")
    print(f"    Analytical: {analytical_u_y*1000:.6f} mm")
    print(f"    Error:      {abs(computed_u_y - analytical_u_y)*1000:.6f} mm")
    print(f"    Relative:   {abs((computed_u_y - analytical_u_y) / analytical_u_y) * 100:.4f}%")
    if element_type != "EulerBernoulliBeamElement3D":
        print(f"    Shear contrib: {shear_contribution*1000:.6f} mm ({shear_contribution/analytical_u_y*100:.2f}%)")
    
    print(f"\n  DOF {tip_theta_z_dof} (theta_z):")
    print(f"    Computed:   {computed_theta_z*180/np.pi:.6f} deg")
    print(f"    Analytical: {analytical_theta_z*180/np.pi:.6f} deg")
    print(f"    Error:      {abs(computed_theta_z - analytical_theta_z)*180/np.pi:.6f} deg")
    print(f"    Relative:   {abs((computed_theta_z - analytical_theta_z) / analytical_theta_z) * 100:.4f}%")
    
    # Check boundary conditions at fixed end (node 0)
    fixed_u_y_dof = 1
    fixed_theta_z_dof = 5
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
        return True, computed_u_y, analytical_u_y
    else:
        print("[FAIL] Some verifications failed:")
        if not u_y_ok:
            print(f"  Tip deflection error: {abs(computed_u_y - analytical_u_y)*1000:.6f} mm > {tol*1000:.3f} mm")
        if not theta_z_ok:
            print(f"  Tip rotation error: {abs(computed_theta_z - analytical_theta_z)*180/np.pi:.6f} deg")
        if not bc_ok:
            print(f"  Boundary conditions not properly fixed")
        return False, computed_u_y, analytical_u_y

if __name__ == "__main__":
    # Find latest results directories
    results_base = PROJECT_ROOT / "post_processing" / "results"
    
    jobs_config = [
        ("job_0000", "EulerBernoulliBeamElement3D"),
        ("job_0001", "TimoshenkoBeamElement3D"),
        ("job_0002", "LevinsonBeamElement3D")
    ]
    
    results = {}
    all_passed = True
    
    for job_name, element_type in jobs_config:
        job_dirs = list(results_base.glob(f"{job_name}_*"))
        if not job_dirs:
            print(f"[ERROR] No results found for {job_name}")
            all_passed = False
            continue
        
        latest_dir = max(job_dirs, key=lambda p: p.stat().st_mtime)
        print(f"\nUsing results directory: {latest_dir.name}")
        
        success, computed, analytical = verify_job_results(job_name, latest_dir, element_type)
        results[job_name] = {
            'element_type': element_type,
            'success': success,
            'computed': computed,
            'analytical': analytical
        }
        if not success:
            all_passed = False
    
    # Comparison summary
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY: All Three Element Formulations")
    print("=" * 70)
    
    if all_passed and len(results) == 3:
        eb = results['job_0000']['computed']
        tim = results['job_0001']['computed']
        lev = results['job_0002']['computed']
        
        print(f"\nTip Deflection Comparison (all negative, showing magnitudes):")
        print(f"  Euler-Bernoulli: {abs(eb)*1000:.6f} mm (bending only)")
        print(f"  Timoshenko:      {abs(tim)*1000:.6f} mm (bending + shear with kappa)")
        print(f"  Levinson:        {abs(lev)*1000:.6f} mm (bending + shear without kappa)")
        print()
        print(f"Expected relationship: |Timoshenko| > |Levinson| > |Euler-Bernoulli|")
        print(f"  (More negative = larger deflection)")
        print()
        
        # Check relationship
        if tim < lev < eb:  # All negative, so tim is most negative
            print("[PASS] Relationship is correct:")
            print(f"  |Timoshenko| ({abs(tim)*1000:.6f} mm) > |Levinson| ({abs(lev)*1000:.6f} mm) > |EB| ({abs(eb)*1000:.6f} mm)")
        else:
            print("[WARNING] Relationship may be unexpected")
            print(f"  Timoshenko: {tim*1000:.6f} mm")
            print(f"  Levinson:   {lev*1000:.6f} mm")
            print(f"  EB:         {eb*1000:.6f} mm")
        
        print()
        print(f"Shear contribution:")
        print(f"  Timoshenko vs EB: {(tim - eb)*1000:.6f} mm ({(tim/eb - 1)*100:.2f}% increase)")
        print(f"  Levinson vs EB:   {(lev - eb)*1000:.6f} mm ({(lev/eb - 1)*100:.2f}% increase)")
        print(f"  Levinson vs Tim:  {(lev - tim)*1000:.6f} mm ({(lev/tim - 1)*100:.2f}% difference)")
    
    print("\n" + "=" * 70)
    if all_passed:
        print("VERIFICATION SUMMARY: All results are CORRECT")
        print("=" * 70)
        sys.exit(0)
    else:
        print("VERIFICATION SUMMARY: Some results have ERRORS")
        print("=" * 70)
        sys.exit(1)

