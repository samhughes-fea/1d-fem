"""
Compare FEM results for all three element types against their respective analytical solutions.
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

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

def compare_job(job_name, results_dir, element_type):
    """Compare FEM results to analytical solution for a single job."""
    print("=" * 70)
    print(f"COMPARISON: {job_name} ({element_type})")
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
    
    # Grid
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
    
    # Read computed results
    U_file = Path(results_dir) / "primary_results" / "global" / "U_global.csv"
    U_df = pd.read_csv(U_file)
    
    # Find tip node
    tip_node = n_nodes - 1
    tip_u_y_dof = 6 * tip_node + 1
    tip_theta_z_dof = 6 * tip_node + 5
    
    # Get computed values
    computed_u_y = U_df[U_df['Global DOF'] == tip_u_y_dof]['Value'].values[0]
    computed_theta_z = U_df[U_df['Global DOF'] == tip_theta_z_dof]['Value'].values[0]
    
    # Analytical solutions based on element type
    print(f"\nInput Parameters:")
    print(f"  P = {P:.1f} N")
    print(f"  L = {L:.3f} m")
    print(f"  E = {E:.2e} Pa")
    print(f"  G = {G:.2e} Pa")
    print(f"  I_z = {I_z:.2e} m^4")
    print(f"  A = {A:.6f} m^2")
    
    if element_type == "EulerBernoulliBeamElement3D":
        analytical_u_y = analytical_cantilever_tip_deflection_eb(P, E, I_z, L)
        analytical_theta_z = analytical_cantilever_rotation(P, E, I_z, L)
        print(f"\nAnalytical Solution (Euler-Bernoulli):")
        print(f"  delta = PL^3/(3EI) = {analytical_u_y*1000:.6f} mm")
        
    elif element_type == "TimoshenkoBeamElement3D":
        kappa = 5.0 / 6.0
        analytical_u_y, delta_bending, delta_shear = analytical_timoshenko_tip_deflection(
            P, E, I_z, L, G, A, kappa
        )
        analytical_theta_z = analytical_cantilever_rotation(P, E, I_z, L)
        print(f"\nAnalytical Solution (Timoshenko):")
        print(f"  Bending: PL^3/(3EI) = {delta_bending*1000:.6f} mm")
        print(f"  Shear:   PL/(kappa*G*A) = {delta_shear*1000:.6f} mm (kappa = {kappa:.4f})")
        print(f"  Total:   delta = {analytical_u_y*1000:.6f} mm")
        
    elif element_type == "LevinsonBeamElement3D":
        analytical_u_y, delta_bending, delta_shear = analytical_levinson_tip_deflection(
            P, E, I_z, L, G, A
        )
        analytical_theta_z = analytical_cantilever_rotation(P, E, I_z, L)
        print(f"\nAnalytical Solution (Levinson):")
        print(f"  Bending: PL^3/(3EI) = {delta_bending*1000:.6f} mm")
        print(f"  Shear:   PL/(G*A) = {delta_shear*1000:.6f} mm (no kappa)")
        print(f"  Total:   delta = {analytical_u_y*1000:.6f} mm")
    
    print(f"\nFEM Results:")
    print(f"  Computed u_y:   {computed_u_y*1000:.6f} mm")
    print(f"  Computed theta_z:   {computed_theta_z*180/np.pi:.6f} deg")
    
    print(f"\nComparison:")
    print(f"  u_y error:      {abs(computed_u_y - analytical_u_y)*1000:.6f} mm")
    print(f"  u_y relative:   {abs((computed_u_y - analytical_u_y) / analytical_u_y) * 100:.4f}%")
    print(f"  theta_z error:      {abs(computed_theta_z - analytical_theta_z)*180/np.pi:.6f} deg")
    print(f"  theta_z relative:   {abs((computed_theta_z - analytical_theta_z) / analytical_theta_z) * 100:.4f}%")
    
    # Check if match
    tol = 1e-6
    u_y_match = abs(computed_u_y - analytical_u_y) < tol
    theta_z_match = abs(computed_theta_z - analytical_theta_z) < tol
    
    print(f"\n{'='*70}")
    if u_y_match and theta_z_match:
        print("[PASS] FEM results match analytical solution!")
    else:
        print("[FAIL] FEM results do NOT match analytical solution")
        if not u_y_match:
            print(f"  u_y mismatch: {abs(computed_u_y - analytical_u_y)*1000:.6f} mm > {tol*1000:.3f} mm")
        if not theta_z_match:
            print(f"  theta_z mismatch: {abs(computed_theta_z - analytical_theta_z)*180/np.pi:.6f} deg")
    print("=" * 70)
    
    return u_y_match and theta_z_match, computed_u_y, analytical_u_y

if __name__ == "__main__":
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
        print(f"\nUsing results: {latest_dir.name}\n")
        
        success, computed, analytical = compare_job(job_name, latest_dir, element_type)
        results[job_name] = {
            'element_type': element_type,
            'success': success,
            'computed': computed,
            'analytical': analytical
        }
        if not success:
            all_passed = False
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: Comparison of All Three Element Types")
    print("=" * 70)
    
    if len(results) == 3:
        eb = results['job_0000']
        tim = results['job_0001']
        lev = results['job_0002']
        
        print(f"\nTip Deflection Results:")
        print(f"  Euler-Bernoulli:")
        print(f"    FEM:       {eb['computed']*1000:.6f} mm")
        print(f"    Analytical: {eb['analytical']*1000:.6f} mm")
        print(f"    Status:     {'[PASS]' if eb['success'] else '[FAIL]'}")
        
        print(f"\n  Timoshenko:")
        print(f"    FEM:       {tim['computed']*1000:.6f} mm")
        print(f"    Analytical: {tim['analytical']*1000:.6f} mm")
        print(f"    Status:     {'[PASS]' if tim['success'] else '[FAIL]'}")
        
        print(f"\n  Levinson:")
        print(f"    FEM:       {lev['computed']*1000:.6f} mm")
        print(f"    Analytical: {lev['analytical']*1000:.6f} mm")
        print(f"    Status:     {'[PASS]' if lev['success'] else '[FAIL]'}")
        
        print(f"\nExpected Relationship (magnitude):")
        print(f"  |Timoshenko| > |Levinson| > |Euler-Bernoulli|")
        print(f"  |{tim['analytical']*1000:.6f}| > |{lev['analytical']*1000:.6f}| > |{eb['analytical']*1000:.6f}|")
        
        if tim['computed'] < lev['computed'] < eb['computed']:
            print(f"  FEM relationship: CORRECT")
        else:
            print(f"  FEM relationship: INCORRECT")
            print(f"    |{abs(tim['computed'])*1000:.6f}| vs |{abs(lev['computed'])*1000:.6f}| vs |{abs(eb['computed'])*1000:.6f}|")
    
    print("\n" + "=" * 70)
    if all_passed:
        print("OVERALL: All element types match their analytical solutions")
    else:
        print("OVERALL: Some element types do NOT match their analytical solutions")
    print("=" * 70)

