"""
Verify that point loads are correctly applied to element and global force vectors.

Checks:
1. Element force vectors contain correct load values
2. Global force vector contains correct load values
3. Comparison between Euler-Bernoulli (working) and Timoshenko/Levinson
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

def find_latest_results(job_name):
    """Find latest results directory for a job."""
    results_base = PROJECT_ROOT / "post_processing" / "results"
    job_dirs = list(results_base.glob(f"{job_name}_*"))
    if not job_dirs:
        return None
    return max(job_dirs, key=lambda p: p.stat().st_mtime)

def read_element_force_vector(results_dir, element_id):
    """Read element force vector from log file."""
    log_file = results_dir / "element_force_vectors" / f"force_element_{element_id}.log"
    if not log_file.exists():
        return None
    
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    # Find "Final Force Vector"
    start_idx = None
    for i, line in enumerate(lines):
        if 'Final Force Vector' in line:
            start_idx = i + 1
            break
    
    if start_idx is None:
        return None
    
    # Read the vector (should be on next line)
    if start_idx < len(lines):
        line = lines[start_idx].strip()
        # Remove brackets and parse
        line = line.replace('[', '').replace(']', '').strip()
        values = [float(x) for x in line.split()]
        if len(values) == 12:
            return np.array(values)
    
    return None

def read_global_force_vector(results_dir):
    """Read global force vector from CSV."""
    f_file = results_dir / "primary_results" / "global" / "F_global.csv"
    if not f_file.exists():
        return None
    
    df = pd.read_csv(f_file)
    if 'Value' in df.columns:
        return df['Value'].values
    return None

def verify_load_application():
    """Verify load application for all three element types."""
    print("=" * 70)
    print("LOAD APPLICATION VERIFICATION")
    print("=" * 70)
    
    # Point load: -500 N at x=2.0m (tip, node 10)
    # Element 9 connects nodes 9-10, so load should be on element 9, node 10
    # Node 10: global DOF 61 = 6*10 + 1 (u_y)
    # Element 9: local DOF 7 = 6*1 + 1 (node 2, u_y)
    
    expected_load = -500.0  # N
    tip_node = 10
    tip_global_dof = 6 * tip_node + 1  # u_y DOF
    tip_element = 9
    tip_local_dof = 7  # node 2 (second node), u_y
    
    print(f"\nExpected Load:")
    print(f"  Point load: {expected_load} N at tip (node {tip_node}, x=2.0m)")
    print(f"  Global DOF: {tip_global_dof} (node {tip_node}, u_y)")
    print(f"  Element: {tip_element} (contains tip node)")
    print(f"  Local DOF: {tip_local_dof} (element node 2, u_y)")
    
    jobs = [
        ("job_0000", "EulerBernoulliBeamElement3D"),
        ("job_0001", "TimoshenkoBeamElement3D"),
        ("job_0002", "LevinsonBeamElement3D"),
    ]
    
    all_ok = True
    
    for job_name, element_type in jobs:
        print(f"\n" + "-" * 70)
        print(f"{job_name} ({element_type})")
        print("-" * 70)
        
        results_dir = find_latest_results(job_name)
        if results_dir is None:
            print(f"  [ERROR] No results found for {job_name}")
            all_ok = False
            continue
        
        print(f"  Results: {results_dir.name}")
        
        # Check element force vector
        Fe = read_element_force_vector(results_dir, tip_element)
        if Fe is None:
            print(f"  [ERROR] Could not read element force vector for element {tip_element}")
            all_ok = False
        else:
            print(f"  Element {tip_element} force vector:")
            print(f"    Fe[{tip_local_dof}] (node 2, u_y): {Fe[tip_local_dof]:.6f} N")
            if abs(Fe[tip_local_dof] - expected_load) < 1e-6:
                print(f"    [PASS] Element force vector correct")
            else:
                print(f"    [FAIL] Expected {expected_load} N, got {Fe[tip_local_dof]:.6f} N")
                print(f"    Error: {abs(Fe[tip_local_dof] - expected_load):.6f} N")
                all_ok = False
        
        # Check global force vector
        F_global = read_global_force_vector(results_dir)
        if F_global is None:
            print(f"  [ERROR] Could not read global force vector")
            all_ok = False
        else:
            print(f"  Global force vector:")
            print(f"    F_global[{tip_global_dof}] (node {tip_node}, u_y): {F_global[tip_global_dof]:.6f} N")
            if abs(F_global[tip_global_dof] - expected_load) < 1e-6:
                print(f"    [PASS] Global force vector correct")
            else:
                print(f"    [FAIL] Expected {expected_load} N, got {F_global[tip_global_dof]:.6f} N")
                print(f"    Error: {abs(F_global[tip_global_dof] - expected_load):.6f} N")
                all_ok = False
    
    # Compare between element types
    print(f"\n" + "-" * 70)
    print("COMPARISON BETWEEN ELEMENT TYPES")
    print("-" * 70)
    
    results = {}
    for job_name, element_type in jobs:
        results_dir = find_latest_results(job_name)
        if results_dir:
            Fe = read_element_force_vector(results_dir, tip_element)
            F_global = read_global_force_vector(results_dir)
            if Fe is not None and F_global is not None:
                results[element_type] = {
                    'Fe': Fe[tip_local_dof],
                    'F_global': F_global[tip_global_dof]
                }
    
    if len(results) == 3:
        eb_fe = results['EulerBernoulliBeamElement3D']['Fe']
        tim_fe = results['TimoshenkoBeamElement3D']['Fe']
        lev_fe = results['LevinsonBeamElement3D']['Fe']
        
        print(f"  Element force vectors (Fe[{tip_local_dof}]):")
        print(f"    Euler-Bernoulli: {eb_fe:.6f} N")
        print(f"    Timoshenko:      {tim_fe:.6f} N")
        print(f"    Levinson:       {lev_fe:.6f} N")
        
        if abs(eb_fe - tim_fe) < 1e-6 and abs(eb_fe - lev_fe) < 1e-6:
            print(f"    [PASS] All element types have identical force vectors")
        else:
            print(f"    [FAIL] Force vectors differ between element types")
            all_ok = False
        
        eb_fg = results['EulerBernoulliBeamElement3D']['F_global']
        tim_fg = results['TimoshenkoBeamElement3D']['F_global']
        lev_fg = results['LevinsonBeamElement3D']['F_global']
        
        print(f"  Global force vectors (F_global[{tip_global_dof}]):")
        print(f"    Euler-Bernoulli: {eb_fg:.6f} N")
        print(f"    Timoshenko:      {tim_fg:.6f} N")
        print(f"    Levinson:       {lev_fg:.6f} N")
        
        if abs(eb_fg - tim_fg) < 1e-6 and abs(eb_fg - lev_fg) < 1e-6:
            print(f"    [PASS] All element types have identical global force vectors")
        else:
            print(f"    [FAIL] Global force vectors differ between element types")
            all_ok = False
    
    print("\n" + "=" * 70)
    if all_ok:
        print("[PASS] Load application verification passed")
    else:
        print("[FAIL] Load application verification failed")
    print("=" * 70)
    
    return all_ok

if __name__ == "__main__":
    verify_load_application()

