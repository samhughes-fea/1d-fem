"""
Compare global stiffness matrices between element types.

Extracts and compares key terms to identify scaling or transformation issues.
"""

import numpy as np
import pandas as pd
import sys
from pathlib import Path
from scipy.sparse import load_npz, csr_matrix

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

def find_latest_results(job_name):
    """Find latest results directory for a job."""
    results_base = PROJECT_ROOT / "post_processing" / "results"
    job_dirs = list(results_base.glob(f"{job_name}_*"))
    if not job_dirs:
        return None
    return max(job_dirs, key=lambda p: p.stat().st_mtime)

def read_global_stiffness(results_dir):
    """Read global stiffness matrix from results."""
    # Try CSV format (primary format)
    k_csv = results_dir / "primary_results" / "global" / "K_global.csv"
    if k_csv.exists():
        df = pd.read_csv(k_csv)
        # Format: 'Row (Global DOF)', 'Column (Global DOF)', 'K Value'
        if 'Row (Global DOF)' in df.columns and 'Column (Global DOF)' in df.columns and 'K Value' in df.columns:
            # Need to know size - check diagnostics or infer from max DOF
            total_dof = int(max(df['Row (Global DOF)'].max(), df['Column (Global DOF)'].max())) + 1
            
            K = csr_matrix((df['K Value'].values, 
                          (df['Row (Global DOF)'].values, df['Column (Global DOF)'].values)),
                          shape=(total_dof, total_dof))
            return K
    
    # Try npz format as fallback
    k_file = results_dir / "primary_results" / "global" / "K_global.npz"
    if k_file.exists():
        return load_npz(str(k_file))
    
    return None

def compare_global_stiffness():
    """Compare global stiffness matrices between element types."""
    print("=" * 70)
    print("GLOBAL STIFFNESS MATRIX COMPARISON")
    print("=" * 70)
    
    # Key DOFs to check
    tip_node = 10
    tip_global_dof_u_y = 6 * tip_node + 1  # u_y
    tip_global_dof_theta_z = 6 * tip_node + 5  # theta_z
    
    jobs = [
        ("job_0000", "EulerBernoulliBeamElement3D"),
        ("job_0001", "TimoshenkoBeamElement3D"),
        ("job_0002", "LevinsonBeamElement3D"),
    ]
    
    stiffness_matrices = {}
    
    for job_name, element_type in jobs:
        print(f"\n{job_name} ({element_type}):")
        results_dir = find_latest_results(job_name)
        if results_dir is None:
            print(f"  [ERROR] No results found")
            continue
        
        K = read_global_stiffness(results_dir)
        if K is None:
            print(f"  [ERROR] Could not read global stiffness matrix")
            continue
        
        stiffness_matrices[element_type] = K
        
        print(f"  K shape: {K.shape}")
        print(f"  K nnz: {K.nnz}")
        print(f"  K norm: {np.linalg.norm(K.data):.6e}")
        
        # Check key terms
        k_tip_tip = K[tip_global_dof_u_y, tip_global_dof_u_y]
        k_tip_rot = K[tip_global_dof_u_y, tip_global_dof_theta_z]
        k_rot_rot = K[tip_global_dof_theta_z, tip_global_dof_theta_z]
        
        print(f"  K[{tip_global_dof_u_y},{tip_global_dof_u_y}] (tip u_y - tip u_y): {k_tip_tip:.6e}")
        print(f"  K[{tip_global_dof_u_y},{tip_global_dof_theta_z}] (tip u_y - tip theta_z): {k_tip_rot:.6e}")
        print(f"  K[{tip_global_dof_theta_z},{tip_global_dof_theta_z}] (tip theta_z - tip theta_z): {k_rot_rot:.6e}")
    
    # Compare between element types
    if len(stiffness_matrices) == 3:
        print(f"\n" + "-" * 70)
        print("COMPARISON")
        print("-" * 70)
        
        K_eb = stiffness_matrices['EulerBernoulliBeamElement3D']
        K_tim = stiffness_matrices['TimoshenkoBeamElement3D']
        K_lev = stiffness_matrices['LevinsonBeamElement3D']
        
        # Compare tip-tip stiffness
        k_eb_tip = K_eb[tip_global_dof_u_y, tip_global_dof_u_y]
        k_tim_tip = K_tim[tip_global_dof_u_y, tip_global_dof_u_y]
        k_lev_tip = K_lev[tip_global_dof_u_y, tip_global_dof_u_y]
        
        print(f"\nTip u_y - u_y stiffness:")
        print(f"  Euler-Bernoulli: {k_eb_tip:.6e}")
        print(f"  Timoshenko:      {k_tim_tip:.6e}")
        print(f"  Levinson:       {k_lev_tip:.6e}")
        
        # Expected: Timoshenko and Levinson should be slightly less stiff (more flexible)
        # due to shear deformation
        ratio_tim = k_tim_tip / k_eb_tip if k_eb_tip != 0 else 0
        ratio_lev = k_lev_tip / k_eb_tip if k_eb_tip != 0 else 0
        
        print(f"  Ratio (Timoshenko/EB): {ratio_tim:.6f}")
        print(f"  Ratio (Levinson/EB):   {ratio_lev:.6f}")
        
        if ratio_tim > 0.9 and ratio_tim < 1.1:
            print(f"  [OK] Timoshenko stiffness is reasonable")
        else:
            print(f"  [WARNING] Timoshenko stiffness ratio is unexpected")
        
        if ratio_lev > 0.9 and ratio_lev < 1.1:
            print(f"  [OK] Levinson stiffness is reasonable")
        else:
            print(f"  [WARNING] Levinson stiffness ratio is unexpected")
        
        # Check if there's a 300x scaling issue
        if abs(ratio_tim - 300) < 10 or abs(ratio_tim - 1/300) < 0.01:
            print(f"  [ERROR] Timoshenko stiffness appears to have 300x scaling issue!")
        if abs(ratio_lev - 500) < 10 or abs(ratio_lev - 1/500) < 0.01:
            print(f"  [ERROR] Levinson stiffness appears to have 500x scaling issue!")
        
        # Compare norms
        norm_eb = np.linalg.norm(K_eb.data)
        norm_tim = np.linalg.norm(K_tim.data)
        norm_lev = np.linalg.norm(K_lev.data)
        
        print(f"\nMatrix norms:")
        print(f"  Euler-Bernoulli: {norm_eb:.6e}")
        print(f"  Timoshenko:      {norm_tim:.6e}")
        print(f"  Levinson:       {norm_lev:.6e}")
        print(f"  Ratio (Timoshenko/EB): {norm_tim/norm_eb:.6f}")
        print(f"  Ratio (Levinson/EB):   {norm_lev/norm_eb:.6f}")
    
    print("\n" + "=" * 70)

if __name__ == "__main__":
    compare_global_stiffness()

