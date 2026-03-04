"""
Compare computed element stiffness matrix with analytical expectations.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

def read_stiffness_from_log(log_file):
    """Read stiffness matrix from log file."""
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    # Find the stiffness matrix
    start_idx = None
    for i, line in enumerate(lines):
        if 'Final Element Stiffness Matrix' in line:
            start_idx = i + 1
            break
    
    if start_idx is None:
        return None
    
    # Read the matrix (12x12)
    matrix_lines = []
    for i in range(start_idx, min(start_idx + 12, len(lines))):
        line = lines[i].strip()
        if line.startswith('[') or line.startswith('['):
            # Remove brackets and split
            line = line.replace('[', '').replace(']', '').strip()
            if line:
                matrix_lines.append(line)
    
    # Parse the matrix
    Ke = np.zeros((12, 12))
    row = 0
    for line in matrix_lines:
        if line:
            values = [float(x) for x in line.split()]
            if len(values) == 12:
                Ke[row, :] = values
                row += 1
            elif len(values) > 0:
                # Might be split across lines
                col_start = 0
                for val in values:
                    if col_start < 12:
                        Ke[row, col_start] = val
                        col_start += 1
                if col_start >= 12:
                    row += 1
    
    return Ke

def analytical_timoshenko_stiffness_terms(E, I_z, L, G, A, kappa):
    """Compute analytical Timoshenko stiffness terms for comparison."""
    EI = E * I_z
    kappa_GA = kappa * G * A
    
    # For a 2-node Timoshenko beam element, the stiffness matrix has both
    # bending and shear contributions
    # This is a simplified check - full matrix is complex
    
    # Bending contribution (similar to Euler-Bernoulli but with rotation-based curvature)
    # Shear contribution adds additional stiffness
    
    # Approximate: K_bending ~ 12*EI/L^3 (similar to EB)
    K_bending = 12 * EI / (L**3)
    
    # Shear contribution: K_shear ~ kappa*GA/L (for shear deformation)
    K_shear = kappa_GA / L
    
    # Combined effect makes Timoshenko stiffer than EB in some terms
    # but more flexible overall due to shear deformation
    
    return K_bending, K_shear

if __name__ == "__main__":
    # Find latest Timoshenko results
    results_base = PROJECT_ROOT / "post_processing" / "results"
    job_dirs = list(results_base.glob("job_0001_*"))
    if not job_dirs:
        print("No Timoshenko results found")
        sys.exit(1)
    
    latest_dir = max(job_dirs, key=lambda p: p.stat().st_mtime)
    log_file = latest_dir / "element_stiffness_matrices" / "stiffness_element_0.log"
    
    print("=" * 70)
    print("COMPARING COMPUTED vs ANALYTICAL STIFFNESS")
    print("=" * 70)
    
    Ke = read_stiffness_from_log(log_file)
    if Ke is None:
        print("Could not read stiffness matrix from log")
        sys.exit(1)
    
    print(f"\nComputed Stiffness Matrix (from log):")
    print(f"  Ke[1, 1] (u_y1-u_y1): {Ke[1, 1]:.6e}")
    print(f"  Ke[1, 5] (u_y1-theta_z1): {Ke[1, 5]:.6e}")
    print(f"  Ke[5, 5] (theta_z1-theta_z1): {Ke[5, 5]:.6e}")
    print(f"  Ke[7, 7] (u_y2-u_y2): {Ke[7, 7]:.6e}")
    print(f"  Ke[11, 11] (theta_z2-theta_z2): {Ke[11, 11]:.6e}")
    
    # Element properties
    L = 0.2  # m
    E = 2.1e11  # Pa
    I_z = 2.08769e-06  # m^4
    G = 8.1e10  # Pa
    A = 0.00131  # m^2
    kappa = 5.0 / 6.0
    
    K_bending, K_shear = analytical_timoshenko_stiffness_terms(E, I_z, L, G, A, kappa)
    
    print(f"\nAnalytical Estimates:")
    print(f"  K_bending ~ 12*EI/L^3: {K_bending:.6e}")
    print(f"  K_shear ~ kappa*GA/L: {K_shear:.6e}")
    
    # Compare
    print(f"\nComparison:")
    print(f"  Ke[1,1] / K_bending: {Ke[1, 1] / K_bending:.2f}x")
    print(f"  Ke[1,1] / K_shear: {Ke[1, 1] / K_shear:.2f}x")
    
    # Check if stiffness is too high
    # For a cantilever with point load, δ = PL^3/(3EI) + PL/(κGA)
    # The effective stiffness at tip is K_eff = P/δ
    # If K is too high, δ will be too small
    
    print(f"\n" + "=" * 70)

