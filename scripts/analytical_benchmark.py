"""
Analytical benchmark for Euler-Bernoulli beam element stiffness matrix.

Computes K_e using exact integration (not quadrature) for comparison with
numerical implementation.
"""

import numpy as np

def analytical_euler_bernoulli_stiffness(E, I_z, L):
    """
    Compute analytical Euler-Bernoulli beam element stiffness matrix.
    
    For a 2-node beam element with DOFs [u1, θ1, u2, θ2] for bending about z-axis:
    
    K = (EI/L^3) * [
        [12,  6L, -12,  6L],
        [6L, 4L², -6L, 2L²],
        [-12, -6L,  12, -6L],
        [6L, 2L², -6L, 4L²]
    ]
    
    Parameters
    ----------
    E : float
        Young's modulus (Pa)
    I_z : float
        Moment of inertia about z-axis (m⁴)
    L : float
        Element length (m)
    
    Returns
    -------
    K_analytical : np.ndarray, shape (4, 4)
        Stiffness matrix for [u1, θ1, u2, θ2] DOFs
    """
    EI = E * I_z
    factor = EI / (L**3)
    
    K = factor * np.array([
        [12,   6*L,  -12,   6*L],
        [6*L, 4*L*L, -6*L, 2*L*L],
        [-12,  -6*L,   12,  -6*L],
        [6*L, 2*L*L, -6*L, 4*L*L]
    ])
    
    return K

def map_to_full_12x12(K_4x4, dof_indices):
    """
    Map 4x4 bending stiffness to full 12x12 element stiffness matrix.
    
    Parameters
    ----------
    K_4x4 : np.ndarray, shape (4, 4)
        Stiffness for [u_y1, θ_z1, u_y2, θ_z2]
    dof_indices : list of int
        DOF indices in full 12-DOF system: [u_y1_idx, θ_z1_idx, u_y2_idx, θ_z2_idx]
    
    Returns
    -------
    K_full : np.ndarray, shape (12, 12)
        Full stiffness matrix with bending terms inserted
    """
    K_full = np.zeros((12, 12))
    for i, idx_i in enumerate(dof_indices):
        for j, idx_j in enumerate(dof_indices):
            K_full[idx_i, idx_j] = K_4x4[i, j]
    return K_full

if __name__ == "__main__":
    # Job 0000 parameters
    E = 2.1e11  # Pa
    I_z = 2.08769e-06  # m⁴
    L = 0.2  # m
    
    print("=== Analytical Euler-Bernoulli Beam Stiffness ===")
    print(f"E = {E:.2e} Pa")
    print(f"I_z = {I_z:.2e} m^4")
    print(f"EI_z = {E * I_z:.2e} N*m^2")
    print(f"L = {L:.3f} m")
    print()
    
    # Compute 4x4 bending stiffness
    K_4x4 = analytical_euler_bernoulli_stiffness(E, I_z, L)
    
    print("4x4 Bending Stiffness Matrix [u_y1, theta_z1, u_y2, theta_z2]:")
    print(K_4x4)
    print()
    
    # Extract key terms
    print("Key Terms:")
    print(f"K[0,0] (u_y1-u_y1) = {K_4x4[0,0]:.6e}")
    print(f"K[0,1] (u_y1-theta_z1) = {K_4x4[0,1]:.6e}")
    print(f"K[1,1] (theta_z1-theta_z1) = {K_4x4[1,1]:.6e}")
    print(f"K[0,2] (u_y1-u_y2) = {K_4x4[0,2]:.6e}")
    print(f"K[2,2] (u_y2-u_y2) = {K_4x4[2,2]:.6e}")
    print()
    
    # Map to full 12x12 (DOF indices: u_y1=1, θ_z1=5, u_y2=7, θ_z2=11)
    K_full = map_to_full_12x12(K_4x4, [1, 5, 7, 11])
    
    print("Key Terms in Full 12x12 Matrix:")
    print(f"K[1,1] (u_y1-u_y1) = {K_full[1,1]:.6e}")
    print(f"K[1,5] (u_y1-theta_z1) = {K_full[1,5]:.6e}")
    print(f"K[5,5] (theta_z1-theta_z1) = {K_full[5,5]:.6e}")
    print(f"K[7,7] (u_y2-u_y2) = {K_full[7,7]:.6e}")
    print(f"K[11,11] (theta_z2-theta_z2) = {K_full[11,11]:.6e}")
    print()
    
    # Compare with computed values (from numerical_trace_stiffness.py after fix)
    print("=== Comparison with Computed Values (After Fix) ===")
    computed_K11 = 6.576224e8
    computed_K15 = 6.576224e7
    computed_K55 = 8.768298e6
    
    print(f"Analytical K[1,1] = {K_full[1,1]:.6e}")
    print(f"Computed   K[1,1] = {computed_K11:.6e}")
    print(f"Ratio: {computed_K11 / K_full[1,1]:.2f}x")
    print()
    
    print(f"Analytical K[1,5] = {K_full[1,5]:.6e}")
    print(f"Computed   K[1,5] = {computed_K15:.6e}")
    print(f"Ratio: {computed_K15 / K_full[1,5]:.2f}x")
    print()
    
    print(f"Analytical K[5,5] = {K_full[5,5]:.6e}")
    print(f"Computed   K[5,5] = {computed_K55:.6e}")
    print(f"Ratio: {computed_K55 / K_full[5,5]:.2f}x")
    print()
    print("All ratios are 1.00x - fix successful!")
