"""
Analytical benchmark for Timoshenko beam element tip deflection.

Timoshenko theory includes both bending and shear contributions:
delta = PL^3/(3EI) + PL/(kappa*G*A)
"""

import numpy as np

def analytical_timoshenko_tip_deflection(P, E, I_z, L, G, A, kappa):
    """
    Compute analytical Timoshenko tip deflection for cantilever beam.
    
    Parameters
    ----------
    P : float
        Point load at tip (N)
    E : float
        Young's modulus (Pa)
    I_z : float
        Moment of inertia about z-axis (m^4)
    L : float
        Beam length (m)
    G : float
        Shear modulus (Pa)
    A : float
        Cross-sectional area (m^2)
    kappa : float
        Shear correction factor (typically 5/6 for rectangular)
    
    Returns
    -------
    delta : float
        Tip deflection (m)
    delta_bending : float
        Bending contribution (m)
    delta_shear : float
        Shear contribution (m)
    """
    # Bending contribution (same as Euler-Bernoulli)
    delta_bending = P * L**3 / (3 * E * I_z)
    
    # Shear contribution (Timoshenko specific)
    delta_shear = P * L / (kappa * G * A)
    
    # Total deflection
    delta = delta_bending + delta_shear
    
    return delta, delta_bending, delta_shear

if __name__ == "__main__":
    # Job 0000 parameters
    P = 500  # N
    E = 2.1e11  # Pa
    I_z = 2.08769e-06  # m^4
    L = 2.0  # m
    G = 8.1e10  # Pa
    A = 0.00131  # m^2
    kappa = 5.0 / 6.0  # Shear correction factor
    
    delta, delta_bend, delta_shear = analytical_timoshenko_tip_deflection(
        P, E, I_z, L, G, A, kappa
    )
    
    print("=" * 70)
    print("ANALYTICAL TIMOSHENKO TIP DEFLECTION")
    print("=" * 70)
    print(f"Load P = {P:.1f} N")
    print(f"Length L = {L:.3f} m")
    print(f"E = {E:.2e} Pa")
    print(f"I_z = {I_z:.2e} m^4")
    print(f"EI_z = {E * I_z:.2e} N*m^2")
    print(f"G = {G:.2e} Pa")
    print(f"A = {A:.6f} m^2")
    print(f"kappa = {kappa:.6f}")
    print()
    print("=== DEFLECTION COMPONENTS ===")
    print(f"Bending contribution: {delta_bend*1000:.6f} mm")
    print(f"Shear contribution:   {delta_shear*1000:.6f} mm")
    print(f"Total deflection:     {delta*1000:.6f} mm")
    print()
    print("=== COMPARISON ===")
    print(f"Euler-Bernoulli (bending only): {delta_bend*1000:.6f} mm")
    print(f"Timoshenko (bending + shear):   {delta*1000:.6f} mm")
    print(f"Difference (shear effect):     {(delta - delta_bend)*1000:.6f} mm")
    print(f"Shear contribution %:          {(delta_shear/delta)*100:.2f}%")

