"""
Verify that the shape function fix resolved the tip deflection issue.
"""

import numpy as np

# Load results
data_old = np.genfromtxt(
    'post_processing/results/job_0000_2025-11-16_17-00-46-814529_pid7088_4840d006/primary_results/global/U_global.csv',
    delimiter=',', skip_header=1, usecols=1
)
data_new = np.genfromtxt(
    'post_processing/results/job_0000_2025-11-16_17-40-08-750345_pid16416_f1cd8117/primary_results/global/U_global.csv',
    delimiter=',', skip_header=1, usecols=1
)

# Theoretical tip deflection: δ = PL³/(3EI)
# P = 500 N, L = 2.0 m, E = 2.1e11 Pa, I_z = 2.08769e-06 m⁴
# EI = 4.38e5 N*m²
theoretical = (500 * (2.0)**3) / (3 * 4.38e5)  # meters
theoretical_mm = theoretical * 1000

print("=" * 70)
print("SHAPE FUNCTION FIX VERIFICATION")
print("=" * 70)
print()
print("=== BEFORE FIX (Old Simulation) ===")
print(f"  Tip deflection: {abs(data_old[61])*1000:.6f} mm")
print(f"  Theoretical:    {theoretical_mm:.6f} mm")
print(f"  Error:          {abs(data_old[61])/theoretical:.2f}x too small")
print()
print("=== AFTER FIX (New Simulation) ===")
print(f"  Tip deflection: {abs(data_new[61])*1000:.6f} mm")
print(f"  Theoretical:    {theoretical_mm:.6f} mm")
print(f"  Error:          {abs(data_new[61])/theoretical:.2f}x")
print()
print("=== IMPROVEMENT ===")
improvement = abs(data_new[61] / data_old[61])
print(f"  Improvement factor: {improvement:.2f}x")
print(f"  Expected:           ~9.4x")
print()
print("=== VERIFICATION ===")
error_percent = abs((abs(data_new[61]) - theoretical) / theoretical) * 100
if error_percent < 5:
    print(f"  SUCCESS: Tip deflection matches theoretical within 5%!")
    print(f"    Error: {error_percent:.2f}%")
elif error_percent < 10:
    print(f"  GOOD: Tip deflection matches theoretical within 10%")
    print(f"    Error: {error_percent:.2f}%")
else:
    print(f"  WARNING: Tip deflection error is {error_percent:.2f}%")
print()
print("=" * 70)

