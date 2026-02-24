# Roark formulae – verification visualisers

This directory contains Roark's beam formulas and scripts to compare them with FEM results.

## Layout (tidied)

| File / folder | Purpose |
|--------------|--------|
| `roarks_formulas_point.py` | Point-load cantilever: V, M, θ, u_y. Wrapper: `roark_point_load_response(x, L, E, I, P, load_type)`. |
| `roarks_formulas_distributed.py` | Distributed loads (UDL, triangular, parabolic): q, V, M, θ, u_y. Wrapper: `roark_distributed_load_response(x, L, E, I, w, load_type)`. |
| `roarks_formulas_visualiser.py` | Roark-only plots (no FEM): load intensities, point and distributed u/θ/V/M. Output: `plots/`. |
| `roark_verification.py` | **FEM vs Roark**: discovers results from `post_processing/results/job_*/...`, compares u_y and θ_z, prints errors, writes `verification/` (overlay PNGs + CSV). |
| `roark_section_forces_verification.py` | **FEM vs Roark V, M**: compares shear V_y and bending moment M_z from `tertiary_results` (nodal or gaussian) to Roark's V(x), M(x). Writes `verification/roark_section_forces_point_loads.png`, `roark_section_forces_distributed_loads.png`, `roark_section_forces_verification_data.csv`. |
| `verification/` | Outputs from `roark_verification.py` and `roark_section_forces_verification.py`: displacement/rotation and section-force overlay PNGs + CSVs. |
| `plots/` | Outputs from `roarks_formulas_visualiser.py`: load intensities and response plots. |
| `new/` | Legacy comparison outputs (optional; can be removed if unused). |

## Job mapping (README_JOBS.md)

- **Point loads (Euler–Bernoulli):** `job_0000` end, `job_0001` midspan, `job_0002` quarter.
- **Distributed loads:** `job_0005` UDL, `job_0006` triangular, `job_0007` parabolic.

`roark_verification.py` uses this mapping and the same discovery pattern as `deflection_tables/deformation_convergence.py` (glob `post_processing/results/job_*/primary_results/global/U_global.csv`).

## Running

From the **project root**:

```text
python post_processing/verification_visualisers/roarks_formulas/roark_verification.py
```

Requires FEM results under `post_processing/results/job_XXXX_.../primary_results/global/U_global.csv` for the jobs above. Prints max and RMS errors (u_y, θ_z) per job and writes plots and CSV under `roarks_formulas/verification/`.

Section forces (V, M) vs Roark:

```text
python post_processing/verification_visualisers/roarks_formulas/roark_section_forces_verification.py
```

Requires tertiary results (nodal or gaussian section forces) for the same jobs. Writes `verification/roark_section_forces_point_loads.png`, `roark_section_forces_distributed_loads.png`, `roark_section_forces_verification_data.csv`.

Roark-only plots (no FEM):

```text
python post_processing/verification_visualisers/roarks_formulas/roarks_formulas_visualiser.py
```

## Checking that FEM agrees with theory

1. Run the relevant jobs (point and/or distributed) so that `U_global.csv` exists for job_0000–0002 and/or job_0005–0007.
2. Run `roark_verification.py`. It evaluates Roark at FEM node positions for errors and on a fine grid for smooth overlay curves (same idea as `deformation_convergence.py`).
3. Inspect console output (max and RMS error per job) and the overlay figures in `verification/`. Use `roark_verification_data.csv` for detailed pointwise comparison.

Formula conventions: point load uses M=0 and θ constant for x > a; distributed uses closed-form V, M and integrated θ, u_y. Beam parameters (E, I_z, P, w) are set in `roark_verification.py` to match the deflection_tables/deformation_convergence and visualiser defaults.
