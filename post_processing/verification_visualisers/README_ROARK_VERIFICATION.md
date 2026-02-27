# Roark formulae – verification visualisers

This directory contains the orchestrator and two subtrees: **roark/** (Roark's beam formulas and FEM comparison scripts) and **grid_convergence_index/** (GCI/Richardson report). Only `run_all_verification_visualisers.py` lives at this level; all other scripts and outputs are under the subtrees.

## Layout

| Path | Purpose |
|------|--------|
| `run_all_verification_visualisers.py` | Orchestrator: runs all verification scripts under `roark/` and `grid_convergence_index/`. |
| **roark/** | |
| `roark/roark_utilities/` | Formula modules: `roarks_formulas_euler_bernoulli_point.py`, `roarks_formulas_euler_bernoulli_distributed.py`, `roarks_formulas_timoshenko_point.py`, `roarks_formulas_timoshenko_distributed.py`. |
| `roark/roark_verification.py` | **FEM vs Roark**: compares u_y and θ_z, writes `roark/deformation_plots/` (overlay PNGs + CSV). |
| `roark/roark_section_forces_verification.py` | **FEM vs Roark V, M**: writes `roark/section_forces_plots/` (overlay PNGs + CSV). |
| `roark/roarks_formulas_visualiser.py` | Roark-only plots (no FEM): writes `roark/deformation_plots/`. |
| `roark/deformation_convergence.py`, `roark/distributed_load_convergence.py` | Convergence plots vs Roark; output in `roark/deformation_plots/`. |
| `roark/deformation_plots/` | Outputs from roark_verification, deformation_convergence, distributed_load_convergence, roarks_formulas_visualiser. |
| `roark/section_forces_plots/` | Outputs from roark_section_forces_verification. |
| **grid_convergence_index/** | |
| `grid_convergence_index/gci_richardson_roark_report.py` | GCI/Richardson report (uses Roark as reference); output in `grid_convergence_index/gci_tables/`. |
| `grid_convergence_index/gci_tables/` | GCI report CSV, LaTeX, PDF. |

## Job mapping (README_JOBS.md)

- **Point loads (Euler–Bernoulli):** `job_0000` end, `job_0001` midspan, `job_0002` quarter.
- **Distributed loads:** `job_0003` UDL, `job_0004` triangular, `job_0005` parabolic.
- **Timoshenko jobs 6–11:** Same load cases as 0–5; verified using `timoshenko_point_load_response`, `timoshenko_distributed_load_response`. Jobs 6–8 = point end/mid/quarter, 9–11 = UDL/triangular/parabolic.

`roark/roark_verification.py` uses the same discovery pattern as `roark/deformation_convergence.py` (glob `post_processing/results/job_*/primary_results/global/U_global.csv`).

## Running

From the **project root** you can run all verification scripts:

```text
python post_processing/verification_visualisers/run_all_verification_visualisers.py
```

Or run individual scripts:

```text
python post_processing/verification_visualisers/roark/roark_verification.py
```

Requires FEM results under `post_processing/results/job_XXXX_.../primary_results/global/U_global.csv`. Writes plots and CSV under `roark/deformation_plots/`.

Section forces (V, M) vs Roark:

```text
python post_processing/verification_visualisers/roark/roark_section_forces_verification.py
```

Writes `roark/section_forces_plots/roark_section_forces_*_euler_bernoulli.png`, `roark_section_forces_*_timoshenko.png`, and `roark_section_forces_verification_data.csv`.

**Section force sign convention (SFD/BMD):** Positive shear → clockwise rotation of a small element. Positive bending moment → sagging (bottom fibre in tension); negative → hogging (top in tension). For downward load on a cantilever (fixed left, free at x=L), V_y &lt; 0 and M_z &lt; 0 (hogging at the support). See the docstring in `roark_section_forces_verification.py` for details.

Roark-only plots (no FEM):

```text
python post_processing/verification_visualisers/roark/roarks_formulas_visualiser.py
```

## Checking that FEM agrees with theory

1. Run the relevant jobs so that `U_global.csv` exists for job_0000–0002 and/or job_0003–0005 (Euler–Bernoulli), and optionally job_0006–0011 (Timoshenko).
2. Run `roark/roark_verification.py`. It evaluates Roark at FEM node positions and on a fine grid for overlay curves.
3. Inspect console output and the overlay figures in `roark/deformation_plots/`. Use `roark_verification_data.csv` there for detailed pointwise comparison.

Beam parameters (E, I_z, P, w) are set in the scripts to match job material/section and deformation_convergence defaults.

## Timoshenko vs Euler–Bernoulli

- **Euler–Bernoulli**: plane sections remain plane and normal; no shear deformation.
- **Timoshenko**: adds shear deflection; V, M, θ unchanged; deflection = u_EB + u_shear.

Use the Timoshenko modules when comparing to a Timoshenko (shear-deformable) FEM beam. Parameters: A, G, and optional k_s (default 5/6 for rectangular section).
