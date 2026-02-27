# Validation visualisers (FEM vs Abaqus)

This directory contains the orchestrator and three subtrees aligned with [verification_visualisers](../verification_visualisers): **deformation/** (u_y, θ_z profiles), **section_forces/** (Vy, Mz), and **grid_convergence_study/** (GCI/Richardson and convergence across meshes). Key comparison outputs use **n128 only** (highest grid density per base job). Only `run_all_validation_visualisers.py`, `run_all_abaqus_jobs.py`, `run_batch_validation.py`, and `job_discovery.py` live at this level; comparison scripts and their outputs are under the subtrees.

## Layout

| Path | Purpose |
|------|--------|
| `run_all_validation_visualisers.py` | Orchestrator: runs deformation, section_forces, and grid_convergence_study scripts. |
| `run_all_abaqus_jobs.py` | Generate Abaqus scripts and run Abaqus for listed or discovered jobs. |
| `run_batch_validation.py` | Batch: run Abaqus (optional), run comparisons, check outputs (default jobs: job_0000_n128, job_0005_n128). |
| `job_discovery.py` | Discover job names from post_processing/results. |
| **deformation/** | |
| `deformation/deformation_comparison.py` | FEM vs Abaqus u_y and θ_z; **n128 only**. Writes `deformation/deformation_plots/`. |
| `deformation/deformation_plots/` | Deformation overlay PNGs and deformation_comparison_errors.csv. |
| **section_forces/** | |
| `section_forces/section_forces_comparison.py` | FEM vs Abaqus Vy, Mz; **n128 only**. Writes `section_forces/section_forces_plots/`. |
| `section_forces/section_forces_plots/` | Section forces overlay PNGs. |
| **grid_convergence_study/** | |
| `grid_convergence_study/gci_richardson_abaqus_report.py` | GCI/Richardson vs Abaqus (n32, 64, 128); reference Abaqus n500. Writes `grid_convergence_study/gci_tables/`. |
| `grid_convergence_study/u_global_largest_mesh_review.py` | Largest-mesh per base job review. Writes to gci_tables/. |
| `grid_convergence_study/csv_to_latex_table.py` | LaTeX table from GCI report CSV. Writes to gci_tables/. |
| `grid_convergence_study/gci_tables/` | GCI CSV, LaTeX table, largest-mesh review CSV. |
| **abaqus/**, **abaqus_results/** | Abaqus script generation, run wrapper, ODB extraction; result dirs per job. |
| **output/** | Review only: abaqus_results_review.csv, abaqus_performance_summary.md, errors log (from review_abaqus_results.py). |

## Job mapping

- **Point loads (Euler–Bernoulli):** job_0000 (end), job_0001 (mid-span), job_0002 (quarter).
- **Distributed loads:** job_0005 (UDL), job_0006 (triangular), job_0007 (parabolic).
- **Mesh levels:** Key plots use FEM n128 vs **Abaqus n500** (converged reference). Grid convergence study uses 6 meshes (n4–n128); GCI report uses Abaqus n500 as reference. To create the n500 Abaqus batch: (1) run the three mesh variant scripts (they include n500), (2) run `python post_processing/validation_visualisers/run_all_abaqus_jobs.py --n500-only`.

## Running

From the **project root**:

```text
python post_processing/validation_visualisers/run_all_validation_visualisers.py
```

Or run individual scripts (see [README_VALIDATION_ABAQUS.md](README_VALIDATION_ABAQUS.md) for Abaqus setup and full workflow):

```text
python post_processing/validation_visualisers/deformation/deformation_comparison.py
python post_processing/validation_visualisers/section_forces/section_forces_comparison.py
python post_processing/validation_visualisers/grid_convergence_study/gci_richardson_abaqus_report.py
```

Requires FEM results under `post_processing/results/job_XXXX_nN_.../` and (for comparison curves) Abaqus **n500 reference** results under `validation_visualisers/abaqus_results/job_XXXX_n500/`. Deformation and section-forces compare FEM n128 to Abaqus n500. To build the n500 batch: run the three mesh variant scripts, then `run_all_abaqus_jobs.py --n500-only`. See [README_VALIDATION_ABAQUS.md](README_VALIDATION_ABAQUS.md) for the full workflow.
