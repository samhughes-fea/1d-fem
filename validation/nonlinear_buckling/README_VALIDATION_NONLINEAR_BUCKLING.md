# Nonlinear buckling validation

FEM vs Abaqus (and analytical Euler) comparison for the imperfect-column
nonlinear buckling benchmark.

## Entrypoints

- `run_nonlinear_buckling_validation.py` — canonical Wave 3 comparator runner
  (status vocabulary: `pass`, `fail`, `missing_artifacts`, `analytical_only_pass`)
- `analytical_nl_buckling.py` — cantilever Euler critical load factor from job inputs
- `gci_richardson_nl_buckling_validation.py` — Richardson/GCI across n-variants vs Abaqus n500

## Expected artifacts

- FEM: `primary_results/nonlinear_buckling_results/continuation_history.csv`
- Abaqus reference (when available):
  `validation/abaqus_results/job_benchmark_nl_buckling_imperfect_column_n500/continuation_history.csv`

## Standard family layout

- `reference_cases/`
- `plots/`
- `tables/`
- `output/`
