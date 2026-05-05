# Nonlinear static validation

This directory hosts Abaqus-oriented validation assets for [`NonlinearStaticSimulationRunner`](../../../../simulation_runner/static/nonlinear_static_simulation.py).

Suite target specification:
- [`NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md`](../../../docs/conventions/NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md)

Recommended substructure:
- `reference_cases/`
- `plots/`
- `tables/`
- `output/`

The first benchmark phase should focus on pinned nonlinear cantilever cases before path-grade FEM-vs-Abaqus comparison metrics are added.

First comparison contract target:
- tip displacement vs load level for [`jobs/job_benchmark_nl_static_cantilever_tip/`](../../../jobs/job_benchmark_nl_static_cantilever_tip/)

Common FEM-side artifact contract:
- `primary_results/nonlinear_static_validation/{job_name}_tip_load_history.csv`

Abaqus-side reference contract:
- `tip_load_history.csv` with `load_step`, `load_factor`, `tip_displacement`

Common suite-level Abaqus artifacts:
- `tip_load_history.csv`
- `U_global.csv`
- `rotation_source.txt`

First comparison tool:
- [`compare_nonlinear_static_tip_history.py`](compare_nonlinear_static_tip_history.py)

Family-local benchmark runner:
- [`run_nonlinear_static_validation.py`](run_nonlinear_static_validation.py)
