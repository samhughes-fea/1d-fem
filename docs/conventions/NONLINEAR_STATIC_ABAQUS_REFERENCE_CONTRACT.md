# Nonlinear static Abaqus reference contract

Benchmark job:

- [`jobs/job_benchmark_nl_static_cantilever_tip/`](../../jobs/job_benchmark_nl_static_cantilever_tip/)

## Purpose

Define the Abaqus-side reference-file contract needed to compare nonlinear-static FEM results against Abaqus for the tip-displacement-vs-load benchmark.

## Required Abaqus reference quantity

For the current benchmark, Abaqus must provide a load-step-wise tip displacement history for the loaded tip node vertical displacement DOF.

## Reference-file contract

Expected Abaqus-side reference file:

- `tip_load_history.csv`

Expected columns:

- `load_step`
- `load_factor`
- `tip_displacement`

## Phase 3 scope

This phase defines the file contract only. It does not yet implement the final FEM-vs-Abaqus comparison script.

## Next phase

The next phase should compare:

- FEM-side [`{job_name}_tip_load_history.csv`](../../processing/static/results/save_nonlinear_static_validation_summary.py)
- Abaqus-side `tip_load_history.csv`

at pinned load levels with documented tolerance.
