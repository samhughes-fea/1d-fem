# Nonlinear buckling benchmark: imperfect Euler-type column

Pinned benchmark job: [`jobs/job_benchmark_nl_buckling_imperfect_column/`](../../jobs/job_benchmark_nl_buckling_imperfect_column/)

## Purpose

This benchmark provides a stable in-repository continuation case for nonlinear buckling with:

- nonlinear Euler-Bernoulli beam elements
- axial compression
- arc-length continuation enabled
- imperfection seeding from the first linear buckling mode

## Pinned repository scalar/path targets

For the current pinned benchmark configuration, the repository-calibrated acceptance targets are:

- first recorded load factor: `1.0`
- first recorded tip displacement: `0.0`
- continuation method string: `load_control` or `arc_length`

These are pinned repository targets for the current benchmark job and test harness, not a final literature-calibrated post-buckling closure.

## Acceptance criteria

1. The run produces [`continuation_history.csv`](../../simulation_runner/buckling/nonlinear_buckling_simulation.py).
2. The run produces [`nonlinear_buckling_summary.json`](../../simulation_runner/buckling/nonlinear_buckling_simulation.py).
3. The continuation history contains at least one converged increment.
4. The first row matches the pinned scalar/path targets within tolerance.
5. The continuation method token is recorded in the history.

## Tolerance policy

Use absolute tolerance `1e-12` on the first recorded load factor and tip displacement for the current benchmark harness.

## Current scope note

This benchmark is now repository-calibrated at the scalar/path level, but it is not yet a literature-tolerance closure benchmark. Promote it further only after a post-buckling reference scalar or path metric is frozen against an external source.
