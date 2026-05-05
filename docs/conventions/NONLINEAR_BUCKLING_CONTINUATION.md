# Nonlinear buckling continuation (MVP contract)

This document defines the initial contract for nonlinear buckling continuation routed through [`NonlinearBucklingSimulationRunner`](../../simulation_runner/buckling/nonlinear_buckling_simulation.py).

## Status

This is an MVP continuation contract with **load-control implemented** and **arc-length predictor support introduced behind settings**.

- Supported: incremental nonlinear equilibrium tracing with structured history output.
- Partially supported: arc-length predictor metadata and predictor step generation.
- Not yet supported: full arc-length corrector loop, branch switching, automatic bifurcation tracking.

## Settings

Use [`[Buckling]`](../../jobs/job_smoke_buckling/simulation_settings.txt) with:

- `nonlinear_buckling = true`
- `continuation_method = load_control` (default)
- `continuation_method = arc_length`
- `num_increments = <int>`
- `load_factors = a, b, c, ...` (optional explicit schedule)
- `arc_length_radius = <float>`
- `arc_length_alpha_scale = <float>`
- `imperfection_mode_index = <int>`
- `imperfection_scale = <float>`
- `imperfection_source = linear_buckling` (planned/currently supported source)
- `line_search = true|false` (optional)
- `line_search_max_backtracks = <int>`
- `line_search_shrink = <float>`

Newton tolerances are reused from [`[Newton]`](../../jobs/job_test_n25_nonlinear_eb/simulation_settings.txt).

## Primary outputs

Successful or partially converged runs write:

- `diagnostics/nonlinear_buckling_diagnostic.log`
- `diagnostics/nonlinear_buckling_summary.json`
- `primary_results/nonlinear_buckling_results/continuation_history.csv`

When `continuation_method = arc_length`, the history includes predictor load-factor metadata even though the corrector still reuses the current nonlinear equilibrium solve.

When imperfection seeding is enabled, the summary JSON records the mode index and scale used to perturb the initial displacement field from the linear buckling mode-shape export under [`primary_results/modal_results/`](../../simulation_runner/spectral/vibration_buckling_backend.py).

## CSV schema

`continuation_history.csv` columns:

- `increment_index`
- `load_factor`
- `converged`
- `iterations_used`
- `residual_norm`
- `tip_dof`
- `tip_displacement`

## Scope note

This workflow traces nonlinear equilibrium of supported beam meshes. It does **not** close the thin-walled lateral–torsional buckling modelling gap documented in [`MODAL_BUCKLING_LTB_VALIDATION.md`](MODAL_BUCKLING_LTB_VALIDATION.md).
