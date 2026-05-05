# Nonlinear static Abaqus validation suite

This document defines the target validation-suite structure for nonlinear static analysis so it can reach the same maturity level as the linear static cantilever validation set.

## Objective

Build a nonlinear-static validation suite that mirrors the linear static six-load-case cantilever family and compares 1D-FEM results against Abaqus reference results across mesh ladders.

## Canonical nonlinear-static load cases

The suite should contain nonlinear-static counterparts to the six linear cantilever verification families.

### Point-load families

1. tip point load
2. midspan point load
3. quarter-point point load

### Distributed-load families

4. UDL
5. triangular distributed load
6. parabolic distributed load

## Required suite capabilities

For each load family, the suite should ultimately support:

- pinned benchmark job root
- mesh ladder variants
- FEM-side nonlinear validation exports
- Abaqus-side reference exports
- family-level FEM-vs-Abaqus comparison outputs
- summary reporting across the family

## Required comparison quantities

At minimum, each family should validate:

- tip displacement vs load factor

## Common Abaqus-side reference artifact contract

For every nonlinear-static benchmark job and mesh variant, the Abaqus side should emit at minimum:

- `tip_load_history.csv`
- `U_global.csv`
- `rotation_source.txt`

### `tip_load_history.csv`

Required columns:
- `frame_index`
- `load_factor`
- `tip_displacement`

This is the suite-level Abaqus counterpart to the FEM-side `tip_load_history.csv`-style validation export.

## Common FEM-side validation artifact contract

For every nonlinear-static benchmark job and mesh variant, the FEM run should write:

- `primary_results/nonlinear_static_validation/{job_name}_tip_load_history.csv`

with columns:

- `load_factor`
- `tip_displacement`

Preferably, later phases should also add:

- selected interior-node displacement vs load factor
- final deformed shape comparison
- reaction-force comparison where practical

## Mesh-ladder target

Recommended mesh levels:

- coarse: `n4`, `n8`
- medium: `n16`, `n32`, `n64`
- fine FEM review: `n128`
- fine Abaqus reference: `n500` where feasible

First planned fine Abaqus reference rollout target:
- [`job_benchmark_nl_static_cantilever_tip_n500`](../../jobs/job_benchmark_nl_static_cantilever_tip_n500/)

## Validation ladder within the suite

Each nonlinear-static load family should move through:

1. pinned acceptance
2. calibrated repository reference
3. external-reference validated against Abaqus

## First calibrated nonlinear-static gate

The first fully calibrated nonlinear-static fine-reference benchmark is:

- [`job_benchmark_nl_static_cantilever_tip_n64`](../../jobs/job_benchmark_nl_static_cantilever_tip_n64/)
- compared against fine Abaqus reference [`job_benchmark_nl_static_cantilever_tip_n500`](../../jobs/job_benchmark_nl_static_cantilever_tip_n500/)

This pair is now the first nonlinear-static family-local benchmark treated as calibrated and tolerance-backed.

### Frozen gate metrics

For the first-pass repository gate, the benchmark passes when:

- `max_load_factor_alignment_error <= 1e-12`
- `max_abs_error <= 1e-5`

These thresholds are implemented in [`evaluate_benchmark_pass_fail()`](../../post_processing/validation_visualisers/static/run_nonlinear_static_validation.py:13).

### First calibrated observed metrics

For the corrected `n64` FEM run against the repaired `n500` Abaqus reference, the observed metrics are:

- `max_load_factor_alignment_error = 0.0`
- `max_abs_error = 2.8889808153266788e-06`
- `mean_abs_error = 2.34892581072062e-06`

This places the first tip-load family inside the frozen tolerance gate and establishes it as the initial calibrated nonlinear-static repository reference.

### Phase-D operational meaning

The suite-level reporting path in [`run_suite()`](../../post_processing/validation_visualisers/static/run_nonlinear_static_validation.py:84) should now surface this calibrated case as:

- `ready = True`
- `passed = True`

Other nonlinear-static families remain outside the calibrated gate until they complete the same fine-reference rollout sequence.

## Additional calibrated distributed-load gates

The same first-pass tolerance gate has now also been satisfied by the first calibrated distributed-load nonlinear-static families.

### Triangular distributed load

Calibrated pair:

- [`job_benchmark_nl_static_triangular_n64`](../../jobs/job_benchmark_nl_static_triangular_n64/)
- compared against fine Abaqus reference [`job_benchmark_nl_static_triangular_n500`](../../jobs/job_benchmark_nl_static_triangular_n500/)

Observed metrics:

- `max_load_factor_alignment_error = 0.0`
- `max_abs_error = 4.8261920479980396e-06`
- `mean_abs_error = 3.6196440464633032e-06`

This places the triangular family inside the same frozen tolerance gate and makes it the second calibrated nonlinear-static family-local repository reference after the tip-load case.

### UDL distributed load

Calibrated pair:

- [`job_benchmark_nl_static_udl_n64`](../../jobs/job_benchmark_nl_static_udl_n64/)
- compared against fine Abaqus reference [`job_benchmark_nl_static_udl_n500`](../../jobs/job_benchmark_nl_static_udl_n500/)

Observed metrics:

- `max_load_factor_alignment_error = 0.0`
- `max_abs_error = 7.296619412120505e-06`
- `mean_abs_error = 5.47246455385932e-06`

This also lies inside the current frozen gate and establishes the UDL family as an additional calibrated nonlinear-static distributed-load benchmark.

## Phase-0 scope

This document is the suite specification only. It does not yet create the six full nonlinear-static job families or mesh ladders.
