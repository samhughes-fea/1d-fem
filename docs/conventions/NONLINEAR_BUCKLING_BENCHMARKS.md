# Nonlinear buckling benchmarks and validation status

This note records the current nonlinear buckling validation posture for the continuation workflow implemented by [`NonlinearBucklingSimulationRunner`](../../simulation_runner/buckling/nonlinear_buckling_simulation.py).

## Scope

The current implementation supports:

- load-control continuation
- arc-length predictor metadata
- imperfection seeding from linear buckling mode shapes

It does **not** yet claim full post-critical arc-length corrector closure or branch switching.

## In-repository benchmark ladder

| Level | Coverage | Status | Reference |
|---|---|---|---|
| Dispatch | [`tests/test_nonlinear_buckling_dispatch_mvp.py`](../../tests/test_nonlinear_buckling_dispatch_mvp.py) | Implemented | runner wiring and artifact generation |
| Continuation smoke | [`tests/test_nonlinear_buckling_continuation_smoke.py`](../../tests/test_nonlinear_buckling_continuation_smoke.py) | Implemented | history CSV, arc-length metadata, imperfection seeding |
| Imperfect column continuation | [`jobs/job_benchmark_nl_buckling_imperfect_column/`](../../jobs/job_benchmark_nl_buckling_imperfect_column/), [`NONLINEAR_BUCKLING_IMPERFECT_COLUMN_BENCHMARK.md`](NONLINEAR_BUCKLING_IMPERFECT_COLUMN_BENCHMARK.md) | Implemented (acceptance level) | pinned job + repository acceptance criteria; not yet literature-calibrated |
| Published post-buckling benchmark | not yet added | Pending | future pinned `jobs/job_*` reference case |

## Recommended benchmark progression

### 1. Imperfect Euler column

Use a straight prismatic column with:

- pinned or cantilever boundary conditions fixed in the job docstring
- axial compression
- imperfection seeded from the first linear buckling mode

Acceptance for the current phase:

- run completes
- continuation history is written
- imperfection metadata is present in [`nonlinear_buckling_summary.json`](../../simulation_runner/buckling/nonlinear_buckling_simulation.py)
- displacement path is non-trivial and reproducible

### 2. Arc-length benchmark

After a full corrector is implemented, promote a benchmark requiring turning-point traversal. Until then, treat arc-length support as predictor-only metadata support.

### 3. Literature-calibrated post-buckling case

Only add a tolerance-bound benchmark when these are pinned:

- geometry
- material
- end conditions
- load definition
- response scalar used for comparison

## Repository policy

- Keep nonlinear buckling benchmark claims narrower than linear buckling claims.
- Carry forward the thin-walled lateral–torsional buckling caveat from [`MODAL_BUCKLING_LTB_VALIDATION.md`](MODAL_BUCKLING_LTB_VALIDATION.md).
- Prefer pinned jobs and explicit response scalars over qualitative screenshots.

## Current recommended CI set

- [`tests/test_nonlinear_buckling_dispatch_mvp.py`](../../tests/test_nonlinear_buckling_dispatch_mvp.py)
- [`tests/test_nonlinear_buckling_continuation_smoke.py`](../../tests/test_nonlinear_buckling_continuation_smoke.py)

These are smoke/regression checks, not final scientific validation closure.
