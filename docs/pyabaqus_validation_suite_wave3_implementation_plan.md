# Wave 3 pyAbaqus validation suite implementation plan and status

## Objective

Record Wave 3 implementation reality and closure evidence for the pyAbaqus validation suite. Wave 3 extends the suite to the nonlinear buckling simulation family, following the pattern established in [`docs/pyabaqus_validation_suite_wave1_implementation_plan.md`](docs/pyabaqus_validation_suite_wave1_implementation_plan.md) and [`docs/pyabaqus_validation_suite_wave2_implementation_plan.md`](docs/pyabaqus_validation_suite_wave2_implementation_plan.md).

Wave 3 families covered here:

- nonlinear buckling

Out of scope for Wave 3:

- full Abaqus-side nonlinear buckling ODB extraction (deferred to operational execution)
- post-buckling trajectory comparison beyond the pre-buckling regime sanity check

## Executive status summary

**Wave 3: CLOSED** as of 2026-05-17. All Phase 3 exit criteria are met.

### Wave 3 closure record

All four Phase 3 exit criteria were verified on 2026-05-17:

1. **Real comparator runner exists** — `run_nonlinear_buckling_validation` in [`validation/nonlinear_buckling/run_nonlinear_buckling_validation.py`](validation/nonlinear_buckling/run_nonlinear_buckling_validation.py) uses `build_runner_result` / `write_runner_summary` from `common_runner_contract.py`.
2. **Metric-based pass/fail semantics** — emits `pass`, `fail`, `missing_artifacts`, `analytical_only_pass`. The `analytical_only_pass` path checks that at least one increment converged and that the maximum FEM load factor remains below 99 % of the Euler critical load factor, using the analytical module in [`validation/nonlinear_buckling/analytical_nl_buckling.py`](validation/nonlinear_buckling/analytical_nl_buckling.py).
3. **Fixture-based behavioural tests** — [`tests/test_validation_visualisers.py`](tests/test_validation_visualisers.py) proves `pass`, `fail`, `missing_artifacts`, `analytical_only_pass`, and the fail-for-mismatched-Abaqus state for the nonlinear buckling family using fixture-style artifact setup. All tests pass under normal pytest without a live Abaqus installation.
4. **No placeholder-only dependencies** — no completion claim depends on scaffold-only summaries or placeholder-only exports. The `analytical_only_pass` path provides a non-Abaqus evidence tier and is itself proven by behavioural tests.

## Wave 3 family

### Nonlinear buckling

#### Canonical benchmark

- [`jobs/job_benchmark_nl_buckling_imperfect_column/`](jobs/job_benchmark_nl_buckling_imperfect_column/)
- Imperfect cantilever column with geometric imperfection; axial compressive load applied incrementally via arc-length/continuation method.

#### Abaqus reference contract

- Reference file: `continuation_history.csv` (columns: `load_factor`, `tip_displacement`, `converged`)
- Abaqus job: `job_benchmark_nl_buckling_imperfect_column_n500`

#### Analytical reference

- [`validation/nonlinear_buckling/analytical_nl_buckling.py`](validation/nonlinear_buckling/analytical_nl_buckling.py)
- Fixed-free (cantilever) Euler critical load: P_cr = π²EI_y / (4L²)
- Benchmark parameters: L = 2 m, E = 2.1 × 10¹¹ Pa, I_y ≈ 3.234 × 10⁻⁷ m⁴
- Benchmark result: P_cr ≈ 41 893 N, λ_cr ≈ 41.89 (P_applied = 1 000 N)
- `analytical_only_pass` tolerance: FEM max load factor < 99 % of λ_cr

#### Runner

- [`validation/nonlinear_buckling/run_nonlinear_buckling_validation.py`](validation/nonlinear_buckling/run_nonlinear_buckling_validation.py)
- Function: `run_nonlinear_buckling_validation(results_root)`
- FEM output path: `{results_root}/primary_results/nonlinear_buckling_results/continuation_history.csv`
- Abaqus reference path: `{ABAQUS_RESULTS_DIR}/job_benchmark_nl_buckling_imperfect_column_n500/continuation_history.csv`

#### Status vocabulary

| Status | Condition |
|---|---|
| `pass` | FEM and Abaqus present, max path interpolation error ≤ 1 × 10⁻⁴ |
| `fail` | Both present but error exceeds tolerance, or FEM present but analytical check fails |
| `missing_artifacts` | FEM CSV absent |
| `analytical_only_pass` | FEM present, Abaqus absent, ≥1 converged increment, max load factor < 99 % of λ_cr |

## Phase 3 exit criteria

Phase 3 should be treated as complete only when all of the following are true at the same time:

1. The nonlinear buckling family has a comparator runner that uses `build_runner_result` / `write_runner_summary` from `common_runner_contract.py`
2. The runner emits explicit metric-based status semantics (`pass`, `fail`, `missing_artifacts`, `analytical_only_pass`)
3. [`tests/test_validation_visualisers.py`](tests/test_validation_visualisers.py) proves all relevant runner states using fixture-style artifact setup
4. No completion claim depends on placeholder-only exports, scaffold-only summaries, or undocumented manual interpretation

If any one of those conditions is unmet, Phase 3 Wave 3 remains open.

## Testing coverage

All Wave 3 behavioural tests are in [`tests/test_validation_visualisers.py`](tests/test_validation_visualisers.py):

| Test | Family | State proved |
|---|---|---|
| `test_nonlinear_buckling_validation_runner_reports_missing_artifacts` | nonlinear_buckling | missing_artifacts |
| `test_nonlinear_buckling_validation_runner_reports_analytical_only_pass` | nonlinear_buckling | analytical_only_pass |
| `test_nonlinear_buckling_validation_runner_reports_fail_when_no_converged_increments` | nonlinear_buckling | fail (analytical check) |
| `test_nonlinear_buckling_validation_runner_reports_pass_for_matching_abaqus_reference` | nonlinear_buckling | pass |
| `test_nonlinear_buckling_validation_runner_reports_fail_for_mismatched_abaqus_reference` | nonlinear_buckling | fail (path mismatch) |
| `test_wave3_plan_declares_closure` | — | closure declaration |

## Design notes

### Why path interpolation for the pass criterion

The nonlinear buckling continuation history is a parametric load-displacement curve. Unlike time-history comparison (transient) or load-aligned histories (NL static), the load factor is the independent variable and differs between FEM and Abaqus meshes. Linear interpolation of the Abaqus tip displacement at each FEM load factor allows a mesh-independent comparison of the equilibrium path.

### Why 99 % Euler margin for analytical_only_pass

The pre-buckling check requires strict confidence that the FEM results are genuinely in the pre-buckling regime. Using 99 % of λ_cr (rather than 100 %) guards against numerical proximity to the bifurcation point, where linearised buckling theory loses reliability.

### Why _MAX_PATH_ERROR_TOL = 1e-4

Post-buckling path agreement is harder to achieve than linear static displacement matching. The looser tolerance (1e-4 vs 1e-5 for NL static) accounts for mesh-dependent path differences in the post-bifurcation regime.

## Deferred to future waves

- Full post-buckling trajectory comparison (Abaqus ODB extraction for n=500 continuation history)
- Multi-mesh GCI for nonlinear buckling families
- Imperfection sensitivity study automation
