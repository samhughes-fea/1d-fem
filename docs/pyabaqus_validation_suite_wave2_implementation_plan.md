# Wave 2 pyAbaqus validation suite implementation plan and status

## Objective

Record Wave 2 implementation reality and closure evidence for the pyAbaqus validation suite. Wave 2 extends the suite to the nonlinear static and transient simulation families, deepening the already-started families identified in [`docs/pyabaqus_validation_suite_plan.md`](docs/pyabaqus_validation_suite_plan.md).

This document follows the same structure established in [`docs/pyabaqus_validation_suite_wave1_implementation_plan.md`](docs/pyabaqus_validation_suite_wave1_implementation_plan.md).

Wave 2 families covered here:

- nonlinear static
- transient

Out of scope for Wave 2:

- nonlinear buckling (Wave 3)
- full Abaqus-side transient ODB extraction (deferred to operational execution)
- multi-mesh GCI for transient families

## Executive status summary

**Wave 2: CLOSED** as of 2026-05-17. All Phase 2 exit criteria are met. Wave 3 planning may proceed.

### Wave 2 closure record

All four Phase 2 exit criteria were verified on 2026-05-17:

1. **Real comparator runners exist** — `run_nonlinear_static_validation` in [`validation/static/run_nonlinear_static_validation.py`](validation/static/run_nonlinear_static_validation.py) and `run_transient_validation` in [`validation/transient/run_transient_validation.py`](validation/transient/run_transient_validation.py) both use `build_runner_result` / `write_runner_summary` from `common_runner_contract.py`.
2. **Metric-based pass/fail semantics** — nonlinear static emits `pass`, `fail`, `missing_artifacts`; transient emits `pass`, `fail`, `missing_artifacts`, `analytical_only_pass`. The transient `analytical_only_pass` path compares FEM maximum tip displacement against the first-mode SDOF analytical steady-state amplitude from [`validation/transient/analytical_transient.py`](validation/transient/analytical_transient.py) within 50 % relative tolerance.
3. **Fixture-based behavioural tests** — [`tests/test_validation_visualisers.py`](tests/test_validation_visualisers.py) proves `pass`, `fail`, and `missing_artifacts` for nonlinear static, and `pass`, `fail`, `missing_artifacts`, and `analytical_only_pass` for transient using fixture-style artifact setup. All tests pass under normal pytest without a live Abaqus installation.
4. **No placeholder-only dependencies** — no completion claim depends on scaffold-only summaries or placeholder-only exports. The `analytical_only_pass` path provides a non-Abaqus evidence tier and is itself proven by behavioural tests.

## Wave 2 families

### Nonlinear static

#### Canonical benchmark

- [`jobs/job_benchmark_nl_static_cantilever_tip/`](jobs/job_benchmark_nl_static_cantilever_tip/)

#### Abaqus reference contract

- [`docs/conventions/NONLINEAR_STATIC_ABAQUS_REFERENCE_CONTRACT.md`](docs/conventions/NONLINEAR_STATIC_ABAQUS_REFERENCE_CONTRACT.md)
- Reference file: `tip_load_history.csv` (columns: `load_factor`, `tip_displacement`)
- Abaqus job: `job_benchmark_nl_static_cantilever_tip_n500`

#### Runner

- [`validation/static/run_nonlinear_static_validation.py`](validation/static/run_nonlinear_static_validation.py)
- Function: `run_nonlinear_static_validation(results_root)`
- Comparison module: [`validation/static/compare_nonlinear_static_tip_history.py`](validation/static/compare_nonlinear_static_tip_history.py)
- Tolerances: `max_abs_error ≤ 1e-5`, `max_load_factor_alignment_error ≤ 1e-12`

#### Status vocabulary

| Status | Condition |
|---|---|
| `pass` | FEM and Abaqus present, errors within tolerance |
| `fail` | FEM and Abaqus present, errors exceed tolerance |
| `missing_artifacts` | FEM or Abaqus CSV absent |

### Transient

#### Canonical benchmark

- [`jobs/job_benchmark_transient_cantilever_multidof/`](jobs/job_benchmark_transient_cantilever_multidof/)
- Sinusoidal tip load: F_y = −1 N at 5 Hz, duration 0.05 s

#### Reference contract

- FEM output: `tip_displacement_history.csv` (columns: `time`, `tip_displacement`)
- Abaqus reference (when available): `transient_reference.csv` (same columns)
- Abaqus job: `job_benchmark_transient_cantilever_multidof_n500`

#### Analytical reference

- [`validation/transient/analytical_transient.py`](validation/transient/analytical_transient.py)
- SDOF model: k = 3EI/L³, m_eff = 33/140·ρAL, ω_n = sqrt(k/m_eff)
- Analytical amplitude: (F₀/k) / |1 − (Ω/ω_n)²|
- For the pinned benchmark: ~6.26 × 10⁻⁶ m
- Tolerance for `analytical_only_pass`: 50 % relative

#### Runner

- [`validation/transient/run_transient_validation.py`](validation/transient/run_transient_validation.py)
- Function: `run_transient_validation(results_root)`

#### Status vocabulary

| Status | Condition |
|---|---|
| `pass` | FEM and Abaqus present, max abs error ≤ 1e-8 |
| `fail` | Both present but exceeds tolerance, or FEM present but analytical check fails |
| `missing_artifacts` | FEM CSV absent |
| `analytical_only_pass` | FEM present, Abaqus absent, FEM max abs within 50 % of SDOF analytical |

## Phase 2 Wave 2 exit criterion

Phase 2 Wave 2 should be treated as complete only when all of the following are true at the same time:

1. nonlinear static and transient each have a comparator runner that uses `build_runner_result` / `write_runner_summary` from `common_runner_contract.py`
2. runners emit explicit metric-based status semantics (`pass`, `fail`, `missing_artifacts`, and `analytical_only_pass` for transient)
3. [`tests/test_validation_visualisers.py`](tests/test_validation_visualisers.py) proves all relevant runner states for both families using fixture-style artifact setup
4. no completion claim depends on placeholder-only exports, scaffold-only summaries, or undocumented manual interpretation

If any one of those conditions is unmet, Phase 2 Wave 2 remains open.

## Testing coverage

All Wave 2 behavioural tests are in [`tests/test_validation_visualisers.py`](tests/test_validation_visualisers.py):

| Test | Family | State proved |
|---|---|---|
| `test_nonlinear_static_validation_runner_reports_missing_artifacts` | nonlinear_static | missing_artifacts |
| `test_nonlinear_static_validation_runner_reports_pass_for_matching_tip_histories` | nonlinear_static | pass |
| `test_nonlinear_static_validation_runner_reports_fail_for_mismatched_tip_histories` | nonlinear_static | fail |
| `test_transient_validation_runner_reports_missing_artifacts` | transient | missing_artifacts |
| `test_transient_validation_runner_reports_analytical_only_pass` | transient | analytical_only_pass |
| `test_transient_validation_runner_reports_fail_when_fem_displacement_unreasonable` | transient | fail |
| `test_transient_validation_runner_reports_pass_for_matching_abaqus_reference` | transient | pass |
| `test_wave2_plan_declares_closure` | — | closure declaration |

## Design notes

### Why no analytical_only_pass for nonlinear static

The elastica large-deflection solution is the nonlinear static analytical reference but it is also used as the primary verification visualiser rather than a simple closed-form benchmark. The tolerance policy for elastica vs FEM is harder to define without knowing the current load magnitude in the benchmark job. Nonlinear static closure is therefore defined as FEM-vs-Abaqus only, with `missing_artifacts` when Abaqus is absent.

### Transient comparison scope

The transient comparison uses tip vertical displacement (`tip_displacement`) as the sole comparison quantity. Full DOF time-history comparison against Abaqus is deferred to operational execution once Abaqus transient ODB extraction is implemented.

## Deferred to Wave 3

- nonlinear buckling runner and contract
- full-field transient DOF history comparison against extracted Abaqus ODB data
- multi-mesh GCI for transient families
