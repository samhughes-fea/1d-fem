# Track 1 — Operational V&V Status

## Summary

All seven analysis families have runner scripts, GCI scripts, and ODB extractors (Waves 1–3 code-complete). This document tracks the operational steps: running FEM jobs, comparing against Abaqus references, and committing results.

---

## What Has Been Done

### Infrastructure fixes (committed: `628452e0`)
- `src/onedfem/simulation_runner/transient/dynamic_simulation.py` — three bugs fixed: `len()` for list elements, float64 cast in `_ensure_sparse_format`, `_write_tip_displacement_history` method added
- `workflow_orchestrator/run_job.py` — transient `dynamic_settings` dict expanded to full 15-key pattern
- `src/onedfem/simulation_runner/api.py` — `_prepare_results_dir` and `_ensure_integration_orders` helpers added; error message wording fixed
- `tests/test_geometric_stiffness_weak_form.py` — **(see K_sigma note below)**

### Transient FEM n-variants — DONE
All 8 variants (base + n4/8/16/32/64/128/500) ran successfully:
- Output: `post_src/onedfem/processing/results/job_benchmark_transient_cantilever_multidof*`
- Key artifact: `primary_results/transient/tip_displacement_history.csv` ✓

### NL buckling FEM n-variants — DONE
All 8 variants ran successfully:
- Output: `post_src/onedfem/processing/results/job_benchmark_nl_buckling_imperfect_column*`
- Key artifact: `primary_results/nonlinear_buckling_results/continuation_history.csv` ✓

### Linear static validation (GCI) — PASS
All base jobs pass vs Abaqus n=500 with <1% error (tip deflection and rotation). GCI indices are well within acceptable bounds.

### NL static runner fix
`validation/static/run_nonlinear_static_validation.py` — `run_suite` was passing `FEM_RESULTS_DIR` directly instead of per-job timestamped subdirectories. Fixed with `_find_latest_result_dir` helper.

`validation/static/compare_nonlinear_static_tip_history.py` — `infer_abaqus_reference_job_name` only remapped `_n64`; fixed to remap any `_nN` suffix to `_n500`.

### K_sigma geometric stiffness — REVERTED
A previous fix removed `M_z/L` and `M_y/L` from `assemble_K_sigma` in `geometric_stiffness.py`. V&V results showed this was wrong: post-fix FEM gave ~0.12 m vs Abaqus 0.376–0.423 m (~70% error) for all NL static families. Pre-fix n=64 results (which survived the file-descriptor crash) gave ~0.40 m vs Abaqus 0.42 m (~5% error).

**Root cause:** For the cantilever under transverse load, N≈0 but M_z is large and negative. Without the `M_z/L` term, `K_sigma ≈ 0`, so `K_T ≈ K_material` — the solver behaves like linear FEM and massively under-predicts deflection. With `s_xy = N + M_z/L`, the negative moment softens K_T, allowing physically correct large deflections.

**Status:** `geometric_stiffness.py` and `test_geometric_stiffness_weak_form.py` have been restored to the pre-fix formulation (`s_xy = Ng + Mzg/L`, `s_xz = Ng + Myg/L`). Tests pass.

---

## Current State of NL Static Re-runs

Jobs must be re-run one family at a time (running all 48 in a single batch hits the OS file-descriptor limit).

| Family | Re-run (post K_sigma revert) | FEM results present |
|---|---|---|
| cantilever_tip | ✅ Done | ✅ |
| midspan_point | ✅ Done | ✅ |
| parabolic | ✅ Done | ✅ |
| quarter_point | ⬜ Not started | — |
| triangular | ⬜ Not started | — |
| udl | ⬜ Not started | — |

---

## Remaining Steps

### Step A — Finish NL static FEM re-runs (3 families remaining)

Status note (2026-06-03): the nonlinear-static adaptive stepping evaluation has now closed for the exact-match quarter-point benchmark pair. The fixed-schedule control run converged through `λ=0.15` and failed at `λ=0.20`, while the adaptive run recovered continuation through cutback/retry logic and reached accepted convergence at `λ=0.175` before failing again at `λ=0.1875`. This closes phase 7 and supports moving into phase 8 using the adaptive controller structure already implemented in [`src/onedfem/simulation_runner/static/nonlinear_static_simulation.py`](src/onedfem/simulation_runner/static/nonlinear_static_simulation.py:1257).

### Pre-Phase-8 Hardening Gate — Nonlinear Static Adaptive Stepping

Phase 7 is functionally closed, but phase 8 should not begin until the nonlinear-static adaptive stepping work is hardened into a repeatable acceptance gate. The required hardening work is regression protection, diagnostics contract definition, solver-observability expectations, and rerun/evidence hygiene.

Reference evidence pair:
- [`post_src/onedfem/processing/results/job_benchmark_nl_static_quarter_point_n128_exact_control_2026-06-01_09-11-11-427239_pid25528_1f36b79d`](post_src/onedfem/processing/results/job_benchmark_nl_static_quarter_point_n128_exact_control_2026-06-01_09-11-11-427239_pid25528_1f36b79d)
- [`post_src/onedfem/processing/results/job_benchmark_nl_static_quarter_point_n128_exact_adaptive_2026-06-01_09-14-42-310743_pid9112_72ce1989`](post_src/onedfem/processing/results/job_benchmark_nl_static_quarter_point_n128_exact_adaptive_2026-06-01_09-14-42-310743_pid9112_72ce1989)

Regression gate:
- [ ] Add a regression test that asserts the control path converges only through `λ=0.15`
- [ ] Assert the control path first fails at `λ=0.20`
- [ ] Assert the adaptive path converges beyond control and reaches `λ=0.175`
- [ ] Assert the adaptive converged history is longer than the control history
- [ ] Assert adaptive recovery includes rejected attempt(s) followed by accepted smaller-step continuation

Diagnostics gate:
- [ ] Adaptive runs produce [`nonlinear_static_tip_trace.csv`](src/onedfem/simulation_runner/static/nonlinear_static_simulation.py:332)
- [ ] Adaptive runs produce [`adaptive_load_stepping_history.csv`](src/onedfem/simulation_runner/static/nonlinear_static_simulation.py:355)
- [ ] Adaptive runs produce [`newton_history.csv`](src/onedfem/simulation_runner/static/nonlinear_static_simulation.py:312)
- [ ] Adaptive runs produce [`SolveCondensedSystem.log`](src/onedfem/processing/static/operations/solver.py:267)
- [ ] Accepted/rejected attempt progression is reconstructable across the CSV artifacts
- [ ] Tip-trace semantics remain converged-path oriented with terminal failure context only

Solver-observability gate:
- [ ] Frontier non-convergence remains visible in [`SolveCondensedSystem.log`](src/onedfem/processing/static/operations/solver.py:267)
- [ ] Retry behavior remains visible in [`SolveCondensedSystem.log`](src/onedfem/processing/static/operations/solver.py:267)
- [ ] Direct fallback remains visible in [`SolveCondensedSystem.log`](src/onedfem/processing/static/operations/solver.py:267) when used
- [ ] End-of-solve residual reporting remains visible without requiring brittle exact numeric matching

Rerun hygiene gate:
- [ ] Use exact-match benchmark names only for gate evidence
- [ ] Run control/adaptive comparison jobs in isolation
- [ ] Select result roots by artifact completeness, not just recency
- [ ] Preserve the designated evidence pair even if newer exploratory runs are created
- [ ] Record future evidence-pair changes explicitly in [`TRACK1_STATUS.md`](TRACK1_STATUS.md)

Phase-8 entry condition:
- [ ] Phase 8 may begin only after the regression gate, diagnostics gate, solver-observability gate, and rerun-hygiene gate are satisfied and the designated evidence pair remains preserved as the nonlinear-static adaptive baseline.

Run one family at a time with `--serial`:

```bash
python -m workflow_orchestrator.run_job --serial --filter nl_static_quarter_point
python -m workflow_orchestrator.run_job --serial --filter nl_static_triangular
python -m workflow_orchestrator.run_job --serial --filter nl_static_udl
```

### Step B — Run NL static validation suite

```bash
python validation/static/run_nonlinear_static_validation.py
```

Expected: 42 ready (6 families × 7 mesh variants), all pass (or close — tolerance is `max_abs_error < 1e-5 m`).

### Step C — Run NL static GCI

```bash
python validation/static/gci_richardson_nl_static_validation.py
```

Expected: monotonic convergence and GCI indices < 5% for all 6 families.

### Step D — Fix transient/NL buckling validator import error

Both `run_transient_validation.py` and `run_nonlinear_buckling_validation.py` crash with `ModuleNotFoundError: No module named 'post_processing'` when invoked from `run_all_validation_visualisers.py`. They need the same `sys.path` injection as the NL static runner (lines 7–14 of `run_nonlinear_static_validation.py`). Check whether they already have it; if not, add it.

### Step E — Run all validators

```bash
python validation/run_all_validation_visualisers.py
```

Expected outcome per family:

| Family | Expected status |
|---|---|
| Linear static | `pass` |
| Harmonic | depends on Abaqus ref |
| Eigen | depends on Abaqus ref |
| Linear buckling | depends on Abaqus ref |
| NL static | `pass` (after Steps A–C) |
| Transient | `analytical_only_pass` (no Abaqus ref yet) |
| NL buckling | `analytical_only_pass` (no Abaqus ref yet) |

### Step F — Run GCI for transient and NL buckling

```bash
python validation/transient/gci_richardson_transient_validation.py
python validation/nonlinear_buckling/gci_richardson_nl_buckling_validation.py
```

These skip gracefully if the Abaqus n=500 reference is missing (expected — no licence).

### Step G — Abaqus n=500 references (requires licence)

```bash
python validation/run_all_abaqus_jobs.py \
    --jobs job_benchmark_transient_cantilever_multidof_n500 \
           job_benchmark_nl_buckling_imperfect_column_n500
```

### Step H — Commit

Stage only tracked result artefacts:

```bash
git add src/onedfem/pre_src/onedfem/processing/element_library/nonlinear/euler_bernoulli/total_lagrangian/utilities/geometric_stiffness.py
git add tests/test_geometric_stiffness_weak_form.py
git add validation/static/run_nonlinear_static_validation.py
git add validation/static/compare_nonlinear_static_tip_history.py
git add validation/abaqus_results/
git add validation/*/plots/gci_tables/
git add validation/static/gci_tables/
git commit -m "fix(nl-static): revert spurious K_sigma removal and repair V&V runner path discovery"
```

---

## Known Outstanding Issues

- **Abaqus references missing** for transient and NL buckling — requires Abaqus licence (Step G).
- **Harmonic/eigen/linear buckling GCI** — these skip because no n-variant FEM results exist for those families (would need to run `job_benchmark_harmonic_sdof_n128` etc. first).
- **`run_all_validation_visualisers.py` exit code 1** — caused by the transient/NL buckling import error (Step D).
