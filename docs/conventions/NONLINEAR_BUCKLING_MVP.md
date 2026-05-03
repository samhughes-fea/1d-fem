# Nonlinear buckling runner (MVP)

## Status

[`NonlinearBucklingSimulationRunner`](../../simulation_runner/buckling/nonlinear_buckling_simulation.py) is an **orchestration shell** only. It does **not** implement arc-length, path-following, or limit-point solvers.

## Dispatch

When **`[Type] buckling`** (or legacy modal buckling) and **`[Buckling]`** contains **`nonlinear_buckling = true`**, [`process_job`](../../workflow_orchestrator/run_job.py) instantiates **`NonlinearBucklingSimulationRunner`** instead of [`LinearBucklingSimulationRunner`](../../simulation_runner/buckling/buckling_simulation.py).

Default is **`nonlinear_buckling = false`**: unchanged **linearized** buckling (generalized eigenproblem on **K** and **K_g** after prestress) via **`LinearBucklingSimulationRunner`** / [`VibrationBucklingBackend`](../../simulation_runner/spectral/vibration_buckling_backend.py).

## Relationship to nonlinear static prestress

**`buckling_prestress = nonlinear_static`** still refers to **how prestress displacements U are obtained** for **linear** buckling. It is **not** the same switch as **`nonlinear_buckling = true`**.

## MVP behaviour

`run()` creates **`diagnostics/nonlinear_buckling_mvp_stub.txt`** under the job results root and logs an informational message. Full solver work is tracked as follow-up product tasks.

## Follow-up (not MVP)

- Primary results layout and `primary_artifacts.json` parity with linear buckling.
- `processing.buckling` (or dedicated) nonlinear equilibrium / extended system kernels.
- Tests on reference instability problems.
