# Nonlinear static Abaqus validation plan

Pinned benchmark job:

- [`jobs/job_benchmark_nl_static_cantilever_tip/`](../../jobs/job_benchmark_nl_static_cantilever_tip/)

## Purpose

Establish the first Abaqus-oriented validation track for [`NonlinearStaticSimulationRunner`](../../simulation_runner/static/nonlinear_static_simulation.py).

## Phase 1 scope

This phase creates:

- a family-local validation home under [`post_processing/validation_visualisers/static/`](../../post_processing/validation_visualisers/static/)
- a pinned nonlinear-static cantilever benchmark job
- initial regression checks that the benchmark assets and directory structure exist

It does **not** yet implement full FEM-vs-Abaqus nonlinear path comparison.

See also the broader suite target in [`NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md`](NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md).

## Fine-reference preflight pair

The first convergence-grade nonlinear-static comparison pair is:

- FEM comparison mesh: [`job_benchmark_nl_static_cantilever_tip_n64`](../../jobs/job_benchmark_nl_static_cantilever_tip_n64/)
- Abaqus fine reference: [`job_benchmark_nl_static_cantilever_tip_n500`](../../jobs/job_benchmark_nl_static_cantilever_tip_n500/)
