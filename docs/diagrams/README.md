# FEM pipeline diagrams

Diagrams in this folder document the finite element model pipeline: job orchestration, parsing, element creation, static/modal simulation flow, and results/post-processing.

## Diagram index

| Diagram | File | Description |
|--------|------|-------------|
| Pipeline overview | [pipeline_overview.md](pipeline_overview.md) | End-to-end flow from job discovery to results and post-processing.
| Static simulation flow | [static_simulation_flow.md](static_simulation_flow.md) | Stage-by-stage linear-static workflow (prepare → assemble → BCs → condense → solve → reconstruct → primary/secondary/tertiary results). Filename is historical (linear-static only; there is no `static_simulation` Python module).
| Data flow | [data_flow.md](data_flow.md) | Job input files and result output directories/artifacts.
| Component structure | [component_structure.md](component_structure.md) | Top-level modules (pre_processing, processing, simulation_runner, post_processing) and dependencies.

## Main entry points

- **Job orchestration:** `workflow_orchestrator/run_job.py` — `main()` discovers jobs under `jobs/`, creates result dirs, runs `process_job()` per job (optionally in parallel).
- **Per-job pipeline:** `workflow_orchestrator/run_job.py` — `process_job()` does parsing → element instantiation → K_e/F_e computation → selects the runner for `simulation_settings.type` → `runner.run()`.
- **Static simulation:** [`simulation_runner/static/linear_static_simulation.py`](../../simulation_runner/static/linear_static_simulation.py) — `LinearStaticSimulationRunner.run()` executes the full linear-static workflow (see [static_simulation_flow.md](static_simulation_flow.md)).
- **Eigen / buckling:** `simulation_runner/eigen/eigen_simulation.py` and `simulation_runner/buckling/buckling_simulation.py` — shared `VibrationBucklingBackend` in `simulation_runner/spectral/vibration_buckling_backend.py`.

Diagrams use [Mermaid](https://mermaid.js.org/) and render in GitHub and most Markdown viewers.
