# FEM pipeline diagrams

Diagrams in this folder document the finite element model pipeline: job orchestration, parsing, element creation, static/modal simulation flow, and results/post-processing.

## Diagram index

| Diagram | File | Description |
|--------|------|-------------|
| Pipeline overview | [pipeline_overview.md](pipeline_overview.md) | End-to-end flow from job discovery to results and post-processing.
| Static simulation flow | [static_simulation_flow.md](static_simulation_flow.md) | Stage-by-stage linear-static workflow (prepare → assemble → BCs → condense → solve → reconstruct → primary/secondary/tertiary results).
| Data flow | [data_flow.md](data_flow.md) | Job input files and result output directories/artifacts.
| Component structure | [component_structure.md](component_structure.md) | Top-level modules (pre_processing, processing, simulation_runner, post_processing) and dependencies.

## Main entry points

- **Job orchestration:** `workflow_orchestrator/run_job.py` — `main()` discovers jobs under `jobs/`, creates result dirs, runs `process_job()` per job (optionally in parallel).
- **Per-job pipeline:** `workflow_orchestrator/run_job.py` — `process_job()` does parsing → element instantiation → K_e/F_e computation → selects static or modal runner → `runner.run()`.
- **Static simulation:** `simulation_runner/static/static_simulation.py` — `StaticSimulationRunner.run()` executes the full linear-static workflow (see [static_simulation_flow.md](static_simulation_flow.md)).
- **Modal simulation:** `simulation_runner/modal/modal_simulation.py` — `ModalSimulationRunner.run()` assembles K/M, applies BCs, solves eigenvalue problem.

Diagrams use [Mermaid](https://mermaid.js.org/) and render in GitHub and most Markdown viewers.
