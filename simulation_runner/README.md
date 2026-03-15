# Simulation runner

Entry points for running linear/nonlinear static, modal, and dynamic FEM simulations. Jobs are typically launched via [workflow_orchestrator/run_job.py](../workflow_orchestrator/run_job.py), which discovers job directories and invokes the appropriate runner.

## Runners

- **Static**: [static/](static/) — Linear static ([LinearStaticSimulationRunner](static/linear_static_simulation.py)) and nonlinear static ([NonlinearStaticSimulationRunner](static/nonlinear_static_simulation.py)). Assembles from pre_processing elements, solves, computes primary/secondary/tertiary results, and saves outputs.
- **Modal**: [modal/](modal/) — Eigenvalue extraction for natural frequencies and mode shapes.
- **Dynamic**: [dynamic/](dynamic/) — Time integration (when used).

## Usage

From the project root, jobs are run via the workflow orchestrator:

```bash
python workflow_orchestrator/run_job.py [options] [job_dirs...]
```

See [workflow_orchestrator/README.md](../workflow_orchestrator/README.md) and [jobs/README_JOBS.md](../jobs/README_JOBS.md) for job layout and options.
