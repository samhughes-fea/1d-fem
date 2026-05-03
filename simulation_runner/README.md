# Simulation runner

Entry points for running FEM simulations aligned with the analysis taxonomy **§1–§5**. Jobs are typically launched via [workflow_orchestrator/run_job.py](../workflow_orchestrator/run_job.py), which discovers job directories and invokes the appropriate runner.

## Canonical simulation `type` (dispatch)

After parsing `simulation_settings.txt`, `settings["type"]` is always one of:

| `type` | Taxonomy | Package / runner |
|--------|----------|------------------|
| `static` | §1 | [static/](static/) — linear ([LinearStaticSimulationRunner](static/linear_static_simulation.py)) or nonlinear ([NonlinearStaticSimulationRunner](static/nonlinear_static_simulation.py)); nonlinear jobs set `_resolved_static_kind == "nonlinear"` (legacy `[Type] static_nonlinear` normalizes here). |
| `eigen` | §2 | [eigen/](eigen/) — `EigenSimulationRunner` (natural frequencies / mode shapes). Legacy `[Type] modal` with `modal.analysis=vibration` normalizes to `eigen`. |
| `transient` | §3 | [transient/](transient/) — `DynamicSimulationRunner` ([transient/dynamic_simulation.py](transient/dynamic_simulation.py)). Legacy `[Type] dynamic` normalizes to `transient`. |
| `harmonic` | §4 | [harmonic/](harmonic/) — `HarmonicSimulationRunner` (frequency sweep); see [harmonic/README.md](harmonic/README.md). |
| `buckling` | §5 | [buckling/](buckling/) — `BucklingSimulationRunner`. Legacy `[Type] modal` with `modal.analysis=buckling` normalizes to `buckling`. |

Optional bracket sections **`[Static]`**, **`[Eigen]`**, **`[Transient]`**, **`[Harmonic]`**, **`[Buckling]`** support `enabled = true` to declare the primary analysis; see [pre_processing/parsing/simulation_settings_parser.py](../pre_processing/parsing/simulation_settings_parser.py) and `finalize_simulation_settings` in [simulation_settings_resolution.py](../pre_processing/parsing/simulation_settings_resolution.py).

## Deprecations (input files)

- **`[Type] modal`** — still parsed; resolves to **`eigen`** or **`buckling`** via `modal.analysis`. Prefer **`[Type] eigen`** or **`buckling`** with **`[Eigen]`** / **`[Buckling]`** sections.
- **`[Modal]`** — deprecated; use **`[Eigen]`** (vibration) or **`[Buckling]`** (linear buckling). Legacy merge-on-read is unchanged.
- **`[Type] dynamic`** — alias for **`transient`**; prefer **`[Type] transient`** and optional **`[Transient]`**.
- **`[Type] static_nonlinear`** — alias for **`static`** with **`_resolved_static_kind == "nonlinear"`**; prefer **`[Type] static`** plus **`[Newton]`** / **`[Nonlinear]`**.

Parser deprecation warnings log at **`WARNING`** unless **`FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS=1`** is set.

Set **`FEM_LEGACY_MODAL_ERROR=1`** to **raise** on legacy **`[Modal]`** / **`[Type] modal`** input (see [**SIMULATION_SETTINGS_TAXONOMY.md**](../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md)).

**`ModalSimulationRunner`** was removed — use **`EigenSimulationRunner`** / **`BucklingSimulationRunner`** via **`run_job`** ([**CHANGELOG**](../docs/CHANGELOG.md) **Removed**).

## Shared spectral backend (eigen + buckling)

- **[spectral/](spectral/)** — **`VibrationBucklingBackend`** ([`vibration_buckling_backend.py`](spectral/vibration_buckling_backend.py)) shared by **`EigenSimulationRunner`** and **`BucklingSimulationRunner`** (not a separate public façade).
- **`simulation_runner.modal`** — **deprecated** one-release import path; **`_vibration_buckling_backend`** / **`modal_diagnostic`** emit **`DeprecationWarning`** and delegate to **`simulation_runner.spectral`** ([CHANGELOG](../docs/CHANGELOG.md)).

## Processing packages (§2 assembly)

- **`processing.eigen`** — canonical global **K** / **M** assembly and penalty BCs (`assembly`, `boundary_conditions`) used by eigen/buckling/harmonic runners.
- **`processing.modal`** — placeholder package only (**`processing/modal/__init__.py`**); legacy **`processing.modal.assembly`**, **`boundary_conditions`**, and **`buckling`** submodules were **removed** — import **`processing.eigen`** and **`processing.buckling`** instead ([CHANGELOG](../docs/CHANGELOG.md) **Removed**).
- **`processing.buckling`** — linear buckling kernels; reuses §2 scatter/BC helpers from **`processing.eigen`**. Runners import **`processing.buckling`** directly.

**Transient:** import **`DynamicSimulationRunner`** from **`simulation_runner.transient.dynamic_simulation`** (the old **`simulation_runner.dynamic`** shim was removed).

## Usage

From the project root, jobs are run via the workflow orchestrator:

```bash
python workflow_orchestrator/run_job.py [options] [job_dirs...]
```

See [workflow_orchestrator/README.md](../workflow_orchestrator/README.md) and [jobs/README_JOBS.md](../jobs/README_JOBS.md) for job layout and options.

## Template for new analysis domains

Mirror the static stack when adding a new primary analysis type:

- Under **`processing/<domain>/`**, use **`operations/`** (assembly, solve helpers), **`results/`** (primary → secondary → tertiary as needed), and **`diagnostics/`** (optional).
- Add a thin orchestration class under **`simulation_runner/<domain>/`** (same role as [static/linear_static_simulation.py](static/linear_static_simulation.py) / [static/nonlinear_static_simulation.py](static/nonlinear_static_simulation.py)): wire parsing outputs to processing stages and I/O only.

## Deferred refactor (assembly sharing)

Extracting shared scatter/assembly utilities between [processing/eigen/assembly.py](../processing/eigen/assembly.py) and [processing/static/operations/assembly.py](../processing/static/operations/assembly.py) is intentionally deferred until static dispatch and job regressions are stable; track as a separate change-set with eigen/static end-to-end checks.
