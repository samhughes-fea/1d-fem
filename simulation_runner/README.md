# Simulation runner

Entry points for running FEM simulations aligned with the analysis taxonomy **§1–§5**. Jobs are typically launched via [workflow_orchestrator/run_job.py](../workflow_orchestrator/run_job.py), which discovers job directories and invokes the appropriate runner.

## Canonical simulation `type` (dispatch)

After parsing `simulation_settings.txt`, `settings["type"]` is always one of:

| `type` | Taxonomy | Package / runner |
|--------|----------|------------------|
| `static` | §1 | [static/](static/) — linear ([LinearStaticSimulationRunner](static/linear_static_simulation.py)) or nonlinear ([NonlinearStaticSimulationRunner](static/nonlinear_static_simulation.py)); nonlinear jobs set `_resolved_static_kind == "nonlinear"` (legacy `[Type] static_nonlinear` normalizes here). |
| `eigen` | §2 | [eigen/](eigen/) — `EigenSimulationRunner` (natural frequencies / mode shapes). Legacy `[Type] modal` with `modal.analysis=vibration` normalizes to `eigen`. |
| `transient` | §3 | [transient/](transient/) — `TransientSimulationRunner` ([transient/dynamic_simulation.py](transient/dynamic_simulation.py)); `DynamicSimulationRunner` remains a deprecated subclass. Legacy `[Type] dynamic` normalizes to `transient`. |
| `harmonic` | §4 | [harmonic/](harmonic/) — `HarmonicSimulationRunner` (frequency sweep); see [harmonic/README.md](harmonic/README.md). |
| `buckling` | §5 | [buckling/](buckling/) — `LinearBucklingSimulationRunner` (linearized **K** / **K_g**); optional **`[Buckling] nonlinear_buckling = true`** dispatches **`NonlinearBucklingSimulationRunner`** (MVP stub). `BucklingSimulationRunner` is a deprecated alias. Legacy `[Type] modal` with `modal.analysis=buckling` normalizes to `buckling`. |

Optional bracket sections **`[Static]`**, **`[Eigen]`**, **`[Transient]`**, **`[Harmonic]`**, **`[Buckling]`** support `enabled = true` to declare the primary analysis; see [pre_processing/parsing/simulation_settings_parser.py](../pre_processing/parsing/simulation_settings_parser.py) and `finalize_simulation_settings` in [simulation_settings_resolution.py](../pre_processing/parsing/simulation_settings_resolution.py).

## Deprecations (input files)

- **`[Type] modal`** — still parsed; resolves to **`eigen`** or **`buckling`** via `modal.analysis`. Prefer **`[Type] eigen`** or **`buckling`** with **`[Eigen]`** / **`[Buckling]`** sections.
- **`[Modal]`** — deprecated; use **`[Eigen]`** (vibration) or **`[Buckling]`** (linear buckling). Legacy merge-on-read is unchanged.
- **`[Type] dynamic`** — alias for **`transient`**; prefer **`[Type] transient`** and optional **`[Transient]`**.
- **`[Type] static_nonlinear`** — alias for **`static`** with **`_resolved_static_kind == "nonlinear"`**; prefer **`[Type] static`** plus **`[Newton]`** / **`[Nonlinear]`**.

Parser deprecation warnings log at **`WARNING`** unless **`FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS=1`** is set.

Set **`FEM_LEGACY_MODAL_ERROR=1`** to **raise** on legacy **`[Modal]`** / **`[Type] modal`** input (see [**SIMULATION_SETTINGS_TAXONOMY.md**](../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md)).

**`ModalSimulationRunner`** was removed — use **`EigenSimulationRunner`** / **`LinearBucklingSimulationRunner`** via **`run_job`** ([**CHANGELOG**](../docs/CHANGELOG.md) **Removed**).

## Shared spectral backend (eigen + buckling)

- **[spectral/](spectral/)** — **`VibrationBucklingBackend`** ([`vibration_buckling_backend.py`](spectral/vibration_buckling_backend.py)) shared by **`EigenSimulationRunner`** and **`LinearBucklingSimulationRunner`** (not a separate public façade). Import **`simulation_runner.spectral`** (or concrete submodules); the old **`simulation_runner.modal`** shim package has been **removed** ([CHANGELOG](../docs/CHANGELOG.md)).

## Processing packages (§2 assembly)

- **`processing.eigen`** — canonical global **K** / **M** assembly and penalty BCs (`assembly`, `boundary_conditions`) used by eigen/buckling/harmonic runners.
- **`processing.modal`** — placeholder package only (**`processing/modal/__init__.py`**); legacy **`processing.modal.assembly`**, **`boundary_conditions`**, and **`buckling`** submodules were **removed** — import **`processing.eigen`** and **`processing.buckling`** instead ([CHANGELOG](../docs/CHANGELOG.md) **Removed**).
- **`processing.buckling`** — linear buckling kernels; reuses §2 scatter/BC helpers from **`processing.eigen`**. Runners import **`processing.buckling`** directly.

**Transient:** import **`TransientSimulationRunner`** from **`simulation_runner.transient.dynamic_simulation`** (the old **`simulation_runner.dynamic`** shim was removed; **`DynamicSimulationRunner`** is a deprecated alias).

## Runtime telemetry and per-stage logs

[`RuntimeMonitorTelemetry`](../processing/static/diagnostics/runtime_monitor_telemetry.py) always writes **`{job_results_root}/diagnostics/RuntimeMonitorTelemetry.log`**: the constructor takes any **descendant** path of the job root (for example `primary_results/` or `diagnostics/`) and normalizes to the parent job folder, then appends `diagnostics/`.

| Runner | Typical `RuntimeMonitorTelemetry(...)` argument |
|--------|-----------------------------------------------|
| §2 eigen / §5 buckling (spectral) | `diagnostics_dir` under job root |
| §3 transient (`TransientSimulationRunner`) | `diagnostics_dir` |
| §4 harmonic | `diagnostics_dir` |
| §1 linear static | `primary_results_dir` in `__init__`, then `diagnostics_dir` in `solve_linear_system_only` / `run()` — effective log path is still **job** `diagnostics/` |
| §1 nonlinear static | `primary_results_dir` on the long-lived monitor, `diagnostics_dir` on the top-level `run()` monitor |

Per-class **`logs/<StageClassName>.log`** files come from [`init_stage_logger`](../processing/common/stage_logging.py) and use **`parent(job_results_dir) / "logs"`** where `job_results_dir` is usually **`primary_results/...`** for that stage — so operators look under **`{job_results_root}/logs/`** for `AssembleGlobalSystem.log`, `PrepareLocalSystem.log`, etc.

## Primary outputs and `job_results_dir` layout (Sections 2–5)

Canonical job root is **`job_results_dir`** (same tree for all types). Below: **primary** paths only; optional formulation-cache post writes **`secondary_results/`** and **`tertiary_results/`** when **`[PostProcessing]`** flags are set (see [RESULTS_DESIGN.md](../processing/static/results/RESULTS_DESIGN.md)).

| `type` | Primary folders / files | Post keys (see [SIMULATION_SETTINGS_TAXONOMY.md](../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md)) |
|--------|-------------------------|-----------------------------------------------------------------------------------------------------|
| `eigen` | `primary_results/modal_results/{job}_frequencies.txt`, `{job}_mode_shapes.txt`; optional secondary `.txt` in same folder | `run_secondary_tertiary_eigen` / `run_secondary_tertiary_modal`, `modal_mode_index`, `modal_amplitude` |
| `buckling` | `primary_results/modal_results/{job}_buckling_load_factors.txt`, `{job}_buckling_mode_shapes.txt` | Same modal flag family + `buckling_displacement`, `buckling_mode_index` |
| `transient` | `primary_results/dynamic_results/{job}_time.txt`, `_displacements.txt`, `_velocities.txt`, `_accelerations.txt` | `run_secondary_tertiary_dynamic`, `dynamic_time_index`, **`dynamic_time_indices`** (comma-separated list; multiple snapshots use `secondary_results/dynamic_post/t_*`) |
| `harmonic` | `primary_results/harmonic_results/` — frequencies, displacement real/imag/abs/phase matrices | `run_secondary_tertiary_harmonic`, `harmonic_frequency_index`, multi-frequency keys |

**Machine-readable index:** successful runs may write **`logs/primary_artifacts.json`** (schema in [workflow_orchestrator/run_manifest.py](../workflow_orchestrator/run_manifest.py) companion field **`paths.primary_artifacts_json`** inside **`logs/run_manifest.json`** when that file exists).

Family READMEs: [eigen/README.md](eigen/README.md), [buckling/README.md](buckling/README.md), [spectral/README.md](spectral/README.md), [harmonic/README.md](harmonic/README.md), [transient/README.md](transient/README.md).

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

Nodal **force** scatter is shared via [`assemble_global_force_vector`](../processing/dynamic/assembly.py) (used by static `AssembleGlobalSystem`, transient, and eigen participation). **Stiffness** scatter between [processing/eigen/assembly.py](../processing/eigen/assembly.py) and static assembly remains intentionally separate until dispatch and regressions are stable.
