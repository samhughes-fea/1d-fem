# `processing.common`

Shared utilities used across analysis families (§1–§5).

## Stage logging

[`stage_logging.init_stage_logger`](stage_logging.py) attaches optional file handlers so each staged operation class can write **`logs/<ClassName>.log`** next to the job tree (typically when `job_results_dir` points under `primary_results/`). Runners that use it include spectral (`PrepareSpectralLocalMatrices`, …), harmonic (`ModifyHarmonicStructuralMatrices`, …), and dynamic (`AssembleDynamicGlobalSystem`, …). See [`simulation_runner/transient/README.md`](../../simulation_runner/transient/README.md) and [`simulation_runner/harmonic/README.md`](../../simulation_runner/harmonic/README.md) for telemetry layout.

Element-level matrix dumps use **`pre_processing.element_library.base_logger_operator.BaseLoggerOperator`** instead (per-element subfolders under the job).

## Telemetry vs diagnostics (quick reference)

| Artifact | Location |
|----------|----------|
| `RuntimeMonitorTelemetry.log` (stage BEGIN/END, resource snapshots) | **`{job_root}/diagnostics/`** — see constructor in [`runtime_monitor_telemetry.py`](../static/diagnostics/runtime_monitor_telemetry.py) |
| `logs/<ClassName>.log` (processing stage classes) | **`{job_root}/logs/`** when stages pass a path under `primary_results/` — see [`stage_logging.py`](stage_logging.py) |

Runner-specific notes: [simulation_runner/README.md](../../simulation_runner/README.md) (Runtime telemetry section).

## Stagewise `run()` vs static `assemble()` / `apply_*()`

Non-static families (**`processing.spectral`**, **`processing.harmonic`**, **`processing.dynamic`**) use a uniform **`run(...)`** method on staged operation classes. Linear static under **`processing.static.operations`** historically uses operation-specific verbs (**`assemble()`**, **`apply_boundary_conditions()`**, etc.). Prefer **`run`** for new non-static stages; do not mass-rename static operations without a dedicated refactor. Details: [`processing/spectral/README.md`](../spectral/README.md).

## Snapshot post-processing (neutral contract)

Non-static families that reuse the static secondary/tertiary **CSV** pipeline do so through a **real displacement snapshot** **`U_global`** plus a **`FormulationResultSet`** built from the job’s **`element_objects`** / **`force_objects`**. The orchestration entry point is **`run_secondary_tertiary_from_formulation_cache`** in `processing.static.results` (see **`RESULTS_DESIGN.md`**). **`processing.common`** hosts small shared helpers such as **`primary_artifact_manifest.write_primary_artifact_manifest`** for machine-readable primary paths; a fuller neutral interface would only move here if multiple families need the same abstraction beyond this manifest.

## Deferred / product backlog

Not implemented here; track if needed:

- **`[Eigen] participation_load_scale`** — scale assembled `F` for participation without editing load files.
- **Transient Rayleigh vs assembled element `C`** — explicit precedence when both are non-empty.
- Expand this README with a **module index** if `processing.common` grows beyond stage logging and small shared helpers.
