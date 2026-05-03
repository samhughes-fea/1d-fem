# Transient dynamics (§3)

**Status:** Newmark integration in [`dynamic_simulation.py`](dynamic_simulation.py) via **`TransientSimulationRunner`** (`DynamicSimulationRunner` is a deprecated alias); matrix path in [`processing/dynamic/operations`](../../processing/dynamic/operations/__init__.py). Taxonomy: [SIMULATION_SETTINGS_TAXONOMY.md](../../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md).

## Telemetry

[`RuntimeMonitorTelemetry`](../../processing/static/diagnostics/runtime_monitor_telemetry.py) records high-level stages under the job **`diagnostics/RuntimeMonitorTelemetry.log`** when `job_results_dir` is set (typical `run_job` wiring).

| Stage name | Meaning |
|------------|---------|
| `AssembleDynamicGlobalSystem` | Global `K`, `M`, optional element `C` (or Rayleigh `C` if assembled). |
| `ModifyDynamicGlobalSystem` | Penalty BCs on `K`/`M`/`C`. |
| `IntegrateTransientSystem` | Newmark time stepping. |

Per-class detail logs are written next to the job tree under **`logs/`** (for example `logs/IntegrateTransientSystem.log`) via [`init_stage_logger`](../../processing/common/stage_logging.py).

## `[Transient]` settings (merged with legacy `[Dynamic]`)

Resolved at runtime with [`effective_transient_config`](../../pre_processing/parsing/simulation_settings_resolution.py): time stepping (`time_step`, `end_time`, `scheme`), `load_scale`, `load_ramp`, optional **`fixed_node_id`** (all DOFs at that node when no `prescribed_displacement_dict`), optional **`force_time_series_file`** (two columns `time` + scalar scale on `F_ref`, or `1 + n_dof` columns for full `F(t)`), **`force_analytic`** (`none`, `sin`, `sin_burst`) with `force_analytic_amplitude`, `force_analytic_frequency_hz`, `force_analytic_phase_rad`, and burst window `force_analytic_t_start` / `force_analytic_t_end`. Rayleigh **`rayleigh_alpha`** / **`rayleigh_beta`** (aliases `rayleigh_m` / `rayleigh_k` in the parser) build **`C = α M + β K`** when element damping matrices are absent. Reference nodal loads for `F(t)` are assembled via [`assemble_global_force_vector`](../../processing/dynamic/assembly.py) (no static assembly dependency in the runner).

Relative `force_time_series_file` paths resolve against the job directory (`settings["job_dir"]` from `run_job`).

## Post-processing

**`[PostProcessing]`:** `run_secondary_tertiary_dynamic`, **`dynamic_time_index`** (non-negative index or Python negative wrap) selects the displacement snapshot passed to the static formulation-cache secondary/tertiary path.

## `job_results_dir` layout

| Location | Role |
|----------|------|
| `diagnostics/RuntimeMonitorTelemetry.log` | Stages `AssembleDynamicGlobalSystem`, `ModifyDynamicGlobalSystem`, `IntegrateTransientSystem`. |
| `logs/*.log` | Per-operation logs (`IntegrateTransientSystem`, …). |
| `logs/primary_artifacts.json` | Optional JSON index of primary time-history files. |
| `primary_results/dynamic_results/` | Newmark outputs (see below). |

## Primary time histories (`primary_results/dynamic_results/`)

| File | Shape / meaning |
|------|-----------------|
| `{job_name}_time.txt` | Column vector of time samples \(t_k\) (s), length **n_steps+1**. |
| `{job_name}_displacements.txt` | **U**: rows = time index, columns = global DOF (same order as assembly). |
| `{job_name}_velocities.txt` | **V**: same layout as **U** (first time derivative, consistent units with Newmark). |
| `{job_name}_accelerations.txt` | **A**: same layout (second time derivative). |

Units follow the underlying model (typically SI: m, rad, s). Link snapshots to post-processing via **`dynamic_time_index`** (single row) or **`dynamic_time_indices`** (comma-separated list); see [RESULTS_DESIGN.md](../../processing/static/results/RESULTS_DESIGN.md).
