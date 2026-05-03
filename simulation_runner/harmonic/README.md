# Harmonic analysis (§4)

**Status:** Direct frequency sweep implemented in [`harmonic_simulation.py`](harmonic_simulation.py); kernels in [`processing/harmonic/frequency_response.py`](../../processing/harmonic/frequency_response.py). Design reference: [HARMONIC_FREQUENCY_DOMAIN.md](../../docs/conventions/HARMONIC_FREQUENCY_DOMAIN.md).

## Job file schema

Parser accepts optional **`[Harmonic]`** keys alongside **`[Type] harmonic`**:

| Key | Meaning |
|-----|---------|
| `enabled` | Must be `true` when using taxonomy activation with `[Harmonic]` (see resolution rules). |
| `frequency_min_hz` | Lower sweep bound (Hz); alias `frequency_min`. |
| `frequency_max_hz` | Upper sweep bound (Hz); alias `frequency_max`. |
| `num_frequency_points` | Number of frequency samples; alias `num_points`. |
| `modal_damping_ratio` | Scalar \(\zeta\) for mass-proportional damping; alias `damping_ratio`. |
| `rayleigh_alpha` | Rayleigh \(\alpha\) on **M** (1/s); alias `rayleigh_m`. |
| `rayleigh_beta` | Rayleigh \(\beta\) on **K** (s); alias `rayleigh_k`. |
| `load_phase_rad` | Global phase (rad) on the assembled load vector; alias `load_phase`. |
| `parallel_frequency_sweep` | `true` to solve frequency samples concurrently. |
| `mp_damping_reference` | `geometric_mean` (default) or `current_frequency` — mass-proportional reference \(\omega\) for **C**. |
| `use_modal_superposition` | `true` for undamped-mode expansion + diagonal modal damping (no Rayleigh on this path). |
| `modal_superposition_num_modes` | Number of modes for the modal path. |
| `prescribed_motion_phase_rad` | Phase (rad) on prescribed-motion DOFs (alias `prescribed_motion_phase`). |
| `harmonic_linear_solver` | `spsolve` (default) or `splu`. |
| `fixed_node_id` | Optional int: clamp all global DOFs at that node for penalty BCs (same [`resolve_penalty_fixed_dofs`](../../processing/boundary_supports/__init__.py) convention as §2/§3). |

If omitted, defaults are applied at runtime (see `effective_harmonic_config` in [`harmonic_simulation.py`](harmonic_simulation.py)).

**`[PostProcessing]`:** `run_secondary_tertiary_harmonic`, `harmonic_frequency_index`, optional **`harmonic_secondary_tertiary_all_frequencies`** — same pattern as modal/dynamic (real part of \(\mathbf{u}\)); see [`HARMONIC_FREQUENCY_DOMAIN.md`](../../docs/conventions/HARMONIC_FREQUENCY_DOMAIN.md).

These populate **`simulation_settings["harmonic"]`** after [`parse_simulation_settings`](../../pre_processing/parsing/simulation_settings_parser.py). See also [SIMULATION_SETTINGS_TAXONOMY.md](../../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md).

### Telemetry and per-stage logs

When `job_results_dir` is set, [`RuntimeMonitorTelemetry`](../../processing/static/diagnostics/runtime_monitor_telemetry.py) writes stage boundaries to **`diagnostics/RuntimeMonitorTelemetry.log`**. Fine-grained logs for harmonic processing stages live under **`logs/`** next to the job tree (same pattern as §3 transient — see [`processing/common/stage_logging`](../../processing/common/stage_logging.py)).

### Outputs and `job_results_dir` tree

| Location | Role |
|----------|------|
| `diagnostics/RuntimeMonitorTelemetry.log` | High-level stage boundaries (assemble → modify → damping → sweep). |
| `logs/*.log` | Per-stage class logs (`ModifyHarmonicStructuralMatrices`, `SolveHarmonicFrequencySweep`, …). |
| `logs/primary_artifacts.json` | Optional JSON index of primary text outputs (schema `schema_version` **1.0**). |
| `primary_results/assembly_modal/`, `assembly_force/` | Intermediate assembly dumps when stages log there. |
| `primary_results/harmonic_results/` | **Primary sweep outputs** (see below). |

Primary matrices under **`primary_results/harmonic_results/`**: **`{job_name}_frequencies_hz.txt`**, **`{job_name}_displacement_real.txt`**, **`{job_name}_displacement_imag.txt`**, **`{job_name}_displacement_abs.txt`**, **`{job_name}_displacement_phase_rad.txt`** — each displacement file is **global DOFs × frequency index** (column `k` is sample `f_k`). Schema and units: [`HARMONIC_FREQUENCY_DOMAIN.md`](../../docs/conventions/HARMONIC_FREQUENCY_DOMAIN.md).

Roadmap for Section 4 extensions: same doc (*Primary output schema* and tertiary backlog).
