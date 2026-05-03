# Eigen / natural vibration (Section 2)

**Runner:** [`EigenSimulationRunner`](eigen_simulation.py) — thin façade over [`VibrationBucklingBackend`](../spectral/vibration_buckling_backend.py) for **`[Type] eigen`**.

## Primary layout (`job_results_dir`)

| Path | Content |
|------|---------|
| `primary_results/modal_results/{job}_frequencies.txt` | Natural frequencies (Hz), one per mode. |
| `primary_results/modal_results/{job}_mode_shapes.txt` | Mode matrix **n_dof × n_modes** (columns are modes). |
| `primary_results/modal_results/{job}_modal_generalized_mass.txt` | Optional secondary when post-processing snapshot is **off** — see [RESULTS_DESIGN.md](../../processing/static/results/RESULTS_DESIGN.md). |
| `primary_results/modal_results/{job}_modal_load_participation.txt` | Optional — mass-normalized load projection. |
| `primary_results/modal_results/{job}_modal_effective_mass_fraction_z.txt` | Optional — directional effective mass fraction for a **global +Z** unit pattern (6 DOF/node meshes only; skipped for 7 DOF/warping). |
| `logs/primary_artifacts.json` | Machine-readable index of primary files (when the runner finishes successfully). |

## `[PostProcessing]`

Aliases **`run_secondary_tertiary_eigen`** / legacy **`run_secondary_tertiary_modal`** drive formulation-cache secondary/tertiary from a mode snapshot. Keys: **`modal_mode_index`**, **`modal_amplitude`**. Taxonomy: [SIMULATION_SETTINGS_TAXONOMY.md](../../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md); behaviour: [RESULTS_DESIGN.md](../../processing/static/results/RESULTS_DESIGN.md).

## Shared implementation

See [simulation_runner/spectral/README.md](../spectral/README.md) for the full staged pipeline diagram.
