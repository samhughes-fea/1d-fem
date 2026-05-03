# Linear buckling (Section 5)

**Linear buckling runner:** [`LinearBucklingSimulationRunner`](buckling_simulation.py) — uses the same [`VibrationBucklingBackend`](../spectral/vibration_buckling_backend.py) as eigen (Section 2), with prestress and **(K, K_g)** solve stages. **`BucklingSimulationRunner`** is a deprecated alias.

**Nonlinear buckling (MVP):** [`NonlinearBucklingSimulationRunner`](nonlinear_buckling_simulation.py) is selected when **`[Buckling] nonlinear_buckling = true`**; it is a wiring stub only — see [`NONLINEAR_BUCKLING_MVP.md`](../../docs/conventions/NONLINEAR_BUCKLING_MVP.md).

Staged pipeline and **static-vs-spectral mapping**: [simulation_runner/spectral/README.md](../spectral/README.md) (*Stage mapping vs linear static*). Primary outputs use **`primary_results/modal_results/`** (same historical path as eigen; see [SIMULATION_SETTINGS_TAXONOMY.md](../../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md)).

## Primary layout (`job_results_dir`)

| Path | Content |
|------|---------|
| `primary_results/modal_results/{job}_buckling_load_factors.txt` | Buckling load factors λ (smallest positive first). |
| `primary_results/modal_results/{job}_buckling_mode_shapes.txt` | Mode shapes (**n_dof × n_modes**). |
| `logs/primary_artifacts.json` | Machine-readable primary index when present. |

Prestress artefacts (linear or nonlinear) live under nested folders inside **`primary_results/`** as implemented in the backend (see spectral README).

## `[PostProcessing]`

Use **`run_secondary_tertiary_buckling`** (alias of the same internal flag as modal/eigen) with **`buckling_displacement`** (`mode` \| `prestress`), **`buckling_mode_index`**, **`modal_amplitude`**. See [RESULTS_DESIGN.md](../../processing/static/results/RESULTS_DESIGN.md) and [SIMULATION_SETTINGS_TAXONOMY.md](../../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md).
