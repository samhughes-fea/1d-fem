# Simulation settings taxonomy (§1–§5)

Job files **`simulation_settings.txt`** map to five canonical analysis families. After [`parse_simulation_settings`](../../pre_processing/parsing/simulation_settings_parser.py), **`settings["type"]`** is always one of: **`static`**, **`eigen`**, **`transient`**, **`harmonic`**, **`buckling`**. Resolution rules live in [`simulation_settings_resolution.py`](../../pre_processing/parsing/simulation_settings_resolution.py) (`finalize_simulation_settings`).

## Sections

| Header | Taxonomy | Role |
|--------|----------|------|
| `[Static]` | §1 | Optional `enabled`; static jobs usually rely on `[Type] static` and `[Newton]` / `[Nonlinear]` for nonlinear behaviour. |
| `[Eigen]` | §2 | `enabled`, `num_modes` (vibration eigenproblem). |
| `[Transient]` | §3 | `enabled`, `time_step`, `end_time`, `scheme`; merges with legacy `[Dynamic]` when present. |
| `[Harmonic]` | §4 | `enabled` plus frequency sweep / damping keys consumed by [`HarmonicSimulationRunner`](../../simulation_runner/harmonic/harmonic_simulation.py). |
| `[Buckling]` | §5 | `enabled`, `num_modes`, `buckling_prestress`, `buckling_load_factor`, `buckling_nonlinear_prestress_twins`. |

Legacy **`[Modal]`** and **`[Dynamic]`** remain supported; prefer **`[Eigen]`** / **`[Buckling]`** / **`[Transient]`**. Deprecation warnings use **`FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS`** (see [`simulation_runner/README.md`](../../simulation_runner/README.md)).

## Legacy type aliases on the `[Type]` line

| Input | Canonical `type` |
|-------|------------------|
| `modal` + default vibration | `eigen` |
| `modal` + `modal.analysis = buckling` | `buckling` |
| `dynamic` | `transient` |
| `static_nonlinear` | `static` + `_resolved_static_kind = nonlinear` |

## Harmonic keys (§4)

Parsed into **`settings["harmonic"]`** and used by **`HarmonicSimulationRunner`** (direct frequency sweep or optional modal superposition). Optional keys (defaults at runtime if omitted): `frequency_min_hz`, `frequency_max_hz`, `num_frequency_points`, `modal_damping_ratio`, `rayleigh_alpha`, `rayleigh_beta`, `load_phase_rad`, `parallel_frequency_sweep`, `mp_damping_reference`, `use_modal_superposition`, `modal_superposition_num_modes`, `prescribed_motion_phase_rad`, `harmonic_linear_solver`. **`[PostProcessing]`** may set **`run_secondary_tertiary_harmonic`**, **`harmonic_frequency_index`**, **`harmonic_secondary_tertiary_all_frequencies`**, **`harmonic_secondary_tertiary_frequency_indices`**, and **`harmonic_secondary_tertiary_displacement_component`** (`real` \| `imag` \| `both`) for formulation-cache stress/strain snapshots from harmonic displacement columns.

See [`simulation_runner/harmonic/README.md`](../../simulation_runner/harmonic/README.md) and [`HARMONIC_FREQUENCY_DOMAIN.md`](HARMONIC_FREQUENCY_DOMAIN.md).

## Processing packages (§2 assembly)

Canonical global **K** / **M** assembly and penalty boundary conditions for eigen, linear buckling prestress-style matrices, and harmonic analysis live under **`processing.eigen`** (`assembly`, `boundary_conditions`). Legacy **`processing.modal.*`** submodule shims have been removed; import from **`processing.eigen`** (see [**CHANGELOG.md**](CHANGELOG.md) **Removed**).

## Legacy `[Type] modal` / `[Modal]` parsing timeline

The parser still accepts legacy **`[Modal]`** blocks and **`[Type] modal`** (normalized to **`eigen`** / **`buckling`**) and emits **warnings** by default (`FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS=1` silences them). **Removal of this compatibility path is not scheduled before v0.5.0.** When a removal release is planned, [**CHANGELOG.md**](CHANGELOG.md) and this document will record the migration expectation (canonical **`[Eigen]`** / **`[Buckling]`** and **`[Type]`** lines) and any change to the silence flag.

Set **`FEM_LEGACY_MODAL_ERROR=1`** (**`true`** / **`yes`** / **`on`**) to **raise** `ValueError` on legacy **`[Modal]`** or **`[Type] modal`** input (CI stricter than warnings-only).

## Terminology: “modal” in docs vs packages

- **Job-file legacy “modal”:** **`[Type] modal`** and **`[Modal]`** sections — deprecated; canonical **`eigen`** / **`buckling`** + **`[Eigen]`** / **`[Buckling]`**.
- **`processing.modal`:** placeholder package only — no **`processing.modal.assembly`** etc.; use **`processing.eigen`** and **`processing.buckling`**.
- **`simulation_runner/spectral/`:** shared **`VibrationBucklingBackend`** for **`EigenSimulationRunner`** / **`BucklingSimulationRunner`** — not the removed **`ModalSimulationRunner`** façade. Legacy imports under **`simulation_runner.modal`** are deprecated ([CHANGELOG](../CHANGELOG.md)).
- **Results keys / folders named “modal”:** e.g. **`run_secondary_tertiary_modal`**, **`modal_results/`** — legacy naming; prefer **`run_secondary_tertiary_eigen`** / **`run_secondary_tertiary_buckling`** aliases in new prose ([**`RESULTS_DESIGN.md`**](../processing/static/results/RESULTS_DESIGN.md)).
