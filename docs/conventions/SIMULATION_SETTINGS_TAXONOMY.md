# Simulation settings taxonomy (§1–§5)

Job files **`simulation_settings.txt`** map to five canonical analysis families. After [`parse_simulation_settings`](../../pre_processing/parsing/simulation_settings_parser.py), **`settings["type"]`** is always one of: **`static`**, **`eigen`**, **`transient`**, **`harmonic`**, **`buckling`**. Resolution rules live in [`simulation_settings_resolution.py`](../../pre_processing/parsing/simulation_settings_resolution.py) (`finalize_simulation_settings`).

## Sections

| Header | Taxonomy | Role |
|--------|----------|------|
| `[Static]` | §1 | Optional `enabled`; static jobs usually rely on `[Type] static` and `[Newton]` / `[Nonlinear]` for nonlinear behaviour. |
| `[Eigen]` | §2 | `enabled`, `num_modes` (vibration eigenproblem). Optional **`dense_threshold`** (int, default **512**): passed to [`SolveGeneralizedEigenproblem`](../../processing/spectral/operations/solve_generalized_eigenproblem.py) to choose dense vs sparse paths inside [`smallest_generalized_eigenpairs`](../../processing/eigen/smallest_generalized_eigenpairs.py). Optional **`fixed_node_id`** (int): clamp all global DOFs at that node for penalty BCs when combined with prescribed data (see [`processing.boundary_supports`](../../processing/boundary_supports/__init__.py)); if `prescribed_displacement_dict` is absent, same anchor semantics as §3 (node clamp or legacy first six DOFs). |
| `[Transient]` | §3 | `enabled`, `time_step`, `end_time`, `scheme`; merges with legacy `[Dynamic]` when present. Authoritative merged dict: **`effective_transient_config`** in [`simulation_settings_resolution.py`](../../pre_processing/parsing/simulation_settings_resolution.py). Optional BC: **`fixed_node_id`**. Loads: **`load_scale`**, **`load_ramp`**, **`force_time_series_file`** (two columns `time` + scalar scale on `F_ref`, or `time` plus one column per global DOF), **`force_analytic`** (`none` \| `sin` \| `sin_burst`) with **`force_analytic_amplitude`**, **`force_analytic_frequency_hz`**, **`force_analytic_phase_rad`**, **`force_analytic_t_start`**, **`force_analytic_t_end`**. Damping when element `C` is empty: **`rayleigh_alpha`** / **`rayleigh_beta`** (parser aliases **`rayleigh_m`** / **`rayleigh_k`**). |
| `[Harmonic]` | §4 | `enabled` plus frequency sweep / damping keys consumed by [`HarmonicSimulationRunner`](../../simulation_runner/harmonic/harmonic_simulation.py). Optional **`fixed_node_id`** for penalty BCs (same helper as §2/§3; see **Harmonic keys** below). |
| `[Buckling]` | §5 | `enabled`, `num_modes`, `buckling_prestress`, `buckling_load_factor`, `buckling_nonlinear_prestress_twins`. Optional **`nonlinear_buckling`** (bool, default **false**): when **true**, [`run_job`](../../workflow_orchestrator/run_job.py) runs [`NonlinearBucklingSimulationRunner`](../../simulation_runner/buckling/nonlinear_buckling_simulation.py) (MVP stub — not the linearized eigen buckling path); see [`NONLINEAR_BUCKLING_MVP.md`](NONLINEAR_BUCKLING_MVP.md). |

Legacy **`[Modal]`** and **`[Dynamic]`** remain supported; prefer **`[Eigen]`** / **`[Buckling]`** / **`[Transient]`**. Deprecation warnings use **`FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS`** (see [`simulation_runner/README.md`](../../simulation_runner/README.md)).

## Legacy type aliases on the `[Type]` line

| Input | Canonical `type` |
|-------|------------------|
| `modal` + default vibration | `eigen` |
| `modal` + `modal.analysis = buckling` | `buckling` |
| `dynamic` | `transient` |
| `static_nonlinear` | `static` + `_resolved_static_kind = nonlinear` |

## Harmonic keys (§4)

Parsed into **`settings["harmonic"]`** and used by **`HarmonicSimulationRunner`** (direct frequency sweep or optional modal superposition). Optional keys (defaults at runtime if omitted): `frequency_min_hz`, `frequency_max_hz`, `num_frequency_points`, `modal_damping_ratio`, `rayleigh_alpha`, `rayleigh_beta`, `load_phase_rad`, `parallel_frequency_sweep`, `mp_damping_reference`, `use_modal_superposition`, `modal_superposition_num_modes`, `prescribed_motion_phase_rad`, `harmonic_linear_solver`, **`fixed_node_id`** (optional int; penalty clamp at that node, consistent with [`resolve_penalty_fixed_dofs`](../../processing/boundary_supports/__init__.py)). **`[PostProcessing]`** may set **`run_secondary_tertiary_harmonic`**, **`harmonic_frequency_index`**, **`harmonic_secondary_tertiary_all_frequencies`**, **`harmonic_secondary_tertiary_frequency_indices`**, and **`harmonic_secondary_tertiary_displacement_component`** (`real` \| `imag` \| `both`) for formulation-cache stress/strain snapshots from harmonic displacement columns.

Eigen vibration secondary output includes **`{job_name}_modal_generalized_mass.txt`** and, when element nodal loads are available, **`{job_name}_modal_load_participation.txt`** (absolute mass-normalized projection of the assembled reference load onto each mode).

See [`simulation_runner/harmonic/README.md`](../../simulation_runner/harmonic/README.md) and [`HARMONIC_FREQUENCY_DOMAIN.md`](HARMONIC_FREQUENCY_DOMAIN.md).

## Processing packages (§2 assembly)

Canonical global **K** / **M** assembly and penalty boundary conditions for eigen, linear buckling prestress-style matrices, and harmonic analysis live under **`processing.eigen`** (`assembly`, `boundary_conditions`). Legacy **`processing.modal.*`** submodule shims have been removed; import from **`processing.eigen`** (see [**CHANGELOG.md**](CHANGELOG.md) **Removed**).

## Legacy `[Type] modal` / `[Modal]` parsing timeline

The parser still accepts legacy **`[Modal]`** blocks and **`[Type] modal`** (normalized to **`eigen`** / **`buckling`**) and emits **warnings** by default (`FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS=1` silences them). **Removal of this compatibility path is not scheduled before v0.5.0.** When a removal release is planned, [**CHANGELOG.md**](CHANGELOG.md) and this document will record the migration expectation (canonical **`[Eigen]`** / **`[Buckling]`** and **`[Type]`** lines) and any change to the silence flag.

Set **`FEM_LEGACY_MODAL_ERROR=1`** (**`true`** / **`yes`** / **`on`**) to **raise** `ValueError` on legacy **`[Modal]`** or **`[Type] modal`** input (CI stricter than warnings-only).

**Policy summary:** legacy job text stays supported until at least **v0.5.0**; strict env rejects it in CI; **`simulation_runner.modal`** import shims are already **gone** (use **`simulation_runner.spectral`**).

## Terminology: “modal” in docs vs packages

- **Job-file legacy “modal”:** **`[Type] modal`** and **`[Modal]`** sections — deprecated; canonical **`eigen`** / **`buckling`** + **`[Eigen]`** / **`[Buckling]`**.
- **`processing.modal`:** placeholder package only — no **`processing.modal.assembly`** etc.; use **`processing.eigen`** and **`processing.buckling`**.
- **`simulation_runner/spectral/`:** shared **`VibrationBucklingBackend`** for **`EigenSimulationRunner`** / **`LinearBucklingSimulationRunner`** — not the removed **`ModalSimulationRunner`** façade. The **`simulation_runner.modal`** Python shim is **removed**; import from **`simulation_runner.spectral`** ([CHANGELOG](../CHANGELOG.md)).
- **Results keys / folders named “modal”:** e.g. **`run_secondary_tertiary_modal`**, **`primary_results/modal_results/`** — legacy **naming** retained for **on-disk stability** (eigen/buckling primary `.txt`, harmonic **`harmonic_modal_basis_dir`**, manifests). Prefer **`run_secondary_tertiary_eigen`** / **`run_secondary_tertiary_buckling`** in new job prose ([**`RESULTS_DESIGN.md`**](../processing/static/results/RESULTS_DESIGN.md)). Renaming the results folder would be a **semver-major** migration (not done here).

## Deferred / follow-up (not parsed today)

These are **not** read from `simulation_settings.txt` yet; reserved for future product work:

| Topic | Intent |
|-------|--------|
| **`[Eigen] participation_load_scale`** | Optional scalar to scale assembled nodal `F` used only for **`modal_load_participation.txt`** without editing load files. |
| **Transient Rayleigh vs element `C`** | When **both** assembled element damping **`C`** (non-empty) and **`rayleigh_alpha` / `rayleigh_beta`** are set, **`TransientSimulationRunner`** keeps **element `C`** and **ignores** Rayleigh for assembly (warning in logs). Combining models would need an explicit policy and tests. |
