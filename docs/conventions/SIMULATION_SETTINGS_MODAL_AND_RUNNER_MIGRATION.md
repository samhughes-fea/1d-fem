# Migration: legacy modal input and runner names

This document complements [`SIMULATION_SETTINGS_TAXONOMY.md`](SIMULATION_SETTINGS_TAXONOMY.md) with a **single key-alias table** and **runner import** notes for downstream tooling.

## Key-alias table (job file text)

| Legacy | Canonical |
|--------|-----------|
| `[Type] modal` with `modal.analysis = vibration` (or default) | `[Type] eigen` + `[Eigen]` |
| `[Type] modal` with `modal.analysis = buckling` | `[Type] buckling` + `[Buckling]` |
| `[Modal]` section keys (`num_modes`, `dense_threshold`, …) | `[Eigen]` or `[Buckling]` as appropriate (see [`simulation_settings_parser.py`](../../pre_processing/parsing/simulation_settings_parser.py)) |
| `run_secondary_tertiary_modal` | `run_secondary_tertiary_eigen` / `run_secondary_tertiary_buckling` |

## Python imports (removed shims)

| Removed | Replacement |
|---------|-------------|
| `simulation_runner.modal` | `simulation_runner.spectral` (`VibrationBucklingBackend`, `log_spectral_diagnostics`) |

## Runner renames (canonical names)

| Deprecated (still importable; warns) | Canonical |
|----------------------------------------|-------------|
| `BucklingSimulationRunner` | `LinearBucklingSimulationRunner` |
| `DynamicSimulationRunner` | `TransientSimulationRunner` |

## External repositories (maintainer checklist)

The repository search for `simulation_runner.modal` is **complete for this tree**. For **other** repos, CI, or forks:

1. Grep for `simulation_runner.modal`, `_vibration_buckling_backend`, `log_modal_diagnostics`.
2. Replace with `simulation_runner.spectral` imports as in the table above.
3. Replace `DynamicSimulationRunner` / `BucklingSimulationRunner` with canonical names when updating dependencies.

## Parser and `processing.modal`

- **`[Modal]`** / **`[Type] modal`** remain supported here with deprecation warnings; **removal** is a **semver-major** follow-up after job migration (see taxonomy **Legacy parsing timeline**).
- **`processing.modal`** remains a **doc-only placeholder**; keep until all external imports are confirmed gone, then remove in a dedicated change set.

## `primary_results/modal_results/`

Eigen and **linear** buckling primary `.txt` outputs remain under **`modal_results/`** for **on-disk stability** (harmonic modal basis paths, manifests). Renaming to e.g. `spectral_results/` is **semver-major** and **not** part of the current migration wave.

## Deferred backlog (tracked; not implemented in this wave)

| Item | Status |
|------|--------|
| Remove **`[Modal]`** / **`[Type] modal`** from the parser | **Breaking** — after all known jobs migrate; keep **`FEM_LEGACY_MODAL_ERROR=1`** strict CI until cutover. |
| Rename **`primary_results/modal_results/`** on disk | **Semver-major** — backend, manifests, harmonic eigen-basis I/O, docs/tests. |
| Rename **`transient/dynamic_simulation.py`** → **`transient_simulation.py`** | **Optional** mechanical follow-up; imports update across the tree. |
| Remove **`processing.modal`** placeholder package | Only after repo-wide external import audit; update deprecation tests / CHANGELOG in a dedicated change set. |
