# Changelog

All notable changes to this project are documented in this file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project uses [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [Unreleased]

### Removed

- **`simulation_runner.modal`** — deprecated import shim (`_vibration_buckling_backend`, `modal_diagnostic`) **deleted**. Use **`simulation_runner.spectral.vibration_buckling_backend.VibrationBucklingBackend`** and **`simulation_runner.spectral.spectral_diagnostics.log_spectral_diagnostics`** instead.
- **`StaticSimulationRunner`** — no longer re-exported from **`simulation_runner.static.linear_static_simulation`** (deprecated module attribute removed). Import **`LinearStaticSimulationRunner`** directly.
- **Breaking:** Public element type strings for thin linear warping aliases are no longer registered: `LinearWarpingEulerBernoulliBeamElement3D`, `LinearWarpingTimoshenkoBeamElement3D`, `LinearWarpingLevinsonBeamElement3D`, `LinearWarpingReddyBeamElement3D`. Use the corresponding baseline type (`LinearEulerBernoulliBeamElement3D`, `LinearTimoshenkoBeamElement3D`, etc.) with `[warping] = 1` in `element.txt` (or `element_dictionary["warping"]`). See [docs/conventions/DEPRECATED_ELEMENT_TYPES.md](conventions/DEPRECATED_ELEMENT_TYPES.md).
- **Breaking:** `NonlinearWarpingTimoshenkoBeamElement3D` is deregistered from `ElementFactory` (implementation was an unimplemented stub). Use `NonlinearTimoshenkoBeamElement3D` for 12-DOF TL Timoshenko until a warping NL Timoshenko is implemented.
- **Breaking:** `NonlinearWarpingEulerBernoulliBeamElement3D` removed — TL EB with Vlasov warping is unified into `NonlinearEulerBernoulliBeamElement3D` with `[warping]` / `element_dictionary["warping"]` (same mesh policy as linear EB).
- **Housekeeping:** Source files for the legacy dedicated classes are deleted: `euler_bernoulli_with_warp/linear_warping_euler_bernoulli_3D.py` and `nonlinear/timoshenko_with_warp/nonlinear_warping_timoshenko_3D.py`. Vlasov operator code remains in `euler_bernoulli_with_warp/utilities/` (used by unified `LinearEulerBernoulliBeamElement3D` / `NonlinearEulerBernoulliBeamElement3D` with warping).

### Added

- **`SectionAnalysis` for every element:**
  - `pre_processing/section_integration/section_analysis.py` holds the 6 × 6 and 4 × 4 section stiffness in the element Voigt order and sign convention (\(\mathbf S = \mathrm{diag}(1,-1,1,1,1,1)\) from the right-handed section analysis), plus \(E\Gamma\), centres and mass.
  - `Element1DBase.section_analysis` gives the thin-walled section assigned in the new optional `element_section.txt`, or the isotropic section from `section.txt` / `material.txt`.
  - `process_job` reads `laminate.txt`, `section_geometry.txt` and `element_section.txt` when present (all three or none).
  - Element types that do not yet use `SectionAnalysis` raise `NotImplementedError` when assigned a thin-walled section.
  - Thin-walled analysis now handles collinear walls: the classical 4 × 4 is available, and the 6 × 6 raises a clear error when a shear direction is rigid.
  - Docs: [element_library/section_analysis.md](element_library/section_analysis.md).
- **Thin-walled section analysis (Saint-Venant level):** `ThinWalledSectionAnalysis` gives the fully coupled 6 × 6 Timoshenko stiffness of open, closed and multi-cell anisotropic thin-walled sections. It uses the Giavotto method on wall mid-lines with Kirchhoff walls, and one cubic warping element per straight wall is exact. It also gives the 4 × 4 classical stiffness, shear and tension centres, torsional warping, \(E\Gamma\) with secondary warping, mass properties and wall strains. Hoop option: `"restrained"` (default) or `"free"`. Also added: `stiffness_about` and `section_analyses_from_files`. Docs: [element_library/thin_walled_section_analysis.md](element_library/thin_walled_section_analysis.md), [proofs/section/saint_venant_thin_walled.md](proofs/section/saint_venant_thin_walled.md).
- **Thin-walled anisotropic section walls** (first step of anisotropic section analysis): `pre_processing/section_integration/thin_walled/` with `PlyMaterial`, `Ply` and `WallLaminate` (classical laminate theory ABD and 5/6 transverse shear `H`), condensed wall laws (`beam_wall_stiffness` 4 × 4, `gbt_wall_stiffness` 5 × 5, general `condensed_wall_stiffness`) and `ThinWalledSection` mid-line geometry with cell counting. New optional job inputs `laminate.txt` and `section_geometry.txt` with parsers `parse_laminate_file` and `parse_section_geometry_file`; they do not affect existing analyses yet. Docs: [element_library/thin_walled_section_walls.md](element_library/thin_walled_section_walls.md), [proofs/section/wall_condensation.md](proofs/section/wall_condensation.md), [conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md](conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md).
- **`TransientSimulationRunner`** (canonical §3 transient orchestration); **`DynamicSimulationRunner`** is a **deprecated** subclass emitting **`DeprecationWarning`** on construction ([`simulation_runner/transient/dynamic_simulation.py`](../simulation_runner/transient/dynamic_simulation.py)).
- **`processing.transient`** — canonical Section 3 processing package; **`processing.dynamic`** now acts as a deprecated compatibility shim while imports migrate.
- **`LinearBucklingSimulationRunner`** (canonical linearized buckling); **`BucklingSimulationRunner`** is a **deprecated** subclass. **`[Buckling] nonlinear_buckling`** (default **false**) selects **`NonlinearBucklingSimulationRunner`** when **true** — MVP diagnostics stub only ([`NONLINEAR_BUCKLING_MVP.md`](conventions/NONLINEAR_BUCKLING_MVP.md)).
- [docs/conventions/SIMULATION_SETTINGS_MODAL_AND_RUNNER_MIGRATION.md](conventions/SIMULATION_SETTINGS_MODAL_AND_RUNNER_MIGRATION.md) — key-alias table, external-repo checklist, deferred **`modal_results/`** rename note.
- `tests/test_nonlinear_buckling_dispatch_mvp.py` — nonlinear buckling stub + `process_job` dispatch wiring.
- **Runner parity:** Section 2–5 runner README tables, spectral pipeline README, **`logs/primary_artifacts.json`** primary index, **`[Eigen] dense_threshold`**, **`[PostProcessing] dynamic_time_indices`**, native **`modal_effective_mass_fraction_z`** (six DOF/node meshes), harmonic parallel sweep **aggregated** per-frequency errors, transient **stability snapshot** log line, family **diagnostic** logs (`transient_run_diagnostic`, `harmonic_run_diagnostic`, `spectral_bc_diagnostic`), and **`run_manifest.paths.primary_artifacts_json`** when present.
- **Simulation settings (taxonomy):** documented **`[Eigen]`** / **`[Transient]`** / **`[Harmonic]`** optional **`fixed_node_id`**; transient **`force_time_series_file`**, **`force_analytic`*** keys, and Rayleigh aliases for assembled **`C`** when element damping is absent. Harmonic structural BCs now honor the same penalty **`fixed_dofs`** resolution as spectral/transient.
- **Eigen secondary:** **`{job_name}_modal_load_participation.txt`** — per-mode scalar \(|\hat{\phi}_j^{\mathsf T} F|\) with mass-normalized \(\hat{\phi}_j\) and \(F\) from element nodal loads (`processing.dynamic.assembly.assemble_global_force_vector`), alongside existing **`modal_generalized_mass`** output.
- [docs/conventions/JOB_INPUT_BEAM_WARPING.md](conventions/JOB_INPUT_BEAM_WARPING.md) — subsection on `prescribed_displacement.txt` with **`CHI`** / **`W`** (7 DOF/node) and global indexing.
- `tests/test_prescribed_displacement_warping.py`, `tests/test_process_job_nl_eb_warp_minimal.py` — parser and minimal end-to-end **`process_job`** coverage for NL EB + `[warping]`.
- [docs/conventions/DEPRECATED_ELEMENT_TYPES.md](conventions/DEPRECATED_ELEMENT_TYPES.md) — manifest of removed/deregistered types and replacements.
- `pre_processing/element_library/removed_element_types.py` — `ensure_element_type_allowed()` raises `ValueError` with migration text and emits `DeprecationWarning` before removal errors.

### Changed

- **Static global `F` assembly** delegates to **`processing.dynamic.assembly.assemble_global_force_vector`** (supports **`complex128`** when any element `F_e` is complex; optional **`local_global_dof_map`**). Transient and eigen participation reuse the same scatter implementation.
- Modal and dynamic simulation runners set **global** `total_dof` from **`mesh_uses_warping_dof`** (7 vs 6 DOF per node), matching static assembly; transient **`process_job`** passes **`element_dictionary`** into **`TransientSimulationRunner`**.
- **`job_to_abaqus_script`** parses **`prescribed_displacement.txt`** with **`dof_per_node`** aligned to the warping mesh when present.
