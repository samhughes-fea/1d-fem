# Changelog

All notable changes to this project are documented in this file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project uses [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [Unreleased]

### Removed

- **Breaking:** Public element type strings for thin linear warping aliases are no longer registered: `LinearWarpingEulerBernoulliBeamElement3D`, `LinearWarpingTimoshenkoBeamElement3D`, `LinearWarpingLevinsonBeamElement3D`, `LinearWarpingReddyBeamElement3D`. Use the corresponding baseline type (`LinearEulerBernoulliBeamElement3D`, `LinearTimoshenkoBeamElement3D`, etc.) with `[warping] = 1` in `element.txt` (or `element_dictionary["warping"]`). See [docs/conventions/DEPRECATED_ELEMENT_TYPES.md](conventions/DEPRECATED_ELEMENT_TYPES.md).
- **Breaking:** `NonlinearWarpingTimoshenkoBeamElement3D` is deregistered from `ElementFactory` (implementation was an unimplemented stub). Use `NonlinearTimoshenkoBeamElement3D` for 12-DOF TL Timoshenko until a warping NL Timoshenko is implemented.
- **Breaking:** `NonlinearWarpingEulerBernoulliBeamElement3D` removed — TL EB with Vlasov warping is unified into `NonlinearEulerBernoulliBeamElement3D` with `[warping]` / `element_dictionary["warping"]` (same mesh policy as linear EB).
- **Housekeeping:** Source files for the legacy dedicated classes are deleted: `euler_bernoulli_with_warp/linear_warping_euler_bernoulli_3D.py` and `nonlinear/timoshenko_with_warp/nonlinear_warping_timoshenko_3D.py`. Vlasov operator code remains in `euler_bernoulli_with_warp/utilities/` (used by unified `LinearEulerBernoulliBeamElement3D` / `NonlinearEulerBernoulliBeamElement3D` with warping).

### Added

- [docs/conventions/JOB_INPUT_BEAM_WARPING.md](conventions/JOB_INPUT_BEAM_WARPING.md) — subsection on `prescribed_displacement.txt` with **`CHI`** / **`W`** (7 DOF/node) and global indexing.
- `tests/test_prescribed_displacement_warping.py`, `tests/test_process_job_nl_eb_warp_minimal.py` — parser and minimal end-to-end **`process_job`** coverage for NL EB + `[warping]`.
- [docs/conventions/DEPRECATED_ELEMENT_TYPES.md](conventions/DEPRECATED_ELEMENT_TYPES.md) — manifest of removed/deregistered types and replacements.
- `pre_processing/element_library/removed_element_types.py` — `ensure_element_type_allowed()` raises `ValueError` with migration text and emits `DeprecationWarning` before removal errors.

### Changed

- Modal and dynamic simulation runners set **global** `total_dof` from **`mesh_uses_warping_dof`** (7 vs 6 DOF per node), matching static assembly; dynamic **`process_job`** passes **`element_dictionary`** into **`DynamicSimulationRunner`**.
- **`job_to_abaqus_script`** parses **`prescribed_displacement.txt`** with **`dof_per_node`** aligned to the warping mesh when present.
