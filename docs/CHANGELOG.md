# Changelog

All notable changes to this project are documented in this file. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project uses [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [Unreleased]

### Added

- **2D section solver** ([element_library/solid_section_analysis.md](element_library/solid_section_analysis.md)): `SolidSectionAnalysis` solves the Giavotto (BECAS) warping problem on a 2D section mesh with 3D anisotropic materials. It gives the coupled 6 × 6 stiffness, the shear and tension centres, `E Γ`, mass and second-order terms of any section shape, and `SectionAnalysis.from_solid` passes them to the beam elements. Jobs assign them with `section_mesh_2d.txt` and `element_section.txt`. Exact for homogeneous anisotropic sections under axial force and bending; within 0.2 % of the thin-walled analysis (`hoop="free"`) for the extension, bending and torsion stiffness of a CFRP I-section, and within 0.8 % for the whole 6 × 6 of a CAS box.
- **2D section meshes and 3D materials** ([element_library/solid_section_mesh.md](element_library/solid_section_mesh.md)): `T6` and `Q8` section elements; orthotropic `SolidMaterial` rotated to its ply angle and wall direction (it reduces to the thin-walled ply stiffness under plane stress); meshes from a rectangle, from a thin-walled section (a strip per wall, a row per ply, tied junctions), from Gmsh 4.1 files and from `section_mesh_2d.txt`. The input of the 2D section solver.
- **Beams in any direction, section twist, curved and pretwisted beams** ([element_library/element_frames.md](element_library/element_frames.md)). The linear Euler-Bernoulli, Timoshenko and Benscoter beams take their axis from the node coordinates and an optional `[twist]` column in `element.txt` (degrees about the axis), and their stiffness, mass, loads, geometric stiffness and wall stresses are rotated to global axes. `precurvature.txt` now makes an element the stress-free arc of that curvature and twist between its nodes, built from condensed straight sub-elements, and a twist rate adds the helical-fibre extension-twist coupling to the section (`SectionAnalysis.with_pretwist`). Checked against a turned cantilever, a quarter-circle cantilever, a helical spring and a pretwisted strip.
- **Breaking:** for those three element types `precurvature.txt` no longer adds the equivalent load `Bᵀ D E₀` on a straight element; the curvature is in the geometry. On meshes with rotated or curved elements, point loads must sit on nodes (global components) and distributed and wall loads are refused. Other element types keep the old meaning of `precurvature.txt`.

- **Shell validation of the second-order and corotational beams** ([post_processing/validation_visualisers/abaqus/shell_validation](../post_processing/validation_visualisers/abaqus/shell_validation/README.md)). Five cases written as Abaqus S4R decks from the same section description the beam uses: channel flexural-torsional buckling, I-beam lateral-torsional buckling, the I-beam's prebuckling path near that moment (NLGEOM), large twist of an I-beam with restrained warping (NLGEOM) and a bend-twist coupled CFRP box under large tip deflection (NLGEOM). `run.py` writes the decks, runs Abaqus, extracts the ODB (`extract_odb.py`), reduces each station to beam displacements by a rigid section fit and compares with cached 1d-fem results against a 3 % target. The same decks run in CalculiX (`--solver calculix`, S8R), which checks the buckling cases and, linearly, the nonlinear decks.

### Removed

- **`simulation_runner.modal`** — deprecated import shim (`_vibration_buckling_backend`, `modal_diagnostic`) **deleted**. Use **`simulation_runner.spectral.vibration_buckling_backend.VibrationBucklingBackend`** and **`simulation_runner.spectral.spectral_diagnostics.log_spectral_diagnostics`** instead.
- **`StaticSimulationRunner`** — no longer re-exported from **`simulation_runner.static.linear_static_simulation`** (deprecated module attribute removed). Import **`LinearStaticSimulationRunner`** directly.
- **Breaking:** Public element type strings for thin linear warping aliases are no longer registered: `LinearWarpingEulerBernoulliBeamElement3D`, `LinearWarpingTimoshenkoBeamElement3D`, `LinearWarpingLevinsonBeamElement3D`, `LinearWarpingReddyBeamElement3D`. Use the corresponding baseline type (`LinearEulerBernoulliBeamElement3D`, `LinearTimoshenkoBeamElement3D`, etc.) with `[warping] = 1` in `element.txt` (or `element_dictionary["warping"]`). See [docs/conventions/DEPRECATED_ELEMENT_TYPES.md](conventions/DEPRECATED_ELEMENT_TYPES.md).
- **Breaking:** `NonlinearWarpingTimoshenkoBeamElement3D` is deregistered from `ElementFactory` (implementation was an unimplemented stub). Use `NonlinearTimoshenkoBeamElement3D` for 12-DOF TL Timoshenko until a warping NL Timoshenko is implemented.
- **Breaking:** `NonlinearWarpingEulerBernoulliBeamElement3D` removed — TL EB with Vlasov warping is unified into `NonlinearEulerBernoulliBeamElement3D` with `[warping]` / `element_dictionary["warping"]` (same mesh policy as linear EB).
- **Housekeeping:** Source files for the legacy dedicated classes are deleted: `euler_bernoulli_with_warp/linear_warping_euler_bernoulli_3D.py` and `nonlinear/timoshenko_with_warp/nonlinear_warping_timoshenko_3D.py`. Vlasov operator code remains in `euler_bernoulli_with_warp/utilities/` (used by unified `LinearEulerBernoulliBeamElement3D` / `NonlinearEulerBernoulliBeamElement3D` with warping).

### Changed

- **Fix: Benscoter torsional warping pole.** `ThinWalledSectionAnalysis` constrained the torsional warping against the `y` and `z` rigid modes while solving it. For a section whose origin is off the shear centre (a channel about its web), this put the pole in the wrong place: `E Γ` came out 16 % low. Only the constant mode is constrained now, and the warping is the same for any origin.

- **Breaking: `LinearReddyBeamElement3D` and `LinearLevinsonBeamElement3D` rewritten** as third-order shear elements derived from the section ([element_library/higher_order_shear.md](element_library/higher_order_shear.md)). The previous elements gave a cantilever tip about 30 times too stiff and did not converge, because their rotation shape functions could not represent a linear rotation.
  - Kinematics: the Timoshenko section motion plus the Saint-Venant shear warping slaved to the shear strain, `eps_xx += Psi(s) . chi` with `chi = gamma' - R psi6`. For a solid rectangle this is Reddy's cubic; for thin walls it is the wall analysis' shear warping.
  - Reddy is variational; Levinson keeps only the classical equilibrium equations (Petrov-Galerkin), which gives Timoshenko deflections for uncoupled sections.
  - Nodes carry two extra DOFs, the slopes `v'` and `w'` (`Q1`, `Q2`), so a mesh with these elements has 8 DOFs per node. Fix `Q1`, `Q2` at a clamp to restrain the shear warping (Reddy).
  - Both elements accept thin-walled sections (`element_section.txt`) and give wall and ply stresses. `[warping] = 1` is refused (use `LinearBenscoterBeamElement3D`), and the `[alpha]` column of `section.txt` is no longer used.
  - `SectionAnalysis` gains `third_order_stiffness` (8 x 8) and `third_order_rate`; `ThinWalledSectionAnalysis` computes them from the walls.
  - Removed: the old `levinson/utilities` and `reddy/utilities` operators, their unit tests and the `scripts/*levinson*` verification scripts.

### Added

- **Corotational large-displacement beams on anisotropic sections** ([element_library/corotational_section_beams.md](element_library/corotational_section_beams.md)): `CorotationalEulerBernoulliBeamElement3D`, `CorotationalTimoshenkoBeamElement3D` and `CorotationalBenscoterBeamElement3D`.
  - They wrap the linear elements with their coupled `SectionAnalysis` in a corotational frame (Battini and Pacoste). Nodal rotations are total rotation vectors.
  - The local energy adds the cubic second-order term (offsets, Wagner, bimoment) and the quartic large-twist term from `SectionAnalysis.axial_points`.
  - The residual is exact (complex-step `dp/dU`).
  - Checked against the elastica, the semicircle roll-up, Wagner torsional stiffening, the nonlinear Vlasov equation with restrained warping, and lateral-torsional stability loss with prebuckling deflection.
  - Example `examples/thin_walled/large_twist`.

- **Second-order geometric stiffness, P-delta and linear buckling for the linear beam elements** ([element_library/second_order.md](element_library/second_order.md)):
  - Euler-Bernoulli, Timoshenko, Benscoter, Reddy and Levinson elements build `K_σ` from the energy of a section that is rigid in its own plane, using its current stresses. The section terms `[N, S_y, S_z, W, V_y, V_z, Q_r]` (`ThinWalledSectionAnalysis.geometric_terms`, `SectionAnalysis.geometric_terms`) bring in the shear-centre and tension-centre offsets, the Wagner term and the bimoment for any reference axis.
  - P-delta static solves: `[Static] second_order = p_delta` (`processing/buckling/p_delta.py`).
  - Checked against Euler, Engesser, flexural-torsional (Timoshenko & Gere), lateral-torsional (uniform moment, monosymmetric Wagner, Michell) and beam-column closed forms. Example `examples/thin_walled/channel_buckling`.
  - The previous `K_σ` of the linear EB and Timoshenko elements is replaced; it had no torsional terms and an `N + M/L` factor. The nonlinear elements keep `GeometricStiffnessOperator` until the large-displacement work.
  - The buckling eigensolver no longer symmetrises a nonsymmetric `K` (Levinson).

- **Wall and ply stress recovery for thin-walled sections:**
  - `stress_recovery.py`: full mid-surface strains from the kept beam strains (`N_ss = 0`, free or restrained hoop bending), wall resultants, ply stresses in wall and ply axes, and interlaminar `tau_xn`, `tau_sn` and `sigma_nn` by integrating 3D equilibrium through the thickness. `SaintVenantRecovery` takes section forces.
  - `element_recovery.py`: wall stresses at any point of a solved Euler-Bernoulli, Timoshenko, Benscoter or GBT element.
  - `process_job` writes `secondary_results/wall_stresses/wall_stresses.csv` for linear static jobs with thin-walled sections.
  - Worked examples `examples/thin_walled/box_spar` and `examples/thin_walled/three_cell_blade` with figures.
  - Docs: [element_library/thin_walled_stress_recovery.md](element_library/thin_walled_stress_recovery.md).
- **GBT warping modes, shear lag and wall loads:**
  - The GBT element now uses plane-section beam kinematics with all warping in warping modes: the Saint-Venant warping span (always carried) plus optional shear-lag modes from a new `[shear_modes]` column. A support with the warping amplitudes fixed keeps the section plane, which gives shear lag and the warping-restraint end effect; with them free, tip loads give the Saint-Venant beam.
  - **Breaking (unreleased element):** `LinearGBTBeamElement3D` is now a Timoshenko subclass and its nodal extras are `[b_1 .. b_p, a_1, a_1', ..]` instead of `[psi, a_1, a_1', ..]`. `SectionAnalysis` gains `gbt_sv_modes` and `gbt_shear_modes`; `gbt_matrices(m, p_extra)` uses the new strain order.
  - Wall surface loads from a new optional `wall_load.txt` (`[element_id] [wall_id] [p_x] [p_s] [p_n]`): GBT elements get them on every DOF; Timoshenko and Benscoter elements get the resultant (`ThinWalledSectionAnalysis.wall_load_resultant`).
  - `section_geometry.txt` takes an optional `[Section_Mesh]` block to set the wall subdivision per section.
- **Cross-section modes and the GBT element:**
  - `CrossSectionModes` builds the conventional (Vlasov) distortional modes of any open, closed or multi-cell thin-walled section from its wall frame: a null space of the Vlasov constraints, then a generalised eigenproblem, ordered by transverse bending.
  - It also gives the coupled GBT section stiffness (Benscoter strains plus \(a''\), \(a'\), \(a\) of each mode) and the section mass.
  - New element `LinearGBTBeamElement3D`: the Benscoter element plus Hermite mode amplitudes. It is driven by a new optional `[modes]` column in `element_section.txt`, which also sets the mesh DOF layout. Diaphragms are `Q` supports.
  - `SectionAnalysis` gains `gbt_modes`, `gbt_stiffness`, `gbt_mass` and `gbt_matrices(m)`.
  - Docs: [element_library/gbt_cross_section_modes.md](element_library/gbt_cross_section_modes.md), [proofs/gbt/single_mode_cantilever.md](proofs/gbt/single_mode_cantilever.md).
- **Mesh-wide nodal DOF layout** (`pre_processing/element_library/dof_layout.py`): the nodal stride is 6 plus the largest number of extra DOFs any element row asks for (`[warping] = 1`, or the new optional `element_dictionary["extra_dofs"]`). All runners, `process_job`, the Abaqus export and the modal \(z\) effective-mass metric use it instead of a 6/7 switch. Elements map their own DOF order onto the stride, so elements with fewer DOFs (trusses, plain beams) sit correctly on wider meshes. `prescribed_displacement.txt` accepts `Q1 … Qn` for extra DOFs (`W` / `CHI` stay as aliases of `Q1`). Docs: [conventions/NODAL_DOF_LAYOUT.md](conventions/NODAL_DOF_LAYOUT.md).
- **Benscoter restrained warping for thin-walled sections:**
  - New element `LinearBenscoterBeamElement3D`: Timoshenko kinematics plus an independent warping amplitude \(\psi\) on the seventh nodal DOF (`[warping] = 1`, thin-walled section required).
  - `ThinWalledSectionAnalysis.benscoter_stiffness` (8 × 8: the Saint-Venant 6 × 6 plus warping deviation \(d = \psi - \phi_x\) and \(\psi'\)) and `warping_inertia`; `SectionAnalysis` carries both in element signs.
  - Selective integration of the \(\{\gamma_{xy}, \gamma_{xz}, d\}\) block avoids warping-shear locking.
  - Docs: [element_library/benscoter_restrained_warping.md](element_library/benscoter_restrained_warping.md), [proofs/benscoter/restrained_torsion.md](proofs/benscoter/restrained_torsion.md).
- **Coupled Euler-Bernoulli and Timoshenko elements:**
  - `LinearEulerBernoulliBeamElement3D` and `LinearTimoshenkoBeamElement3D` build `D` from a thin-walled `SectionAnalysis` when one is assigned. Euler-Bernoulli uses the shear-rigid classical 4 × 4; Timoshenko uses the full coupled 6 × 6.
  - The section mass matrix includes mass-centre offsets.
  - `assemble_timoshenko_K0` now under-integrates the full 2 × 2 bending and shear blocks, so the \(D_{34}\) shear-shear coupling is reduced with the shear terms. The result is unchanged for diagonal `D`.
  - Thin-walled sections with the warping DOF are refused for now.
  - Proof: [proofs/timoshenko/coupled_selective_integration.md](proofs/timoshenko/coupled_selective_integration.md).
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

### Fixed

- **Warping meshes (7 DOF per node):** `Element1DBase.assemble_global_dof_indices` now returns indices in the element's own order `[n1 0..5, n2 0..5, n1 W, n2 W]`. Before, it returned them node-major while every 14-DOF beam builds `K_e` with the two warping DOFs last, so `process_job` mixed the warping DOF of node 1 with `u_x` of node 2. Element-level results were unaffected.

### Changed

- **Static global `F` assembly** delegates to **`processing.dynamic.assembly.assemble_global_force_vector`** (supports **`complex128`** when any element `F_e` is complex; optional **`local_global_dof_map`**). Transient and eigen participation reuse the same scatter implementation.
- Modal and dynamic simulation runners set **global** `total_dof` from **`mesh_uses_warping_dof`** (7 vs 6 DOF per node), matching static assembly; transient **`process_job`** passes **`element_dictionary`** into **`TransientSimulationRunner`**.
- **`job_to_abaqus_script`** parses **`prescribed_displacement.txt`** with **`dof_per_node`** aligned to the warping mesh when present.
