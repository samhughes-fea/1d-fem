# Job input: beam elements, Vlasov warping, and `section.txt`

This note ties together `ElementParser`, `SectionParser`, and the runtime helpers in `pre_processing/element_library/beam_warping.py`.

## Single policy object (linear and nonlinear beams)

At element construction, linear and Total Lagrangian beam classes resolve the same bundle via
`beam_warping_policy(element_dictionary, row_index, type_string, gamma_from_section_row)` →
`BeamWarpingPolicy` with:

- `mesh_allocates_chi_dof` — global 7 DOF/node mesh if any element enables warping (`mesh_uses_warping_dof`);
- `warping_stiffness_on` — whether **this** row assembles `E·Γ` on `D[6,6]` (`element_warping_stiffness_on`);
- `gamma_section` — Γ from `section_array` (0 if `section.txt` has no `[Gamma]` column);
- `gamma_effective` — Γ used in `D` (zeroed when stiffness is off for the row).

Legacy helpers `mesh_uses_warping_dof`, `element_warping_stiffness_on`, and `effective_warping_gamma` remain; the policy is the recommended reference for documentation and tests.

## `element.txt`: element type and optional columns

- **Standard row:** 11 columns per `ElementParser` (element id, two nodes, type string, six integration orders).
- **Optional 12th column:** either `[curvature]` (legacy scalar, not used by current elements) **or** `[warping]` (0/1 per element), disambiguated by the subheader.
- **Optional 13th column:** when both `[curvature]` and `[warping]` appear, order is **fixed**: curvature then warping.

**Preferred authoring:** use baseline beam type strings (e.g. `LinearTimoshenkoBeamElement3D`, `LinearEulerBernoulliBeamElement3D`, `NonlinearEulerBernoulliBeamElement3D`) and set `[warping]` to **1** where Vlasov warping DOFs and/or stiffness are needed. Removed legacy aliases (`LinearWarping*`, `NonlinearWarpingEulerBernoulliBeamElement3D`) are listed in [`DEPRECATED_ELEMENT_TYPES.md`](DEPRECATED_ELEMENT_TYPES.md). When the `[warping]` column is **absent**, `beam_warping.py` can still infer warping from a type name containing `"Warping"` (legacy); explicit `[warping]` is recommended.

**Nonlinear EB + warping:** use **`NonlinearEulerBernoulliBeamElement3D`** with `[warping]` — the implementation selects 12 vs 14 local DOFs from the same policy as linear EB.

## `section.txt`: column tiers (6 / 8 / 10 / 11)

As documented in `SectionParser`:

| Tier | Trailing columns |
|------|------------------|
| 6 | `[element_id]` … `[J_t]` |
| 8 | … `[kappa]` `[alpha]` |
| 10 | … `[y_sc]` `[z_sc]` |
| 11 | … `[Gamma]` |

Bulk upgrades of older section files can use `scripts/migrate_section_txt_to_11_columns.py`.

### When to set `Gamma`, shear centre, `kappa`, `alpha`

- **`Gamma`:** Use **> 0** when you want Vlasov **warping stiffness** (`E·Γ` on the warping part of the material matrix) **and** warping is enabled for that element (`[warping]` / `element_warping_stiffness_on`). If `[warping]` is off for the row, `effective_warping_gamma` treats Γ as zero for stiffness regardless of the stored section value. Leave **0** for runs without warping stiffness or when you only need mesh/DOF policy without Γ stiffness.

- **`y_sc`, `z_sc`:** Non-zero when the formulation uses shear-centre offsets from `section_array` (e.g. Timoshenko paths). **0** means centroid-centred coupling (no offset).

- **`kappa` / `alpha`:** Shear correction and higher-order coefficients; populate per section physics (defaults from migration scripts may be placeholders).

### Thin-walled / open sections

Default Timoshenko `κ` (often `5/6`) and a scalar `κGA` diagonal are **not** a substitute for full thin-wall shear flow or plate theory. For accuracy-sensitive thin-walled beams, validate against handbook solutions or higher-fidelity models; use non-zero `y_sc`, `z_sc`, and Γ when LTB / warping matters. The beam line still uses resultant **Voigt** `D`; see `FORMULATION_DOCSTRING_STANDARDS.md`.

## `prescribed_displacement.txt` when the mesh uses warping (7 DOF/node)

If any row in `element.txt` turns on **Vlasov warping** (`[warping]` / `element_dictionary["warping"]` or a legacy type name containing `"Warping"`), the global mesh allocates **seven** degrees of freedom per node. The static job runner then parses `prescribed_displacement.txt` with **`dof_per_node = 7`** (see `parse_prescribed_displacement` in `pre_processing/parsing/prescribed_displacement_parser.py`).

- **Local DOF indices 0–5** are the same as for a standard 3D beam: `UX`, `UY`, `UZ`, `RX`, `RY`, `RZ`.
- **Local index 6** is the warping intensity **χ**, accepted in the file as **`CHI`** or **`W`** (alias).

**Global DOF index** for a node `n` and local index `i`:

`global_dof = n * 7 + i` (when the warping mesh is active).

**Boundary conditions:** prescribe **`CHI`** / **`W`** wherever the model needs a warping restraint (e.g. fixed root). For some loadings and element combinations, an additional restraint on **χ** at a free end may be required to remove a null mode in the warping subspace; the reference job `jobs/job_smoke_nl_eb_warp` fixes **χ** at the root and tip. See also `tests/test_warping_integration_benchmark.py` for a minimal element-level discussion.

## Linear modal buckling with warping (7 DOF per node)

When `simulation_settings["type"]` is `modal` and `modal.analysis` is `buckling`, linear EB/Timoshenko elements with `[warping]` allocate **χ** at nodes. The buckling eigenproblem uses elastic **K** and stress geometric stiffness **K_σ** assembled from element `linear_geometric_stiffness_matrix`. For warping-enabled beams, **K_σ** embeds the same **12×12** beam-column geometric stiffness on the standard translation/rotation DOFs as in the nonlinear TL warping tangent (rows/cols for **χ** are zero in **K_σ**). This captures classical beam-column buckling on the line; thin-walled **lateral–torsional** modes that rely on full χ–bending coupling in **K_σ** are only partially represented—validate critical loads against benchmarks when Γ-driven behaviour dominates.

**Workflow alternative:** run buckling with **[warping] off** (or a 6 DOF/node mesh) for a conservative column check, then enable warping for detailed static/post results.

**Regression / benchmarks:** [`tests/test_modal_buckling_warping_benchmark.py`](../../tests/test_modal_buckling_warping_benchmark.py) checks Euler column band (warping-off mesh) and confirms **14×14** `K_σ` equals the **12×12** beam-column block on the first twelve DOFs when χ degrees of freedom are inactive in the displacement vector (embedding behaviour).

## See also

- [`MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md`](MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md) — future nonlinear prestress for linear buckling (design only)
- [`DEPRECATED_ELEMENT_TYPES.md`](DEPRECATED_ELEMENT_TYPES.md) — removed type strings and migrations
- [`BEAM_SHEAR_CORRECTION_AND_THINWALL.md`](BEAM_SHEAR_CORRECTION_AND_THINWALL.md) — Timoshenko \(\kappa\) defaults vs thin-walled sections
- `pre_processing/parsing/element_parser.py`
- `pre_processing/parsing/section_parser.py`
- `pre_processing/element_library/beam_warping.py`
