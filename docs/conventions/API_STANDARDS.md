# API standards — 1D FEM (beams / frames)

Applied conventions for this repository: **weak-form assembly**, **tensor layout**, **linear `utilities/` modules**, **nonlinear TL extensions**, and minimal Python style rules.

**Element and utility docstrings:** Use **NumPy / numpydoc** section layout and **plain-text, code-like math** in Python source (e.g. ``K_e += B.T @ D @ B * w_g * detJ``), not LaTeX-only exposition. Full checklist, Voigt tables, TL utility roles, and copy-paste weak-form lines: [FORMULATION_DOCSTRING_STANDARDS.md](FORMULATION_DOCSTRING_STANDARDS.md) (*Docstring gold standard*, *Plain-text math*, *Nonlinear Total Lagrangian*).

---

## 1. Stiffness: B, D, J

At each Gauss point (natural coordinate ξ, physical `x` along the element):

- **Jacobian:** `detJ = dx/dξ` (e.g. `L/2` for a straight 2-node map). Element code exposes this as `jacobian` / `jacobian_determinant`.
- **Strain–displacement:** `ε = B U_e` with `B` in **physical** `x` where required. [`StrainDisplacementOperator`](../../pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/utilities/B_matrix.py) provides `natural_coordinate_form` (B̃) and `physical_coordinate_form` (B).
- **Constitutive:** Stress resultants from `S = D @ ε` with `D` from [`MaterialStiffnessOperator`](../../pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/utilities/D_matrix.py). Row ordering of `ε` and `S` must match (e.g. Timoshenko-style six components).
- **Quadrature:** `K_e += Bᵀ D B w_g detJ` (sum over Gauss points). **Selective integration** (e.g. Timoshenko shear) is still the same pattern with different rules on **sub-rows** of `B`.

**Generality:** Warping or higher-order theories may use **7×14** (or other) `B`; document shapes in the element docstring — do not assume `(6, 12)` only in prose.

---

## 2. Loads: N and f

- **Distributed:** `F_dist = ∫ Nᵀ q detJ dξ` with `q` interpolated at Gauss `x` ([`LoadInterpolationOperator`](../../pre_processing/element_library/linear/beam/first_order_shear_deformation_theory/timoshenko/utilities/interpolate_loads.py)).
- **Point:** `F_point = N(x_p)ᵀ P`.
- Shape functions `N`, `dN_dξ`, `d2N_dξ2` come from [`get_shape_function_operator`](../../pre_processing/element_library/shape_function_registry.py) or theory-local [`shape_functions.py`](../../pre_processing/element_library/linear/beam/zero_order_shear_deformation_theory/euler_bernoulli/utilities/shape_functions.py).

---

## 3. Linear elements: `utilities/` layout

Under each beam-theory package, e.g. `pre_processing/element_library/linear/beam/<sdft_order>/<theory>/utilities/` (and similarly for bar/truss under `linear/bar/`, `linear/truss/`):

| File | Role |
|------|------|
| `B_matrix.py` | Strain–displacement operator(s) |
| `D_matrix.py` | Section/material stiffness |
| `interpolate_loads.py` | Distributed loads at Gauss points |
| `shape_functions.py` | `N(ξ)` and derivatives |

**Re-exports:** If logic is shared (e.g. warping Euler–Bernoulli reusing the straight EB `D_matrix`), the file must still exist and the module docstring must state the dependency.

**Optional extras:** e.g. [`local_frame.py`](../../pre_processing/element_library/linear/bar/utilities/local_frame.py) for bar/truss — not counted against the four-file rule.

---

## 4. Nonlinear (geometric, Total Lagrangian)

**Composition:** Nonlinear elements **import** linear **`D_matrix`**, **`B_matrix`** / shape derivatives, **load** utilities, and registry where applicable.

**Weak-form Gauss rule (mandatory):** `F_int`, `K_mat`, and `K_σ` must be assembled as **sums over Gauss points** with weights `w_g` and `detJ`, consistent with [FORMULATION_DOCSTRING_STANDARDS.md](FORMULATION_DOCSTRING_STANDARDS.md) § *Weak-form assembly*. In particular, **`K_σ`** must not rely on undocumented closed-form matrices; it must match a stated integral (e.g. stress-weighted products of bending **slopes** `∂h/∂x` summed with `w_g |J|`).

**Nonlinear `utilities/` (typical):**

| Module | Role |
|--------|------|
| `green_lagrange_strain.py` | `E(u)`, `B_lin`, `B_nl` |
| `stress_resultant.py` | `S = D @ E`, extract `N, M_y, M_z` |
| `geometric_stiffness.py` | `K_σ` via **Gauss sum** of section forces with `dN/dx` (bending planes) |

Other families (e.g. GEBT) may use a **different** utility set; document per family.

---

## 5. Simulation runners and jobs

- **Entry point:** [`run_job.py`](../../workflow_orchestrator/run_job.py) builds elements and dispatches by **canonical** `simulation_settings["type"]`: **`static`**, **`eigen`**, **`transient`**, **`harmonic`**, **`buckling`**. Legacy **`[Type]`** lines **`modal`**, **`dynamic`**, **`static_nonlinear`** still parse and normalize (see [`SIMULATION_SETTINGS_TAXONOMY.md`](SIMULATION_SETTINGS_TAXONOMY.md)).
- **§2 global matrices:** Use **`processing.eigen.assembly`** / **`processing.eigen.boundary_conditions`** for **`assemble_global_matrices`** / **`apply_boundary_conditions`**. The legacy submodule paths **`processing.modal.assembly`**, **`processing.modal.boundary_conditions`**, and **`processing.modal.buckling`** have been **removed** (see [**CHANGELOG.md**](../CHANGELOG.md) **Removed**); **`processing.modal`** is a doc-only placeholder package. Buckling kernels live under **`processing.buckling`**. Smallest-modal (**K**, **M**) pencils for §2-style eigen bases use **`processing.eigen.smallest_generalized_eigenpairs`** (shared with harmonic modal superposition and **`VibrationBucklingBackend`**).
- **Linear buckling (§5):** When `type` is **`buckling`** (legacy: `modal` + `modal.analysis = buckling`), [`BucklingSimulationRunner`](../../simulation_runner/buckling/buckling_simulation.py) / shared modal pipeline builds prestress displacements via **`linear_static`** (default) or **`nonlinear_static`** (`buckling.buckling_prestress` or legacy `modal.buckling_prestress`), then assembles elastic **`K`** and stress-based geometric stiffness **`K_σ`**, then solves the linear buckling eigenproblem. **`nonlinear_static`** requires nonlinear beam elements and a converged nonlinear prestress solve; see [`MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md`](MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md). **Warping meshes (7 DOF per node)** are supported for linear EB/Timoshenko: **`K_σ`** embeds the same **12×12** beam-column geometric stiffness on standard translation/rotation DOFs as in the nonlinear TL warping path (χ rows/cols do not receive extra geometric stiffness terms). See [`JOB_INPUT_BEAM_WARPING.md`](JOB_INPUT_BEAM_WARPING.md) (*Linear buckling with warping*).
- **Co-rotational tangent:** Optional `[Nonlinear]` / `nonlinear` dict key **`corotational_tangent_mode`**: `finite_difference` (default, consistent tangent via force probing) or `elastic_material` (analytic `Tᵀ K_local T` without spin stiffness — see [`large_rotation_vs_total_lagrangian.md`](../element_library/large_rotation_vs_total_lagrangian.md)). Passed only when instantiating [`CorotationalBeamElement3D`](../../pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_3D.py). End-to-end Newton parity for both modes is smoke-tested in [`tests/test_corotational_nonlinear_job_tangent_modes.py`](../../tests/test_corotational_nonlinear_job_tangent_modes.py).
- **GESDB:** Optional **`gesdb_kernel`** (`tl_locked` default, or `native` engineering/native weak form), and **`gesdb_tl_fallback`** (`true` / `false`): when `true`, [`GeometricallyExactShearDeformableBeam3D`](../../pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/geometrically_exact_shear_deformable_beam_3D.py) uses the parent TL strain hook only. Weak-form reference: [`gesdb_weak_form.md`](../element_library/gesdb_weak_form.md).
- **Buckling nonlinear prestress twins:** Optional **`buckling_nonlinear_prestress_twins`** (legacy: **`modal.buckling_nonlinear_prestress_twins`**) (`true` / `false`): with **`buckling_prestress = nonlinear_static`**, build nonlinear twin elements for the prestress solve while keeping linear elements for **`K`** and **`K_σ`** assembly (see [`MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md`](MODAL_BUCKLING_NONLINEAR_PRESTRESS_DESIGN.md)).
- **Runners** use **`run()`** as the main public step (not a generic `execute()` facade).
- **Settings:** [`simulation_settings_parser.py`](../../pre_processing/parsing/simulation_settings_parser.py); optional `[Newton]` for nonlinear static (`tolerance`, `max_iterations`, `tolerance_delta_u`, optional `relative_tolerance`, `relative_reference` = `first_residual` or `external_force`). Newton convergence uses the condensed RHS norm `‖F_cond‖` versus `tolerance + relative_tolerance × reference_scale`, not the full-vector `‖F_ext−F_int‖`. Optional `[Nonlinear]`: `num_increments` (uniform load factors `λ = 1/n … 1`), optional comma-separated `load_factors`, optional `line_search`, `line_search_max_backtracks`, `line_search_shrink` for discrete residual minimization along the Newton direction.
- **Logs (job outputs):** [`configure_child_logging`](../../workflow_orchestrator/run_job.py) attaches the root logger to `logs/process_job.log` and stdout so orchestrator and **`simulation_runner`** INFO lines (e.g. load increments and Newton iterations) appear together. Per-global-Newton metrics are also appended to **`logs/newton_history.csv`**. Each condensed linear solve appends to **`logs/SolveCondensedSystem.log`** (banner line per solve). TL Euler–Bernoulli / Timoshenko nonlinear elements have **no inner Newton loop** at the element—only Gauss quadrature sums in `internal_force_vector` / `tangent_stiffness_matrix`. **`primary_summary.csv`** includes optional nonlinear columns when present: `newton_iterations_total`, `newton_converged`, `load_increments_completed`.

**Research-oriented artifacts (join keys and provenance):**

| File | Role |
|------|------|
| `logs/newton_history.csv` | One row per **global Newton iteration** after that iteration’s update (columns include `load_increment_index`, `load_factor`, `newton_iter`, `norm_R_full`, `norm_F_cond`, `threshold`, `norm_delta_u`, `alpha`, `residual_ok`). **`residual_ok`** is `1` when the **condensed** residual norm `‖F_cond‖` satisfies the Newton tolerance test (`atol + rtol × reference`) **after** that iteration—it does **not** by itself mean “Newton finished”; convergence also requires `‖ΔU‖` below `tolerance_delta_u` where enforced. |
| `logs/inner_solve_history.csv` | One row per **condensed linear solve** (`K_cond δu = F_cond`): same **`load_increment_index`**, **`load_factor`**, **`newton_iter`** as `newton_history.csv` for nonlinear jobs (linear static leaves these blank). Columns include **`solver`**, **`preconditioner`**, **`iterations`** (`-1` if the chosen solver was **direct** from the start), **`final_relative_residual`** (iterative: last callback residual over first; direct: `‖r‖/‖F‖`), **`condensed_residual_norm`**, **`fallback_superlu`** (`1` if an iterative attempt fell back to SuperLU). Merge with `newton_history.csv` on `(load_increment_index, newton_iter, load_factor)`. |
| `logs/run_manifest.json` | Written at successful job end: **`wall_time_sec`**, **`git_commit`** (best-effort), Python / **NumPy** / **SciPy** versions, **`simulation_settings_txt_sha256`**, resolved **`simulation_settings`** dict, and paths to **`newton_history.csv`**, **`inner_solve_history.csv`**, **`primary_summary.csv`** when present. |
| `scripts/collect_run_metrics.py` | Optional: scans `post_processing/results/*/`, reads manifests and summaries, writes **`post_processing/runs_summary.csv`** for batch comparisons. |

Condensed-solve log banners include the same NR context when provided: `inc=…`, `nr=…`, `lambda=…`. Finer **per-increment energy / line-search dumps** remain a possible future **`[Nonlinear] diagnostic_level`** hook (not implemented by default).

---

## 6. Python style (short)

- Prefer **keyword-only** arguments after `*` in public constructors where practical.
- Use **`TYPE_CHECKING`** and quoted forward references for circular imports.
- **Validate early** (`ValueError` for bad input, `RuntimeError` for illegal call order).

---

## 7. Mixed meshes and factory-wide behaviour

- Assembly and runners must not assume a **single** element type unless documented.
- **Modal/dynamic:** every element in the mesh must implement `element_mass_matrix()` or the job must fail with a **clear** message listing unsupported types.

---

*See also [FORMULATION_DOCSTRING_STANDARDS.md](FORMULATION_DOCSTRING_STANDARDS.md) and [TESTING_STANDARDS.md](TESTING_STANDARDS.md).*
