# Results pipeline design

Short reference for where result types live and how optional/unused computers are handled.

## Resolution levels

- **Global**: primary only (K, F, U, R_global, R_residual, etc.).
- **Gaussian**: secondary (strain, stress, energy density); tertiary (section forces, principal stress).
- **Nodal**: secondary (projected strain, stress, energy density).
- **Elemental**: primary (K_e, F_e, U_e, R_e, R_residual_e); tertiary (total strain energy, integrated section forces).

## Unused / optional result computers

These modules exist but are **not** called by any orchestrator. The design is explicit so that future wiring or removal is clear.

| Module | Location | Status | Notes |
|--------|----------|--------|--------|
| **ComputeTotalStrainEnergyPerElement** | `compute_secondary/total_strain_energy.py` | **Deprecated** | Total strain energy per element is computed in the **tertiary** pipeline by `ComputeIntegratedElementalResults` (same formula, current API). This module expected a legacy `element_dictionary` and dict-like Gaussian results; do not use. |
| **ComputeNodalForce** (and **ComputeGlobalNodalForce**, **ComputeNodalReactionBalance**) | `compute_secondary/nodal_force.py` | **Optional / unwired** | Computes F_int = K_e @ U_e (nodal internal forces). Uses only primary data (K_e, U_e). If needed for equilibrium checks or reporting, wire into the **primary** pipeline (e.g. after disassembly) or a dedicated step; not part of secondary. |

## Summary savers

Each level has a summary CSV written in the run:

- **Primary**: `SavePrimaryResultsSummary` → `primary_results/primary_summary.csv`
- **Secondary**: `SaveSecondaryResultsSummary` → `secondary_results/secondary_summary.csv`
- **Tertiary**: `SaveTertiaryResultsSummary` → `tertiary_results/tertiary_summary.csv`

See `save_*_container.py` in this directory.

## Formulation cache: shape functions

The formulation cache (ElementObject / ForceObject and their per-Gauss `gauss_data`) is the single source of truth for how K_e and F_e were assembled. When ElementObject / ForceObject are used in the results pipeline, they **must** store N (and dN/dξ for stiffness) at each Gauss point so the cache is the single source of truth and projection can be element-consistent. Result computers (secondary and tertiary) are **formulation-agnostic**:

- They support **any number n** of Gauss points per element; they iterate over `gauss_data` and do not assume a fixed quadrature rule.
- They use **whatever shape functions** (and derivatives) were used to assemble the element: N(ξ), dN/dξ, B, D, J, etc. are read from the cache at each Gauss point.

So result computers always see the same formulation that built the element; no duplicate shape-function logic is required in the results pipeline. For a short **checklist when adding a new element** that returns ElementObject/ForceObject, see the Phase 2 subsection “Adding a new element (checklist)” in `docs/plans/formulation_cache_shape_functions_implementation_plan.md`.

### Warping (7 DOF per node) and cached `shape_functions`

Beam elements with Vlasov warping use an enlarged **`B`** (e.g. 7×14) while **`shape_functions`** on each `StiffnessGaussPointData` remain the **baseline beam operator values** `N(ξ)`, `dN/dξ` from the element’s shape-function operator (same 12 rows × 6 columns style as the non-warping Timoshenko/Euler case). Tertiary **`NodalSectionForcesProjector`** uses those blocks for GP→nodal projection; χ DOFs do not add extra rows to the cached `N` slice.

### Strict validation (optional)

After building `FormulationResultSet`, static runners call `validate_shape_functions_populated`. By default missing `shape_functions` logs a **warning**. Set environment variable **`FEM_FORMULATION_CACHE_STRICT_SHAPE=1`** (or `true` / `yes` / `on`) to **raise** `ValueError` on the first missing entry (CI / debugging).

### Modal, dynamic, and harmonic pipelines (optional formulation-cache post-processing)

By default, modal, dynamic, and harmonic jobs **do not** run the static secondary/tertiary CSV pipeline. When **`simulation_settings["post_processing"]["run_secondary_tertiary_modal"]`** is **true** — parser aliases **`run_secondary_tertiary_eigen`** and **`run_secondary_tertiary_buckling`** set the same internal key — when **`run_secondary_tertiary_dynamic`** is **true**, or when **`run_secondary_tertiary_harmonic`** is **true**, the runner builds a **`FormulationResultSet`** from the job’s **`element_objects`** / **`force_objects`** (same wiring as static), calls **`validate_shape_functions_populated`** with **`FEM_FORMULATION_CACHE_STRICT_SHAPE`** semantics, and drives **`SecondaryResultsOrchestrator`** + tertiary savers using a **snapshot displacement** **`U_global`**:

| Analysis | Snapshot **`U_global`** | Settings |
|----------|-------------------------|----------|
| Modal vibration | Column **`mode_shapes[:, modal_mode_index] × modal_amplitude`** | Modal BCs use the penalty method on full-length vectors; **`eigsh`** modes are already **`n_dof`** long. |
| Modal buckling | **`buckling_displacement = mode`**: column **`buckling_modes[:, buckling_mode_index] × modal_amplitude`**; **`prestress`**: prestress **`U`** after the linear/nonlinear static prestress solve, scaled by **`modal_amplitude`** | **`buckling_mode_index`** defaults to **`modal_mode_index`** when omitted. |
| Dynamic (Newmark) | Row **`U[dynamic_time_index, :]`** or multiple rows from **`dynamic_time_indices`** | **`dynamic_time_index`** default **`-1`** (last time step); negative indices count from the end of the displacement history. When **`dynamic_time_indices`** is set (comma-separated list), each listed row is passed to the formulation-cache pipeline; multiple rows write under **`secondary_results/dynamic_post/t_<row>/`** (and matching **`tertiary_results/`** subtree) so outputs do not overwrite each other. |
| Harmonic (§4) | Snapshot displacement for secondary/tertiary comes from **`real(U[:, k])`**, **`imag(U[:, k])`**, or **both** (separate exports) per **`harmonic_secondary_tertiary_displacement_component`** (`real` \| `imag` \| `both`). Columns **`k`** are chosen like modal: default **`harmonic_frequency_index`**, **`harmonic_secondary_tertiary_frequency_indices`**, or **`harmonic_secondary_tertiary_all_frequencies`** (writes under **`secondary_results/harmonic_post/freq_*`**; imaginary-part snapshots use **`freq_*_imag`**). | Priority when multiple keys are set: explicit **indices** > **all_frequencies** > single **harmonic_frequency_index** (default **`0`**). Primary **`harmonic_results/`** real/imag/abs/phase matrices remain the source of truth for magnitude/phase workflows. |

Outputs mirror static layout: **`secondary_results/`**, **`tertiary_results/`** under the job results root.

#### Eigen secondary (undamped vibration) — `primary_results/modal_results/`

When **`run_secondary_tertiary_modal`** (parser aliases **`run_secondary_tertiary_eigen`** / **`run_secondary_tertiary_buckling`**) is **false**, §2 vibration still writes lightweight **secondary** scalars as text files alongside frequencies and mode shapes:

| File | Content | Notes |
|------|---------|-------|
| `{job_name}_modal_generalized_mass.txt` | One value per computed mode: \( \phi_j^{\mathsf T} M_{\mathrm{mod}} \phi_j \) using the **penalty-modified** mass from the eigen BC step. | Uses the same \(M_{\mathrm{mod}}\) as the generalized eigenproblem after BCs; not a rigid-body physical modal mass. |
| `{job_name}_modal_load_participation.txt` | One value per mode: \( \lvert \hat{\phi}_j^{\mathsf T} F \rvert \) with \(\hat{\phi}_j = \phi_j / \sqrt{\phi_j^{\mathsf T} M_{\mathrm{mod}} \phi_j}\); \(F\) is the **assembled element nodal load** from `force_objects` / `F_e` scatter ([`assemble_global_force_vector`](../../dynamic/assembly.py)), identical in spirit to the transient reference load. | All zeros if `force_objects` is missing or \(\lVert F \rVert\) is negligible. **Not** building-code seismic modal participation (no prescribed rigid-body vector or code normalization); treat as a **load–mode alignment** indicator only. |
| `{job_name}_modal_effective_mass_fraction_z.txt` | One value per mode: \((\phi_j^{\mathsf T} M r)^2 / ((r^{\mathsf T} M r)(\phi_j^{\mathsf T} M \phi_j))\) with \(r\) a **unit pattern on translational +Z** DOFs (six DOFs per node: activates index `uz` per node). | **Native metric** (implemented in `processing.eigen.metrics`); omitted for unsupported DOF layouts (e.g. warping / seven DOF per node until a pattern is defined). |

See also [SIMULATION_SETTINGS_TAXONOMY.md](../../../docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md) and [simulation_runner/transient/README.md](../../../simulation_runner/transient/README.md).

Parser: **`[PostProcessing]`** (also **`[post_processing]`**, spacing-insensitive) with keys **`run_secondary_tertiary_modal`** (aliases **`run_secondary_tertiary_eigen`**, **`run_secondary_tertiary_buckling`**), **`run_secondary_tertiary_dynamic`**, **`run_secondary_tertiary_harmonic`**, **`modal_mode_index`**, **`modal_amplitude`**, **`buckling_displacement`** (`mode` \| `prestress`), **`buckling_mode_index`**, **`dynamic_time_index`**, **`dynamic_time_indices`** (comma-separated row indices into **`U`**; when present, selects one or more snapshots), **`harmonic_frequency_index`**, **`harmonic_secondary_tertiary_all_frequencies`**, **`harmonic_secondary_tertiary_frequency_indices`** (comma-separated column indices), **`harmonic_secondary_tertiary_displacement_component`** (`real` \| `imag` \| `both`).

#### Harmonic tertiary — complex-valued stress recovery (design backlog)

Today’s pipeline feeds **real** **`U_global`** snapshots into the same linear strain recovery used for static jobs (`real(U)`, optionally **`imag(U)`** as a second pass). **Future extension:** true harmonic stress harmonics could follow either **(i)** complex-valued element strains/stresses derived from \(\mathbf{u}=\mathbf{u}_r+i\mathbf{u}_i\) with linear hooks evaluated separately on \(\mathbf{u}_r\) and \(\mathbf{u}_i\), or **(ii)** amplitude/phase reconstruction from two phases (e.g. \(0\) and \(\pi/2\)) without changing element kernels. Choice affects **`SecondaryResultsOrchestrator`** / tertiary contracts — gate behind **`[PostProcessing]`** keys after **`RESULTS_DESIGN`** + parity tests.

#### Wiring checklist

1. Every **`ElementObject`** / **`ForceObject`** must populate **`shape_functions`** (and stiffness **`shape_derivatives`**) at Gauss points used for recovery — see [`formulation_cache_shape_functions_implementation_plan.md`](../../../docs/plans/formulation_cache_shape_functions_implementation_plan.md).
2. **`NodalResultProjector`** prefers cache **`shape_functions`** when complete; otherwise it uses the Lagrange fallback in `compute_secondary/nodal_result_projector.py`.
3. Run **`tests/test_formulation_cache_shape_functions.py`** and modal/dynamic post-processing tests with **`FEM_FORMULATION_CACHE_STRICT_SHAPE=1`** where CI parity matters.

### B2 shape function coefficients (optional)

ElementObject may optionally carry **B2** coefficient arrays for N, dN/dxi, and d2N/dxi2 in monomial basis (shape (12, 6, 4) for 1D beams). These allow evaluating shape functions at arbitrary xi after save/load without the element class. Implemented for Euler-Bernoulli 3D; see `docs/plans/b2_shape_function_coefficients_eb_implementation_plan.md`. Generic evaluator: `evaluate_shape_functions_from_coefficients` in `pre_processing/element_library/utilities/shape_function_coefficient_evaluator.py`.
