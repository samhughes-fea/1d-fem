# Implementation plan: B2 shape function coefficients (Euler-Bernoulli only)

**Goal**: Store polynomial coefficients for N(ξ), dN/dξ, and d²N/dξ² on `ElementObject` for Euler-Bernoulli elements, plus a generic evaluator, so that post-processing can evaluate shape functions at arbitrary ξ after save/load without the element class. B1 (callable) remains supported; B2 adds serializable data.

**Scope**: Euler-Bernoulli 3D beam only. Timoshenko and Levinson are out of scope; they can be added later using the same format if applicable.

**Estimated effort**: 2–3 days.

### Visualization goal (why we need N(ξ) at arbitrary points)

In post-processing we want to plot:

1. **Discrete raw data as markers** – Gauss point values at their physical locations.
2. **Nodal values as markers** – with distinct styles for *projected* vs *explicit* nodal values.
3. **Continuous line in between** – The space between these points is drawn as a **continuous thin black line (linewidth 0.25) defined by shape-function interpolation**: u(ξ) = N(ξ) @ nodal_values (or the appropriate field), evaluated on a fine grid of ξ. This is the true FEM-interpolated field.

We do **not** rely on a generic curve fitter (numpy/matplotlib polyfit, spline, etc.) for the continuous curve. The continuous line is exactly what the element’s shape functions define. Storing N(ξ) (B1 callable or B2 coefficients) is what makes this possible.

**Marker convention** (to be used consistently in resolution/plotting code):

| Data type | Marker | Description |
|-----------|--------|-------------|
| Gauss point values | Small **solid circle** | Raw discrete values at Gauss point locations. |
| Projected nodal values | **Hollow circle** | Nodal values obtained by projection from Gauss (e.g. NodalResultProjector). |
| Explicit nodal values | **Solid triangle** | Nodal values that are not projections (e.g. primary solution at nodes). |
| Continuous interpolant | Thin black line (linewidth 0.25) | Shape-function interpolation between nodes; not a curve fit. |

This convention is implemented in `post_processing/graphical_visualisers/resolution_plotting_utils.py` (defaults and `nodal_data_type` parameter) and used by the deformation, strain, stress, and energy density profile visualisers.

---

## 1. Current state

- **ElementObject** ([`gauss_point_data.py`](pre_processing/element_library/gauss_point_data.py)) has optional `evaluate_shape_functions` (B1 callable). Euler-Bernoulli 3D sets it in [`euler_bernoulli_3D.py`](pre_processing/element_library/euler_bernoulli/euler_bernoulli_3D.py). Callables are not serializable.
- **EB shape functions** ([`euler_bernoulli/utilities/shape_functions.py`](pre_processing/element_library/euler_bernoulli/utilities/shape_functions.py)): N(ξ) shape (n_points, 12, 6). All entries are polynomials in ξ: linear (axial, torsion) or Hermite cubic (bending). Max degree 3. Explicit formulas in code (no numerical fitting needed).
- **Post-processing** ([`resolution_plotting_utils.py`](post_processing/graphical_visualisers/resolution_plotting_utils.py)): uses B1 callable when present via `_CacheShapeFunctionAdapter`; otherwise builds `ShapeFunctionOperator` from element_type + length.
- **Save formulation** ([`save_formulation_container.py`](processing_OOP/static/results/save_formulation_container.py)): saves K_e, F_e, B, D, J, xi, weight per GP. No load path exists yet in the repo.

---

## 2. Coefficient format (design)

**Convention**: Monomial basis in ξ. For each (dof_index, component_index) with 0 ≤ dof_index < 12, 0 ≤ component_index < 6, store coefficients for **ξ⁰, ξ¹, ξ², ξ³** (degree 3 max for EB).

- **N_coefficients**: `np.ndarray` shape `(12, 6, 4)`. `N_coefficients[dof, comp, k]` = coefficient of ξ^k in N_dof,comp(ξ). So `N(ξ) = sum_k N_coefficients[dof, comp, k] * ξ**k`.
- **dN_dxi_coefficients**: same shape `(12, 6, 4)` for dN/dξ (degree 2 max for EB cubics; higher indices zero).
- **d2N_dxi2_coefficients**: same shape `(12, 6, 4)` for d²N/dξ² (degree 1 max; higher indices zero).

**Storage on ElementObject**: Add optional fields (all or none set for a given element):

- `shape_function_N_coefficients: Optional[np.ndarray] = None`  # (12, 6, 4)
- `shape_function_dN_dxi_coefficients: Optional[np.ndarray] = None`  # (12, 6, 4)
- `shape_function_d2N_dxi2_coefficients: Optional[np.ndarray] = None`  # (12, 6, 4)

Alternatively a single dataclass or single array with a layout convention; the above keeps the plan simple and matches EB’s (N, dN_dξ, d2N_dξ2) API.

**EB mapping**: From [`shape_functions.py`](pre_processing/element_library/euler_bernoulli/utilities/shape_functions.py):

- (0,0), (6,0): linear 0.5(1±ξ) → [0.5, ∓0.5, 0, 0].
- (1,1), (7,1): N1 = 0.5 - 0.75ξ + 0.25ξ³ → [0.5, -0.75, 0, 0.25]; N3 = 0.5 + 0.75ξ - 0.25ξ³ → [0.5, 0.75, 0, -0.25].
- (5,5), (11,5): N2 = (L/8)(1 - ξ - ξ² + ξ³), N4 = -(L/8)(1 + ξ - ξ² - ξ³); depend on L, so coefficients are element-length-dependent and must be filled when building ElementObject (L is known there).
- Copy/negation rules: (2,2)=(1,1), (8,2)=(7,1); (4,4)=-(5,5), (10,4)=-(11,5); (3,3)=(0,0), (9,3)=(6,0). Export can fill the array by computing coefficients for the independent entries then copying.

---

## 3. Implementation tasks

### Phase 1: Data model and generic evaluator

| # | Task | File(s) | Notes |
|---|------|---------|--------|
| 1.1 | Add optional `shape_function_N_coefficients`, `shape_function_dN_dxi_coefficients`, `shape_function_d2N_dxi2_coefficients` to `ElementObject`. Document shape `(12, 6, 4)` and monomial convention. | [`gauss_point_data.py`](pre_processing/element_library/gauss_point_data.py) | Keep B1 `evaluate_shape_functions`; B2 is additive. |
| 1.2 | Implement generic evaluator: `evaluate_shape_functions_from_coefficients(N_coeffs, dN_coeffs, d2N_coeffs, xi) -> (N, dN_dξ, d2N_dξ2)`. Input `xi` shape (n_points,); output shapes (n_points, 12, 6). Use `np.polyval` or explicit loop over dof/comp. | New module e.g. [`pre_processing/element_library/utilities/shape_function_coefficient_evaluator.py`](pre_processing/element_library/utilities/shape_function_coefficient_evaluator.py) | Single function or small class; no element-type branching. |
| 1.3 | Unit test: for a known coefficient array (e.g. EB N1(ξ) = 0.5 - 0.75ξ + 0.25ξ³), call evaluator at several ξ and assert match to hand values or to `ShapeFunctionOperator.natural_coordinate_form`. | `tests/test_shape_function_coefficient_evaluator.py` or extend formulation cache tests | Validates evaluator correctness. |

### Phase 2: Euler-Bernoulli export

| # | Task | File(s) | Notes |
|---|------|---------|--------|
| 2.1 | In EB `element_stiffness_matrix()`, after building `gauss_cache`, build three arrays (12, 6, 4) from the known polynomials. Use `element_length` for N2, N4 terms. Fill copy/negation entries from the independent entries. | [`euler_bernoulli_3D.py`](pre_processing/element_library/euler_bernoulli/euler_bernoulli_3D.py) | Extract coefficients from formulas in `shape_functions.py`; no numerical fitting. |
| 2.2 | Pass the three coefficient arrays into `ElementObject(..., shape_function_N_coefficients=..., shape_function_dN_dxi_coefficients=..., shape_function_d2N_dxi2_coefficients=...)`. | Same | Keep existing B1 callable; add B2 fields for EB. |
| 2.3 | Unit test: build EB `ElementObject`, assert coefficient arrays are not None; at 3–5 sample ξ, compare `evaluate_shape_functions_from_coefficients(elem_obj.shape_function_*_coefficients, xi)` to `elem_obj.evaluate_shape_functions(xi)` (or to operator). | `tests/test_formulation_cache_shape_functions.py` or new test file | Ensures EB export matches operator. |

### Phase 3: Post-processing and serialization

| # | Task | File(s) | Notes |
|---|------|---------|--------|
| 3.1 | In `get_shape_function_operator_for_element`: when `element_object` has B2 coefficients (all three arrays non-None), return an adapter that uses `evaluate_shape_functions_from_coefficients` with those arrays. Prefer B2 over B1 when both present (so that after load, coefficients path is used). Order: B2 coefficients → B1 callable → create operator from type+length. | [`resolution_plotting_utils.py`](post_processing/graphical_visualisers/resolution_plotting_utils.py) | Extend `_CacheShapeFunctionAdapter` or add a second adapter that takes coefficients; or one adapter that accepts either callable or (N_coeffs, dN_coeffs, d2N_coeffs). |
| 3.2 | In `SaveFormulationContainer` (or equivalent), when saving per-element data, save the three coefficient arrays for each element that has them (e.g. one CSV or NPY per array per element, or one file with element index). Document layout in header or companion doc. | [`save_formulation_container.py`](processing_OOP/static/results/save_formulation_container.py) | Only save if arrays are present. |
| 3.3 | **(Optional / when load exists)** When a load path for formulation cache is implemented, read the coefficient arrays and attach to reconstructed `ElementObject`s. | Future load module | Can be done in a follow-up when load is added. |

### Phase 4: Documentation and checklist

| # | Task | File(s) | Notes |
|---|------|---------|--------|
| 4.1 | Document B2 coefficient format (shape, convention, EB-only for now) in `gauss_point_data.py` docstring and/or in `RESULTS_DESIGN.md` or this plan. | [`gauss_point_data.py`](pre_processing/element_library/gauss_point_data.py), [`RESULTS_DESIGN.md`](processing_OOP/static/results/RESULTS_DESIGN.md) | So future elements (Timoshenko, Levinson) can adopt same format if desired. |
| 4.2 | Add “B2 coefficients (EB)” to the formulation cache / “Adding a new element” checklist if present: optional; EB implements it for serializable arbitrary-point evaluation. | [`formulation_cache_shape_functions_implementation_plan.md`](docs/plans/formulation_cache_shape_functions_implementation_plan.md) or README | Optional checklist item. |

---

## 4. File touch list

| File | Change |
|------|--------|
| `pre_processing/element_library/gauss_point_data.py` | Add three optional coefficient fields to `ElementObject`; docstrings. |
| `pre_processing/element_library/utilities/shape_function_coefficient_evaluator.py` | **New.** Generic evaluator from (N_coeffs, dN_coeffs, d2N_coeffs, xi) → (N, dN_dξ, d2N_dξ2). |
| `pre_processing/element_library/euler_bernoulli/euler_bernoulli_3D.py` | Build and pass coefficient arrays when creating `ElementObject`. |
| `post_processing/graphical_visualisers/resolution_plotting_utils.py` | Prefer B2 coefficients when present; use generic evaluator in adapter path. |
| `processing_OOP/static/results/save_formulation_container.py` | Save the three coefficient arrays per element when present. |
| `tests/test_formulation_cache_shape_functions.py` or new test file | Tests for evaluator and for EB ElementObject B2 export vs operator. |
| `processing_OOP/static/results/RESULTS_DESIGN.md` | Short note on B2 coefficient format (optional). |

---

## 5. Success criteria

- EB `ElementObject` can carry optional B2 coefficient arrays; evaluator returns (N, dN_dξ, d2N_dξ2) matching `ShapeFunctionOperator` at arbitrary ξ.
- Post-processing uses B2 when coefficients are present (e.g. after load), else B1 callable, else operator from type+length.
- Coefficient arrays for EB are saved with the formulation; load path can be added later.
- All existing formulation cache and EB tests pass; new tests cover evaluator and EB B2 export.

---

## 6. Dependencies and order

1. Phase 1.1 and 1.2 can be done in parallel (data model + evaluator).
2. Phase 1.3 (evaluator test) depends on 1.2.
3. Phase 2 depends on 1.1 and 1.2 (need evaluator to validate export).
4. Phase 3.1 depends on 1.2 and 2; Phase 3.2 depends on 1.1.
5. Phase 4 can be done alongside or after Phase 3.
