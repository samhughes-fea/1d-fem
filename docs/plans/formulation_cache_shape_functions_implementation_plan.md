# Implementation plan: Shape functions in formulation cache

**Goal**: Ensure every element type that contributes to the formulation cache stores the shape functions (and their evaluated magnitudes at Gauss points) used to build K_e and F_e, so the cache is the single source of truth and downstream (e.g. nodal projection) can use element-consistent data.

---

## 1. Current state

| Component | Role |
|-----------|------|
| **StiffnessGaussPointData** | Per-GP data for stiffness: `xi`, `weight`, `B_matrix`, `D_matrix`, `jacobian`, optional `shape_functions`, `shape_derivatives`. |
| **ForceGaussPointData** | Per-GP data for force: `xi`, `weight`, `shape_functions`, `jacobian`, optional `distributed_load`. |
| **ElementObject** | Per element: `element_id`, `element_type`, `K_e`, `gauss_data: List[StiffnessGaussPointData]`. |
| **ForceObject** | Per element: `element_id`, `element_type`, `F_e`, `gauss_data: List[ForceGaussPointData]`. |

**Element types that return ElementObject / ForceObject** (via `ElementFactory.ELEMENT_CLASS_MAP` and used in static pipeline):

- **EulerBernoulliBeamElement3D** – already fills `shape_functions` and `shape_derivatives` in `StiffnessGaussPointData`, and `shape_functions` in `ForceGaussPointData`.
- **TimoshenkoBeamElement3D** – same.
- **LevinsonBeamElement3D** – same.

**Element types that currently return `np.ndarray`** (not ElementObject) and are not in the static formulation-cache path:

- co_tide_beam_C, co_tide_beam_ML, tidal_benchmark_blade_C, tidal_benchmark_blade_ML (return type hint `np.ndarray` from `element_stiffness_matrix()`; not in `ELEMENT_CLASS_MAP` for the static runner used here).

---

## 2. Principles

1. **Single source of truth**: The formulation cache (ElementObject / ForceObject and their `gauss_data`) is the only place that defines what shape functions and constituent tensors (B, D, J) were used to build K_e and F_e.
2. **Per-element, per-Gauss**: Shape function values N(ξ) and derivatives dN/dξ (where used) are stored at each Gauss point for that element’s formulation.
3. **Mandatory for new elements**: Any new element type that returns ElementObject / ForceObject must populate `shape_functions` (and `shape_derivatives` where applicable) in every Gauss point record.
4. **Result computers are formulation-agnostic**: Result computers (secondary and tertiary) can handle **any number n** of Gauss points per element and have access to **whatever shape functions** (and derivatives) were used to assemble that element. They iterate over `gauss_data` and use the stored N, dN/dξ, B, D, J, etc.; they do not assume a fixed quadrature rule or element type.

---

## 3. Implementation tasks

### Phase 1: Contract, validation, and documentation

| # | Task | Owner / notes |
|---|------|----------------|
| 1.1 | **Formalise the contract** in `pre_processing/element_library/gauss_point_data.py`: in docstrings for `StiffnessGaussPointData` and `ElementObject`, state that when used in the formulation cache for the results pipeline, `shape_functions` (and `shape_derivatives` for stiffness) must be populated at every Gauss point. | Docstring update |
| 1.2 | **Add validation** when building or using the formulation cache: e.g. in `FormulationResultSet` or in the static runner when constructing it, loop over all `element_objects` and their `gauss_data`; if any `StiffnessGaussPointData.shape_functions` is `None`, log a warning (and optionally raise in strict mode). Same idea for `ForceObject.gauss_data` and `shape_functions`. | New helper in `processing_OOP.static.results.containers.formulation_results` or in the runner; call once when `FormulationResultSet` is built. |
| 1.3 | **Document the requirement** in `processing_OOP/static/results/RESULTS_DESIGN.md`: add a short subsection “Formulation cache: shape functions” stating that ElementObject/ForceObject must store N (and dN/dξ where relevant) at each Gauss point so that the cache is the single source of truth and projection can be element-consistent. | RESULTS_DESIGN.md |
| 1.4 | **Add a unit test** (or extend an existing formulation test) that, for at least one element type (e.g. Euler-Bernoulli 3D), builds an ElementObject and a ForceObject and asserts that every entry in `gauss_data` has non-None `shape_functions` (and for stiffness, non-None `shape_derivatives` where the element type uses them). | e.g. under `tests/` or inside `pre_processing/element_library/` |

### Phase 2: Future element types and optional projector change

| # | Task | Owner / notes |
|---|------|----------------|
| 2.1 | **Onboard any new element type** that returns ElementObject/ForceObject: when implementing `element_stiffness_matrix()` / `element_force_vector()`, ensure each `StiffnessGaussPointData` and `ForceGaussPointData` is built with `shape_functions` (and stiffness with `shape_derivatives` if used). Treat this as a checklist item in code review or a short “Adding a new element” note in the element library. | Per-element PRs; optional: add a short “Adding a new element” section in `pre_processing/element_library/README` or in this plan. |
| 2.2 | **(Optional) Element-consistent nodal projection**: In `NodalResultProjector`, when `gp.shape_functions` is present for all Gauss points of an element, use an element-consistent extrapolation (e.g. values_at_nodes = N_gauss^{-1} @ values_at_gauss using the stored N at Gauss points and N at node natural coords) instead of generic Lagrange in ξ. Fall back to current Lagrange interpolation when `shape_functions` is missing. | `processing_OOP/static/results/compute_secondary/nodal_result_projector.py` |

#### Adding a new element (checklist)

When implementing a new element type that returns **ElementObject** (from `element_stiffness_matrix()`) or **ForceObject** (from `element_force_vector()`):

1. Build each **StiffnessGaussPointData** with `shape_functions` and `shape_derivatives` set at **every** Gauss point.
2. Build each **ForceGaussPointData** with `shape_functions` set at **every** Gauss point.
3. Treat this as **mandatory** for the results pipeline; validation is performed by `validate_shape_functions_populated()` when the formulation cache is built.

See the contract in `pre_processing/element_library/gauss_point_data.py` (StiffnessGaussPointData, ForceGaussPointData, ElementObject docstrings).

### Phase 3: Optional extensions (out of scope for minimal plan)

- **N(ξ) at node locations**: If needed for extrapolation, store or compute N at node natural coordinates (e.g. ξ = -1, +1) per element type; or compute on the fly in the projector from the element’s shape function API when available.
- **Co-tide / tidal_benchmark**: If these are later integrated into the static formulation-cache pipeline, they must be changed to return ElementObject/ForceObject with full `gauss_data` including `shape_functions` (and `shape_derivatives` for stiffness), following the same contract as EB/Timoshenko/Levinson.

---

## 4. File touch list

| File | Change |
|------|--------|
| `pre_processing/element_library/gauss_point_data.py` | Docstring updates (StiffnessGaussPointData, ElementObject, ForceGaussPointData, ForceObject) stating that shape_functions/shape_derivatives must be set when used in the results formulation cache. |
| `processing_OOP/static/results/containers/formulation_results.py` | Optional: add `validate_shape_functions_populated(element_objects, force_objects)` and call from runner or from FormulationResultSet. |
| `processing_OOP/static/results/RESULTS_DESIGN.md` | New subsection “Formulation cache: shape functions” (contract + where they are stored). |
| `processing_OOP/static/results/compute_secondary/nodal_result_projector.py` | Phase 2 only: use `gp.shape_functions` when present for extrapolation; else keep current Lagrange in ξ. |
| New or existing test file | Phase 1: test that EB (and optionally Timoshenko/Levinson) ElementObject and ForceObject have shape_functions (and shape_derivatives for stiffness) non-None in all gauss_data. |

---

## 5. Success criteria

- **Phase 1**: Contract is documented; validation runs when the formulation cache is built (warning or strict); RESULTS_DESIGN.md explains the requirement; at least one test enforces non-None shape_functions (and shape_derivatives where applicable) for one element type.
- **Phase 2**: Any new element returning ElementObject/ForceObject populates shape data; optionally, NodalResultProjector uses stored shape functions when available.
- **Ongoing**: No new element type is added to the formulation-cache pipeline without populating shape_functions (and shape_derivatives for stiffness) at every Gauss point.

---

## 6. Dependencies and order

1. Phase 1.1 and 1.3 can be done in parallel (docstrings + RESULTS_DESIGN.md).
2. Phase 1.2 (validation) can be implemented next; it depends only on the existing structure of ElementObject/ForceObject.
3. Phase 1.4 (test) can be done in parallel or after 1.1.
4. Phase 2.2 (projector) depends on 1.2 so that we can assume shape_functions are present when strict validation is on; otherwise it can fall back to Lagrange when missing.
