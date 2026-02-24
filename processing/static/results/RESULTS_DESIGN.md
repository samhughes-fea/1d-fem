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

### B2 shape function coefficients (optional)

ElementObject may optionally carry **B2** coefficient arrays for N, dN/dxi, and d2N/dxi2 in monomial basis (shape (12, 6, 4) for 1D beams). These allow evaluating shape functions at arbitrary xi after save/load without the element class. Implemented for Euler-Bernoulli 3D; see `docs/plans/b2_shape_function_coefficients_eb_implementation_plan.md`. Generic evaluator: `evaluate_shape_functions_from_coefficients` in `pre_processing/element_library/utilities/shape_function_coefficient_evaluator.py`.
