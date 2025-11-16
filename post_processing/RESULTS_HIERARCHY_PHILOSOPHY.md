# Results Hierarchy Philosophy Documentation

## Executive Summary

This document analyzes the results directory hierarchy philosophy, which is based on **native resolution** - the smallest resolution where each quantity is first computed. The hierarchy supports both native results and their projections/integrations to other resolutions.

## Core Philosophy

### Native Resolution Principle

The directory hierarchy (`global/`, `elemental/`, `nodal/`, `gaussian/`) represents the **native resolution** where quantities are **first computed** (smallest resolution). This is distinct from where they can be projected or integrated.

### Resolution Transformations

Results can be transformed to other resolutions using cached shape functions:
- **Nodal projections**: Interpolate/extrapolate from Gaussian to nodes
- **Element-level sums/integrations**: Integrate quantities over element domain via quadrature

### Clear Distinction Required

The directory structure must distinguish between:
1. **Native resolution results**: Computed directly at their native resolution
2. **Projected results**: Interpolated/extrapolated from native using shape functions
3. **Integrated results**: Summed/integrated over element domain via quadrature

---

## Terminology Definitions

### Global Resolution
- **Definition**: System-wide quantities computed at the assembled system level
- **Native quantities**: 
  - Displacements (`U_global`) - solved from `K_global @ U_global = F_global`
  - Reactions (`R_global`) - computed from `R_global = K_global @ U_global - F_global`
  - Global matrices (`K_global`, `F_global`) - assembled from element contributions
- **Shape**: `(total_dofs,)` for vectors, `(total_dofs, total_dofs)` for matrices

### Elemental Resolution
- **Definition**: Per-element quantities, either:
  - Disassembled from global (e.g., `U_e = U_global[dof_indices]`)
  - Integrated from Gaussian (e.g., total strain energy per element)
- **Native quantities**:
  - Element displacements (`U_e`) - disassembled from global
  - Element reactions (`R_e`) - disassembled from global
  - Total strain energy per element - integrated from Gaussian energy density
- **Shape**: `List[element] -> (n_dofs_per_element,)` or `List[element] -> float`

### Gaussian Resolution
- **Definition**: Quantities computed at integration points (Gauss points) within elements
- **Native quantities**:
  - Strain (`ε`) - computed via `ε = B @ U_e` at each Gauss point
  - Stress (`σ`) - computed via `σ = D @ ε` at each Gauss point
  - Energy density (`w`) - computed via `w = 0.5 * ε^T @ σ` at each Gauss point
- **Shape**: `List[element] -> List[gauss_point] -> (n_components,)`

### Nodal Resolution
- **Definition**: Quantities interpolated/extrapolated to nodes
- **Native quantities**: None (all nodal results are projections)
- **Projected quantities**:
  - Nodal strain - projected from Gaussian using shape functions
  - Nodal stress - projected from Gaussian using shape functions
  - Nodal energy density - projected from Gaussian using shape functions
- **Shape**: `(n_nodes, n_components)` for tensors, `(n_nodes,)` for scalars

---

## Current Directory Structure Analysis

### Primary Results Structure

```
primary_results/
├── global/
│   ├── K_global.csv          # Native: Global stiffness matrix
│   ├── K_mod.csv             # Native: Modified stiffness matrix
│   ├── K_cond.csv            # Native: Condensed stiffness matrix
│   ├── F_global.csv          # Native: Global force vector
│   ├── F_mod.csv             # Native: Modified force vector
│   ├── F_cond.csv            # Native: Condensed force vector
│   ├── U_global.csv          # Native: Global displacement vector
│   ├── U_cond.csv            # Native: Condensed displacement vector
│   ├── R_global.csv          # Native: Global reaction vector
│   └── R_residual.csv        # Native: Global residual vector
├── elemental/
│   ├── element_stiffness/    # Native: Element stiffness matrices (formulation)
│   ├── external_force/       # Native: Element force vectors (formulation)
│   ├── deformation/          # Derived: Element displacements (disassembled from global)
│   ├── reaction_force/      # Derived: Element reactions (disassembled from global)
│   └── residual/             # Derived: Element residuals (disassembled from global)
└── formulation/
    ├── element_types.csv
    ├── stiffness/            # Cached element stiffness matrices
    └── force/                # Cached element force vectors
```

**Analysis**: 
- ✅ Global results correctly stored at native resolution
- ⚠️ Elemental results mix native (formulation) and derived (disassembled) quantities
- ✅ Formulation data correctly separated

### Secondary Results Structure

```
secondary_results/
└── secondary_results/        # ❌ ISSUE: Nested directory
    ├── gaussian/
    │   ├── strain/           # ✅ Native: Strain at Gauss points
    │   ├── stress/           # ✅ Native: Stress at Gauss points
    │   └── energy_density/   # ✅ Native: Energy density at Gauss points
    ├── nodal/
    │   ├── nodal_strain.csv      # ✅ Projected: Strain at nodes
    │   ├── nodal_stress.csv     # ✅ Projected: Stress at nodes
    │   └── nodal_strain_energy_density.csv  # ✅ Projected: Energy at nodes
    └── elemental/
        └── element_strain_energy.csv  # ✅ Integrated: Total energy per element
```

**Analysis**:
- ✅ Gaussian results correctly stored at native resolution
- ✅ Nodal results correctly identified as projections
- ✅ Elemental results correctly identified as integrations
- ❌ **CRITICAL ISSUE**: Nested `secondary_results/secondary_results/` directory

### Root Cause of Nested Directory

**Location**: `processing_OOP/static/results/save_secondary_container.py:46`

```python
self.secondary_dir = self.save_dir / "secondary_results"
```

**Problem**: `save_dir` is already `secondary_results_dir` (set in `static_simulation.py:83`), so this creates:
- `secondary_results_dir` = `{job_results_dir}/secondary_results`
- `self.secondary_dir` = `{job_results_dir}/secondary_results/secondary_results`

**Fix**: Should be `self.secondary_dir = self.save_dir` (since `save_dir` is already the secondary results directory)

---

## Computation-to-Storage Mapping

### Primary Results

| Quantity | Native Resolution | Computation | Storage Location | Status |
|----------|------------------|-------------|------------------|--------|
| `U_global` | Global | Solved from `K_global @ U_global = F_global` | `primary_results/global/U_global.csv` | ✅ Correct |
| `R_global` | Global | `R_global = K_global @ U_global - F_global` | `primary_results/global/R_global.csv` | ✅ Correct |
| `U_e` | Elemental (derived) | `U_e = U_global[dof_indices]` | `primary_results/elemental/deformation/` | ✅ Correct |
| `R_e` | Elemental (derived) | `R_e = R_global[dof_indices]` | `primary_results/elemental/reaction_force/` | ✅ Correct |
| `K_e` | Elemental (formulation) | Element formulation | `primary_results/elemental/element_stiffness/` | ✅ Correct |
| `F_e` | Elemental (formulation) | Element formulation | `primary_results/elemental/external_force/` | ✅ Correct |

### Secondary Results

| Quantity | Native Resolution | Computation | Storage Location | Status |
|----------|------------------|-------------|------------------|--------|
| `ε` (strain) | Gaussian | `ε = B @ U_e` at each Gauss point | `secondary_results/gaussian/strain/` | ✅ Correct |
| `σ` (stress) | Gaussian | `σ = D @ ε` at each Gauss point | `secondary_results/gaussian/stress/` | ✅ Correct |
| `w` (energy density) | Gaussian | `w = 0.5 * ε^T @ σ` at each Gauss point | `secondary_results/gaussian/energy_density/` | ✅ Correct |
| `ε_nodal` | Nodal (projected) | Interpolated from Gaussian | `secondary_results/nodal/nodal_strain.csv` | ✅ Correct |
| `σ_nodal` | Nodal (projected) | Interpolated from Gaussian | `secondary_results/nodal/nodal_stress.csv` | ✅ Correct |
| `w_nodal` | Nodal (projected) | Interpolated from Gaussian | `secondary_results/nodal/nodal_strain_energy_density.csv` | ✅ Correct |
| `U_total` (per element) | Elemental (integrated) | `∫ w dΩ` via quadrature | `secondary_results/elemental/element_strain_energy.csv` | ✅ Correct |

---

## Identified Issues

### 1. Nested Directory Structure (CRITICAL)

**Issue**: `secondary_results/secondary_results/` nested directory

**Location**: `processing_OOP/static/results/save_secondary_container.py:46`

**Current Code**:
```python
self.secondary_dir = self.save_dir / "secondary_results"
```

**Root Cause**: `save_dir` parameter is already `secondary_results_dir`, so adding `/secondary_results` creates nesting.

**Fix**: Change to:
```python
self.secondary_dir = self.save_dir  # save_dir is already secondary_results_dir
```

### 2. Terminology Clarity

**Issue**: The distinction between "elemental" as:
- Disassembled quantities (derived from global)
- Integrated quantities (summed from Gaussian)
- Formulation quantities (computed during element setup)

**Recommendation**: Consider subdirectories or naming conventions:
- `elemental/disassembled/` - for U_e, R_e (from global)
- `elemental/integrated/` - for total strain energy (from Gaussian)
- `elemental/formulation/` - for K_e, F_e (element setup)

### 3. Missing Documentation

**Issue**: No clear documentation of which quantities are native vs. projected vs. integrated.

**Recommendation**: Add README files in each resolution directory explaining:
- What quantities are stored
- Whether they are native, projected, or integrated
- How they were computed

---

## Recommended Directory Structure

### Proposed Structure (After Fixes)

```
primary_results/
├── global/                    # Native: System-wide quantities
│   ├── matrices/             # K_global, K_mod, K_cond
│   └── vectors/              # U_global, F_global, R_global, etc.
├── elemental/
│   ├── disassembled/        # Derived: U_e, R_e (from global)
│   └── formulation/         # Native: K_e, F_e (element setup)
└── formulation/             # Cached: B, D matrices, Gauss point data

secondary_results/
├── gaussian/                # Native: Quantities at Gauss points
│   ├── strain/
│   ├── stress/
│   └── energy_density/
├── nodal/                   # Projected: Quantities at nodes
│   ├── strain/
│   ├── stress/
│   └── energy_density/
└── elemental/               # Integrated: Quantities per element
    └── total_strain_energy/
```

### Alternative: Flat Structure with Metadata

```
secondary_results/
├── strain/
│   ├── gaussian/           # Native resolution
│   ├── nodal/              # Projected resolution
│   └── elemental/          # Integrated resolution (if applicable)
├── stress/
│   ├── gaussian/           # Native resolution
│   ├── nodal/              # Projected resolution
│   └── elemental/          # Integrated resolution (if applicable)
└── energy_density/
    ├── gaussian/           # Native resolution
    ├── nodal/              # Projected resolution
    └── elemental/          # Integrated resolution
```

**Recommendation**: Keep current structure (quantity-first) as it aligns with the native resolution philosophy. The nested directory issue is a bug, not a design flaw.

---

## Tertiary Results Structure

### Current Structure (Resolution-First)

```
tertiary_results/
├── gaussian/                    # ✅ Native: Gaussian resolution
│   ├── section_forces/
│   │   └── section_forces_elem_*.csv
│   └── principal_stress/
│       └── principal_stress_elem_*.csv
├── elemental/                   # ✅ Native: Elemental resolution (integrated)
│   ├── total_strain_energy.csv
│   └── integrated_section_forces.csv
└── tertiary_summary.csv         # Summary statistics
```

**Analysis**:
- ✅ Section forces correctly stored at native Gaussian resolution
- ✅ Principal stresses correctly stored at native Gaussian resolution
- ✅ Integrated elemental results correctly stored at native elemental resolution
- ✅ **Structure consistent**: Uses resolution-first organization (aligned with primary/secondary)

### Tertiary Results Mapping

| Quantity | Native Resolution | Computation | Storage Location | Status |
|----------|------------------|-------------|------------------|--------|
| Section Forces | Gaussian | Direct mapping from stress tensor | `tertiary_results/gaussian/section_forces/` | ✅ Correct |
| Principal Stresses | Gaussian | Eigenvalue decomposition of stress | `tertiary_results/gaussian/principal_stress/` | ✅ Correct |
| Von Mises Stress | Gaussian | `√(½[(σ1-σ2)² + (σ2-σ3)² + (σ3-σ1)²])` | Summary CSV only | ✅ Correct |
| Max Shear Stress | Gaussian | `(σ1 - σ3) / 2` | Summary CSV only | ✅ Correct |
| Failure Index | Gaussian | `σ_vm / (σ_yield / SF)` | Summary CSV only | ✅ Correct |
| Total Strain Energy | Elemental (integrated) | `∫ w dΩ` via quadrature | `tertiary_results/elemental/total_strain_energy.csv` | ✅ Correct |
| Integrated Section Forces | Elemental (integrated) | Weighted average over element | `tertiary_results/elemental/integrated_section_forces.csv` | ✅ Correct |

**Note**: See `TERTIARY_RESULTS_PHILOSOPHY.md` for detailed analysis.

---

## Summary of Native Resolutions

### By Quantity Type

| Quantity Type | Native Resolution | Can Project To | Can Integrate To |
|--------------|-------------------|----------------|------------------|
| Displacements | Global | Elemental (disassembly) | N/A |
| Reactions | Global | Elemental (disassembly) | N/A |
| Strain | Gaussian | Nodal (interpolation) | Elemental (integration) |
| Stress | Gaussian | Nodal (interpolation) | Elemental (integration) |
| Energy Density | Gaussian | Nodal (interpolation) | Elemental (integration) |
| Total Strain Energy | Elemental | N/A | Global (sum) |
| Section Forces | Gaussian | Nodal (if needed) | Elemental (integrated) |
| Principal Stresses | Gaussian | Nodal (if needed) | N/A |
| Von Mises Stress | Gaussian | Nodal (if needed) | N/A |
| Total Strain Energy (Tertiary) | Elemental | N/A | Global (sum) |
| Integrated Section Forces | Elemental | N/A | N/A |

### By Resolution Level

| Resolution | Native Quantities | Projected Quantities | Integrated Quantities |
|-----------|-------------------|---------------------|---------------------|
| **Global** | U_global, R_global, K_global, F_global | None | None |
| **Elemental** | K_e, F_e (formulation) | U_e, R_e (disassembled) | Total strain energy (secondary), Total strain energy (tertiary), Integrated section forces (from Gaussian) |
| **Gaussian** | ε, σ, w (energy density), Section forces, Principal stresses, Von Mises stress, Max shear stress, Failure index | None | None |
| **Nodal** | None | All (projected from Gaussian) | None |

---

## Implementation Recommendations

### Immediate Fixes

1. **Fix nested directory** in `save_secondary_container.py:46`
2. **Add README files** explaining native vs. projected vs. integrated
3. **Update docstrings** to clarify native resolution for each quantity

### Future Enhancements

1. **Metadata files**: Add JSON/YAML files describing:
   - Native resolution of each quantity
   - Projection method used (if applicable)
   - Integration method used (if applicable)

2. **Validation**: Add checks to ensure:
   - Native results exist before projections
   - Projections use cached shape functions
   - Integrations use proper quadrature

3. **Documentation**: Create visual diagrams showing:
   - Computation flow from native to projected/integrated
   - Data dependencies between resolutions

---

## Conclusion

The current results hierarchy philosophy is **sound** and correctly implements the native resolution principle. The main issues are:

1. ✅ **Philosophy is correct**: Native resolution is properly identified
2. ✅ **Storage is mostly correct**: Quantities stored at their native resolution
3. ❌ **Nested directory bug**: Needs immediate fix
4. ⚠️ **Documentation gap**: Needs clarification of native vs. projected vs. integrated

The structure supports both native results and their transformations (projections and integrations), which aligns with the intended philosophy.

