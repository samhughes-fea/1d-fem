# Tertiary Results Philosophy Documentation

## Executive Summary

Tertiary results are **highly derived engineering quantities** computed from secondary results (stress/strain at Gauss points). They represent design-critical quantities used for verification, failure analysis, and engineering interpretation.

## Core Philosophy

### Tertiary Results Definition

Tertiary results are **post-processed engineering quantities** that:
1. Are computed from secondary results (stress/strain at Gauss points)
2. Represent design-critical engineering quantities
3. Are used for failure analysis, design verification, and engineering interpretation
4. May exist at multiple resolutions (Gaussian, elemental)

### Native Resolution Principle

Unlike primary and secondary results, tertiary results have **multiple native resolutions** depending on the quantity:
- **Gaussian resolution**: Section forces, principal stresses, Von Mises stress (computed directly from stress at Gauss points)
- **Elemental resolution**: Integrated quantities (total strain energy, integrated section forces)

---

## Terminology Definitions

### Tertiary Results Categories

#### 1. Stress-Derived Quantities (Gaussian Resolution)
- **Section Forces** `[N, Vy, Vz, T, My, Mz]`: Stress resultants at Gauss points
  - Native Resolution: **Gaussian** (computed from stress at Gauss points)
  - Computation: Direct mapping from stress tensor components
  - Shape: `List[element] -> List[gauss_point] -> np.ndarray(6,)`

- **Principal Stresses** `[σ1, σ2, σ3]`: Eigenvalues of stress tensor
  - Native Resolution: **Gaussian** (computed from stress at Gauss points)
  - Computation: Eigenvalue decomposition of stress tensor
  - Shape: `List[element] -> List[gauss_point] -> np.ndarray(3,)`

- **Von Mises Stress**: Equivalent stress for failure criteria
  - Native Resolution: **Gaussian** (computed from stress at Gauss points)
  - Computation: `σ_vm = √(½[(σ1-σ2)² + (σ2-σ3)² + (σ3-σ1)²])`
  - Shape: `List[element] -> List[gauss_point] -> float`

- **Maximum Shear Stress**: Maximum shear stress component
  - Native Resolution: **Gaussian** (computed from stress at Gauss points)
  - Computation: `τ_max = (σ1 - σ3) / 2`
  - Shape: `List[element] -> List[gauss_point] -> float`

- **Failure Index**: Material failure criterion
  - Native Resolution: **Gaussian** (computed from Von Mises and yield strength)
  - Computation: `FI = σ_vm / (σ_yield / SF)`
  - Shape: `List[element] -> List[gauss_point] -> float`

#### 2. Integrated Elemental Quantities (Elemental Resolution)
- **Total Strain Energy**: Energy stored per element
  - Native Resolution: **Elemental** (integrated from Gaussian energy density)
  - Computation: `U_e = ∫_Ω w(x) dΩ ≈ ∑_g w_g ⋅ w(x_g) ⋅ |J(x_g)|`
  - Shape: `List[element] -> float`
  - Units: Joules (J)

- **Integrated Section Forces**: Average section forces per element
  - Native Resolution: **Elemental** (integrated/averaged from Gaussian)
  - Computation: Weighted average over element length
  - Shape: `List[element] -> np.ndarray(6,)`
  - Units: `[N, N, N, N⋅m, N⋅m, N⋅m]`

---

## Current Directory Structure Analysis

### Tertiary Results Structure (Resolution-First)

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
- ✅ **Structure consistent**: Uses resolution-first organization (aligned with primary/secondary results)

---

## Computation-to-Storage Mapping

### Tertiary Results

| Quantity | Native Resolution | Computation | Storage Location | Status |
|----------|------------------|-------------|------------------|--------|
| Section Forces | Gaussian | Direct mapping from stress tensor | `tertiary_results/gaussian/section_forces/` | ✅ Correct |
| Principal Stresses | Gaussian | Eigenvalue decomposition of stress | `tertiary_results/gaussian/principal_stress/` | ✅ Correct |
| Von Mises Stress | Gaussian | `√(½[(σ1-σ2)² + (σ2-σ3)² + (σ3-σ1)²])` | Summary CSV only | ✅ Correct |
| Max Shear Stress | Gaussian | `(σ1 - σ3) / 2` | Summary CSV only | ✅ Correct |
| Failure Index | Gaussian | `σ_vm / (σ_yield / SF)` | Summary CSV only | ✅ Correct |
| Total Strain Energy | Elemental (integrated) | `∫ w dΩ` via quadrature | `tertiary_results/elemental/total_strain_energy.csv` | ✅ Correct |
| Integrated Section Forces | Elemental (integrated) | Weighted average over element | `tertiary_results/elemental/integrated_section_forces.csv` | ✅ Correct |

---

## Key Observations

### 1. Multiple Native Resolutions

**Observation**: Tertiary results have quantities at two native resolutions:
- **Gaussian**: Section forces, principal stresses, Von Mises, etc. (computed from stress at Gauss points)
- **Elemental**: Total strain energy, integrated section forces (integrated from Gaussian)

**Status**: ✅ **Correct** - Each quantity is stored at its native resolution

### 2. Derived vs. Native Quantities

**Observation**: Some tertiary quantities are stored individually (section forces, principal stresses), while others are only in summary (Von Mises, max shear, failure index).

**Rationale**: 
- Section forces and principal stresses are **core engineering quantities** used for design
- Von Mises, max shear, and failure index are **derived metrics** primarily used for summary/analysis

**Status**: ✅ **Acceptable** - Design choice to reduce file count while preserving critical data

### 3. Integration Philosophy

**Observation**: Total strain energy and integrated section forces are integrated from Gaussian quantities.

**Formula**:
- Total strain energy: `U_e = ∫_Ω w(x) dΩ ≈ ∑_g w_g ⋅ w(x_g) ⋅ |J(x_g)|`
- Integrated section forces: Weighted average over element length

**Status**: ✅ **Correct** - Properly integrated using Gauss quadrature

### 4. Directory Structure Consistency

**Comparison with Secondary Results**:
- Secondary: `secondary_results/gaussian/`, `secondary_results/nodal/`, `secondary_results/elemental/`
- Tertiary: `tertiary_results/gaussian/`, `tertiary_results/elemental/`

**Observation**: Tertiary now uses **resolution-first** organization (gaussian, elemental), consistent with secondary results.

**Status**: ✅ **Consistent** - Aligned with primary/secondary results organization pattern

---

## Identified Issues

### 1. Directory Structure Inconsistency ✅ RESOLVED

**Issue**: Tertiary results previously used quantity-first organization, while secondary results used resolution-first organization.

**Previous Structure** (before refactor):
```
tertiary_results/
├── section_forces/          # Quantity-first
├── principal_stress/        # Quantity-first
└── elemental/               # Resolution-first (mixed!)
```

**Current Structure** (after refactor):
```
tertiary_results/
├── gaussian/                # Resolution-first
│   ├── section_forces/
│   └── principal_stress/
└── elemental/               # Resolution-first
```

**Secondary Structure** (for comparison):
```
secondary_results/
├── gaussian/                # Resolution-first
│   ├── strain/
│   ├── stress/
│   └── energy_density/
├── nodal/                   # Resolution-first
└── elemental/               # Resolution-first
```

**Resolution**: ✅ **RESOLVED** - Structure has been refactored to resolution-first:
```
tertiary_results/
├── gaussian/                # Resolution-first
│   ├── section_forces/
│   └── principal_stress/
└── elemental/               # Resolution-first
    ├── total_strain_energy.csv
    └── integrated_section_forces.csv
```

**Status**: ✅ **Consistent** - Now aligned with secondary results philosophy

### 2. Missing Resolution-Level Organization ✅ RESOLVED

**Issue**: Gaussian-level tertiary results were previously organized by quantity, not by resolution.

**Previous**: `tertiary_results/section_forces/` (quantity-first)
**Current**: `tertiary_results/gaussian/section_forces/` (resolution-first)

**Resolution**: Structure has been refactored to use resolution-first organization, consistent with secondary results.

**Status**: ✅ **Resolved** - Structure now follows resolution-first organization

### 3. Summary-Only Quantities

**Issue**: Von Mises stress, max shear stress, and failure index are only saved in summary CSV, not as individual files.

**Current**: Only in `tertiary_summary.csv`
**Alternative**: Could save as individual files per element (like section forces)

**Recommendation**: 
- **Keep current approach** if these are primarily used for summary/analysis
- **Add individual files** if detailed per-element analysis is needed

**Status**: ✅ **Acceptable** - Design choice, not a bug

---

## Recommended Directory Structure

### ✅ IMPLEMENTED: Resolution-First Structure (Option B)

**Status**: This structure has been implemented and is now the current structure.

**Pros**:
- ✅ Consistent with secondary results philosophy
- ✅ Clear resolution distinction
- ✅ Better for programmatic access
- ✅ Aligned with primary/secondary results organization

**Cons**:
- Slightly more nesting (acceptable trade-off for consistency)

**Current Structure**:
```
tertiary_results/
├── gaussian/                    # Resolution-first
│   ├── section_forces/
│   │   └── section_forces_elem_*.csv
│   └── principal_stress/
│       └── principal_stress_elem_*.csv
├── elemental/                   # Resolution-first
│   ├── total_strain_energy.csv
│   └── integrated_section_forces.csv
└── tertiary_summary.csv
```

### Option A: Quantity-First (Previous Structure - Deprecated)

**Note**: This structure was previously used but has been replaced with the resolution-first structure for consistency.

**Previous Structure** (for reference):
```
tertiary_results/
├── section_forces/              # Gaussian resolution (quantity-first)
│   └── section_forces_elem_*.csv
├── principal_stress/            # Gaussian resolution (quantity-first)
│   └── principal_stress_elem_*.csv
├── elemental/                   # Elemental resolution
│   ├── total_strain_energy.csv
│   └── integrated_section_forces.csv
└── tertiary_summary.csv
```

**Why it was changed**:
- Inconsistent with secondary results
- Mixed organization (quantity-first for Gaussian, resolution-first for elemental)
- Less intuitive for engineering users

```
tertiary_results/
├── gaussian/                    # Native: Gaussian resolution
│   ├── section_forces/
│   │   └── section_forces_elem_*.csv
│   ├── principal_stress/
│   │   └── principal_stress_elem_*.csv
│   ├── von_mises_stress/        # Optional: if saving individually
│   └── max_shear_stress/        # Optional: if saving individually
├── elemental/                   # Native: Elemental resolution
│   ├── total_strain_energy.csv
│   └── integrated_section_forces.csv
└── tertiary_summary.csv
```

**Recommendation**: **Option A (Current)** - Keep quantity-first organization for tertiary results because:
1. Tertiary results are engineering quantities, not field quantities
2. Users typically search by quantity (section forces, principal stress) rather than resolution
3. Simpler structure is more intuitive
4. The inconsistency is acceptable given the different nature of tertiary vs. secondary results

---

## Summary of Native Resolutions

### By Quantity Type

| Quantity Type | Native Resolution | Can Project To | Can Integrate To |
|--------------|-------------------|----------------|------------------|
| Section Forces | Gaussian | Nodal (if needed) | Elemental (integrated) |
| Principal Stresses | Gaussian | Nodal (if needed) | N/A |
| Von Mises Stress | Gaussian | Nodal (if needed) | N/A |
| Max Shear Stress | Gaussian | Nodal (if needed) | N/A |
| Failure Index | Gaussian | Nodal (if needed) | N/A |
| Total Strain Energy | Elemental | N/A | Global (sum) |
| Integrated Section Forces | Elemental | N/A | N/A |

### By Resolution Level

| Resolution | Native Quantities | Derived Quantities | Integrated Quantities |
|-----------|-------------------|-------------------|---------------------|
| **Gaussian** | Section forces, Principal stresses | Von Mises, Max shear, Failure index | None |
| **Elemental** | Total strain energy, Integrated section forces | None | None |

---

## Implementation Recommendations

### Immediate Actions

1. ✅ **Documentation**: Add this philosophy document
2. ⚠️ **Consistency**: Decide on quantity-first vs. resolution-first (recommend keeping current)
3. ✅ **Clarification**: Update docstrings to clarify native resolution for each quantity

### Future Enhancements

1. **Optional Individual Files**: Consider saving Von Mises, max shear, and failure index as individual files if needed
2. **Nodal Projections**: Add nodal projections of section forces and principal stresses if needed for visualization
3. **Metadata**: Add JSON/YAML files describing computation methods and native resolutions

---

## Conclusion

The tertiary results structure is **functionally correct** and stores quantities at their native resolutions. The main observation is:

1. ✅ **Philosophy is sound**: Quantities stored at native resolution
2. ✅ **Computation is correct**: Proper integration and derivation methods
3. ⚠️ **Structure is inconsistent**: Uses quantity-first vs. secondary's resolution-first
4. ✅ **Design choices are reasonable**: Summary-only quantities acceptable

The current structure prioritizes **usability for engineering analysis** over strict consistency with secondary results, which is a reasonable design choice given the different nature of tertiary results (engineering quantities vs. field quantities).

