# Bar and Truss 3D formulation

## Stiffness

- **Bar-3D** and **Truss-3D** assemble the element stiffness as **K_e = ∫ Bᵀ D B dx** (same as beam formulations).
- Integration uses **Gauss–Legendre** quadrature. Default quadrature order is 1 (exact for constant B); it can be set explicitly or derived from `integration_orders` (axial, torsion; truss also uses shear orders) when not provided.

## Strain and stress components

- **Bar:** 2 components (axial strain, torsion; stress conjugate: N_axial, T).
- **Truss:** 3 components (axial, transverse shear, torsion; stress conjugate: N_axial, V_trans, T).
- **Beam (EB / Timoshenko / Levinson):** 6 components (ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x and resultants N, Vy, Vz, T, My, Mz).

Post-processing (nodal projection, CSV export, section forces, principal stress) supports variable-length strain/stress: 2 (bar), 3 (truss), or 6 (beam). Section force and principal-stress tertiary results are only defined for 6-component beam stress; bar/truss get placeholder outputs for compatibility.

## Optional formulation features (implemented)

- B2 shape-function coefficients for consistency with beam formulation cache.
- Point loads via shape-function operator (all 6 load components).
- Distributed loads via multi-GP quadrature and `LoadInterpolationOperator`.
- Public operator names: `shape_function_operator`, `strain_displacement_operator`, `material_stiffness_operator`.

## Utilities

Bar and truss each have a `utilities` package: shape functions, B-matrix, D-matrix, load interpolation, and local frame (direction cosines; truss adds transverse). The L-matrix (local–global transformation) is no longer stored on the element; stiffness is assembled from B and D only. Reference L and K_local can be built from `build_L_matrix_*` in utilities for tests (e.g. K_e equivalence vs Lᵀ K_local L).
