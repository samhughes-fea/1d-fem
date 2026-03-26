# Formulation and tensor docstrings (checklist)

Use this checklist for **element modules** (`*_3D.py`) and **`utilities/`** (`B_matrix.py`, `D_matrix.py`, `interpolate_loads.py`, `shape_functions.py`). Goal: **generalisability** (any theory in the factory) and **completeness** (shapes and assumptions explicit).

## Element class or module docstring

1. **Identity:** Number of nodes, DOF per node, ordering of `U_e` (match [`Element1DBase`](../../pre_processing/element_library/element_1D_base.py) / job convention).
2. **Shapes:** `K_e`, `F_e`; per Gauss point `B`, `D`, strain vector `ε` (or `E`), stress resultant packing.
3. **Kinematics:** Strain definitions and **which axis / frame** (straight local `x`, curved `s`, etc.).
4. **Constitutive:** What `D` contains (EA, EI, GA, GJ, κ, …); **integration** (full / reduced / selective).
5. **Limits:** Reductions (e.g. κ₀→0); links to [`docs/proofs/`](../../docs/proofs) or [`docs/element_library/`](../../docs/element_library).
6. **Public methods:** Return types (`ElementObject`, `MassObject`, …) and what `element_stiffness_matrix`, `element_force_vector`, `element_mass_matrix`, `tangent_stiffness_matrix` compute.

## `B_matrix` / `D_matrix` / loads / shapes utilities

- Input/output **array shapes** and **index meaning** (row *i* of `B` = which strain measure).
- Natural vs physical form (`B̃` vs `B`) and **Jacobian** `dx/dξ`, `detJ`.
- For loads: how `q` maps to `F_e` (`∫ Nᵀ q detJ dξ`).

## Nonlinear (Total Lagrangian family)

- State **linear** modules composed (`D`, shape operator, load interpolation).
- Document **green_lagrange**, **stress resultant**, **geometric stiffness** roles and **TL assumptions** (moderate rotation, etc.).
