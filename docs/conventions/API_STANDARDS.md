# API standards — 1D FEM (beams / frames)

Applied conventions for this repository: **weak-form assembly**, **tensor layout**, **linear `utilities/` modules**, **nonlinear TL extensions**, and minimal Python style rules.

**Element and utility docstrings:** Use **NumPy / numpydoc** section layout and **plain-text, code-like math** in Python source (e.g. ``K_e += B.T @ D @ B * w_g * detJ``), not LaTeX-only exposition. Full checklist, Voigt tables, TL utility roles, and copy-paste weak-form lines: [FORMULATION_DOCSTRING_STANDARDS.md](FORMULATION_DOCSTRING_STANDARDS.md) (*Docstring gold standard*, *Plain-text math*, *Nonlinear Total Lagrangian*).

---

## 1. Stiffness: B, D, J

At each Gauss point (natural coordinate ξ, physical `x` along the element):

- **Jacobian:** `detJ = dx/dξ` (e.g. `L/2` for a straight 2-node map). Element code exposes this as `jacobian` / `jacobian_determinant`.
- **Strain–displacement:** `ε = B U_e` with `B` in **physical** `x` where required. [`StrainDisplacementOperator`](../../pre_processing/element_library/linear/timoshenko/utilities/B_matrix.py) provides `natural_coordinate_form` (B̃) and `physical_coordinate_form` (B).
- **Constitutive:** Stress resultants from `S = D @ ε` with `D` from [`MaterialStiffnessOperator`](../../pre_processing/element_library/linear/timoshenko/utilities/D_matrix.py). Row ordering of `ε` and `S` must match (e.g. Timoshenko-style six components).
- **Quadrature:** `K_e += Bᵀ D B w_g detJ` (sum over Gauss points). **Selective integration** (e.g. Timoshenko shear) is still the same pattern with different rules on **sub-rows** of `B`.

**Generality:** Warping or higher-order theories may use **7×14** (or other) `B`; document shapes in the element docstring — do not assume `(6, 12)` only in prose.

---

## 2. Loads: N and f

- **Distributed:** `F_dist = ∫ Nᵀ q detJ dξ` with `q` interpolated at Gauss `x` ([`LoadInterpolationOperator`](../../pre_processing/element_library/linear/timoshenko/utilities/interpolate_loads.py)).
- **Point:** `F_point = N(x_p)ᵀ P`.
- Shape functions `N`, `dN_dξ`, `d2N_dξ2` come from [`get_shape_function_operator`](../../pre_processing/element_library/shape_function_registry.py) or theory-local [`shape_functions.py`](../../pre_processing/element_library/linear/euler_bernoulli/utilities/shape_functions.py).

---

## 3. Linear elements: `utilities/` layout

Under each `pre_processing/element_library/linear/<theory>/utilities/`:

| File | Role |
|------|------|
| `B_matrix.py` | Strain–displacement operator(s) |
| `D_matrix.py` | Section/material stiffness |
| `interpolate_loads.py` | Distributed loads at Gauss points |
| `shape_functions.py` | `N(ξ)` and derivatives |

**Re-exports:** If logic is shared (e.g. curved Timoshenko reusing straight Timoshenko `D`), the file must still exist and the module docstring must state the dependency.

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

- **Entry point:** [`run_job.py`](../../workflow_orchestrator/run_job.py) builds elements and dispatches by `simulation_settings["type"]` (e.g. `static`, `static_nonlinear`, `modal`, `dynamic`).
- **Runners** use **`run()`** as the main public step (not a generic `execute()` facade).
- **Settings:** [`simulation_settings_parser.py`](../../pre_processing/parsing/simulation_settings_parser.py); optional `[Newton]` for nonlinear static.

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
