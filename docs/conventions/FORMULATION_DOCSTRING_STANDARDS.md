# Formulation and tensor docstrings (checklist)

Use this checklist for **element modules** (`*_3D.py`) and **`utilities/`** (linear: `B_matrix.py`, `D_matrix.py`, `interpolate_loads.py`, `shape_functions.py`; nonlinear TL: `green_lagrange_strain.py`, `stress_resultant.py`, `geometric_stiffness.py`, etc.). Goal: **generalisability** (any theory in the factory) and **completeness** (shapes and assumptions explicit).

## Docstring gold standard (library-wide)

Every element formulation should document **both** of the following in Python docstrings (NumPy / numpydoc style — see below):

1. **Discrete weak form** — Which Gauss sums (or provably equivalent operations) the code uses; natural coordinate **xi** in `[-1, 1]`; **`detJ = dx/dxi`**; selective or reduced integration where used.
2. **Physical meaning** — Voigt ordering of strains and stress resultants (or a pointer to this document **plus** theory-specific rows); what **`D`** contains; local frame (**`x`**, curved **`s`**, κ₀, etc.).

Weak form without Voigt meaning is hard to verify; physics without the stated quadrature pattern is ambiguous about what runs in the loop.

## NumPy docstring structure (Python source)

Follow the [NumPy docstring guide](https://numpydoc.readthedocs.io/en/latest/format.html).

- **Module docstring:** One-line summary; blank line; extended summary using **consistent bullet headings** so every `*_3D.py` file scans the same way: **Tensors**, **Weak forms**, **Kinematics**, **Constitutive**, **Quadrature**, **Public API** (omit headings that do not apply).
- **Class docstring:** Summary line; extended summary (theory + weak-form linkage); **Parameters**, **Attributes**, **Raises** as needed; put long weak-form or implementation detail in **Notes**; **See Also** for `docs/element_library/` or proofs; **Examples** only when they clarify usage.
- **Method docstrings:** **Parameters**, **Returns**, optional **Notes**.

Do not put the only statement of the weak form solely in **Returns** without context in the extended summary or **Notes**.

## Plain-text math in Python docstrings

Docstrings are read in IDEs and `help()`, not compiled as LaTeX. Prefer readable, code-aligned notation:

- Stiffness increment: ``K_e += B.T @ D @ B * w_g * detJ`` (summed over Gauss points `g`).
- Strain and stress: ``eps = B @ U_e``, ``S = D @ eps`` (linear); TL variants use ``E``, ``B_lin``, ``B_nl`` as implemented.
- Loads: ``F_dist += w_g * N.T @ q * detJ``; point loads: ``F_point = N(x_p).T @ P`` (or equivalent at natural coordinate).
- Name the reference coordinate explicitly: **xi in [-1, 1]** (Unicode ξ is optional in prose, not required on the canonical line).
- Use ASCII names in prose (`eps_x`, `kappa_y`, `M_y`) or plain words; avoid large LaTeX display blocks inside triple-quoted Python strings.
- For **D** structure, prefer a short diagonal list or words over big matrices.

The **Markdown tables and LaTeX below** in *this* conventions file remain the reference for documentation generators and reviews; Python docstrings should still follow the plain-text rules above when restating the same facts.

## Continuum and discrete strain conventions (infinitesimal vs Green–Lagrange)

This library’s **linear** 1D beam elements use **infinitesimal strain** (displacement gradients \(\ll 1\)); **Total Lagrangian (TL) nonlinear** beam elements in `nonlinear/.../euler_bernoulli` and `nonlinear/.../timoshenko` use a **Green–Lagrange–type** strain vector in the **reference (initial) configuration** and work-conjugate **2nd Piola–Kirchhoff** resultants in Voigt form. **Other** nonlinear element types (Updated Lagrangian stubs, co-rotational, GEBT placeholders) are **not** required to use the same \((\mathbf{E}, \mathbf{S})\) pair; their docstrings must state their own strain and stress measures.

### Infinitesimal strain (linear track)

- **Strain tensor:** \(\varepsilon_{ij} = \tfrac{1}{2}\bigl(\partial u_i/\partial x_j + \partial u_j/\partial x_i\bigr)\). **Discrete:** \(\mathbf{u} = \mathbf{N}\mathbf{d}\), generalized engineering strain vector \(\boldsymbol{\varepsilon} = \mathbf{B}\mathbf{d}\) (Voigt \((6,)\) for 12-DOF beams per table below).
- **Stress:** Cauchy \(\sigma_{ij}\) (true stress on the current configuration in 3D theory). In the **beam resultants** used here, \(\mathbf{S} = \mathbf{D}\,\boldsymbol{\varepsilon}\) in Voigt (axial, moments, shears, torque) — the same packing as linear `B_matrix` / `D_matrix`.
- **Stiffness:** Material stiffness \(\mathbf{K} = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,\mathrm{d}x\) (Gauss form with \(|J| = L/2\)) is **constant** w.r.t. \(\mathbf{d}\); one linear solve.

### Green–Lagrange strain (TL nonlinear track)

- **Strain:** \(\mathbf{E} = \tfrac{1}{2}(\mathbf{F}^\top\mathbf{F} - \mathbf{I})\) with \(F_{ij} = \delta_{ij} + \partial u_i/\partial X_j\) in 3D; the code uses a **reduced beam vector** \(\mathbf{E}\) (Voigt) with \(\mathbf{E} = \mathbf{E}_\mathrm{lin} + \mathbf{E}_\mathrm{nl}\) from `GreenLagrangeStrainOperator`.
- **Stress:** **2nd Piola–Kirchhoff** \(\mathbf{S}\) conjugate to \(\mathbf{E}\); **St. Venant–Kirchhoff–style** beam reduction \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) with the same \(\mathbf{D}\) structure as linear elasticity but evaluated on \(\mathbf{E}\) (see `stress_resultant`).
- **Tangent:** Schematically \(\mathbf{K}_T(\mathbf{d}) = \mathbf{K}_L + \mathbf{K}_{NL}(\mathbf{d}) + \mathbf{K}_\sigma(\mathbf{d})\). In this repository’s TL Timoshenko implementation: \(\mathbf{K}_T = \mathbf{K}_0 + \mathbf{K}_\delta + \mathbf{K}_\sigma\) where \(\mathbf{K}_0\) is the linear-type material stiffness (selective assembly), \(\mathbf{K}_\delta\) captures the incremental material tangent from \(\mathbf{B}_\mathrm{tot}=\mathbf{B}_\mathrm{lin}+\mathbf{B}_\mathrm{nl}\), and \(\mathbf{K}_\sigma\) is the **geometric** stiffness from current section forces. TL Euler–Bernoulli uses the analogous \(\mathbf{B}_\mathrm{lin}\)/\(\mathbf{B}_\mathrm{nl}\) split and often denotes material tangent \(\mathbf{K}_\mathrm{mat}\) in element docs. **Newton–Raphson** rebuilds \(\mathbf{K}_T\) each iteration.

### Voigt vs tensor concepts (straight 12-DOF beam)

| Voigt row | Symbol | Infinitesimal (linear) meaning | TL nonlinear: Green–Lagrange vector \(\mathbf{E}\) |
|-----------|--------|--------------------------------|---------------------------------------------------|
| 0 | \(\varepsilon_x\) / \(E_x\) | Axial infinitesimal strain | Axial \(E\) with Green–Lagrange quadratics + coupling |
| 1–2 | \(\kappa_y\), \(\kappa_z\) | Bending curvatures | Same rows + nonlinear supplements from `GreenLagrangeStrainOperator` |
| 3–4 | \(\gamma_{xy}\), \(\gamma_{xz}\) | Shear (zero for EB) | EB: zero; Timoshenko: linear part + optional NL shear terms |
| 5 | \(\phi_x\) | Twist rate | Typically linear in \(\theta_x\) |

---

## `utilities/` operator docstrings — governing equations checklist

Every **`utilities/`** operator that participates in strain or stiffness assembly must, in its **module** and **class** docstrings, make the following explicit within the first extended summary (so reviewers see theory without opening the parent element):

### Linear beam (`linear/.../utilities/B_matrix.py`, `D_matrix.py`)

| Operator | Governing equations to state |
|----------|------------------------------|
| **StrainDisplacementOperator** | \(\varepsilon_{ij}\) (infinitesimal); discrete \(\boldsymbol{\varepsilon} = \mathbf{B}\,\mathbf{U}_e\); row packing vs Voigt table above; **reference:** small strain on the undeformed beam axis in local frame. |
| **MaterialStiffnessOperator** | Linear elastic \(\mathbf{S} = \mathbf{D}\,\boldsymbol{\varepsilon}\); \(\mathbf{D}\) entries (EA, EI, GA, GJ, \(\kappa\)); **constant** contribution to \(\mathbf{K} = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,\mathrm{d}x\). |

### Nonlinear TL beam (`nonlinear/.../utilities/green_lagrange_strain.py`, `stress_resultant.py`, `geometric_stiffness.py`)

| Operator | Governing equations to state |
|----------|------------------------------|
| **GreenLagrangeStrainOperator** | \(\mathbf{E} = \mathbf{E}_\mathrm{lin} + \mathbf{E}_\mathrm{nl}\); \(\mathbf{B}_\mathrm{lin}\), \(\mathbf{B}_\mathrm{nl}\); reference configuration \(X\); relation to `strain_linear_part` / `strain_nonlinear_part` / `nonlinear_strain_displacement_gradient`. |
| **StressResultantOperator** | \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) (2PK-style resultants in Voigt); which components drive \(\mathbf{K}_\sigma\) (typically \(N\), \(M_y\), \(M_z\)). |
| **GeometricStiffnessOperator** | \(\mathbf{K}_\sigma\) as Gauss sum of stress-weighted products of shape-function derivatives; **not** a substitute for \(\mathbf{K}_L\) / \(\mathbf{K}_0\). |

**Naming:** In code and docstrings, prefer **`E`** / **`E_lin`** / **`E_nl`** for Green–Lagrange, **`eps`** or **`epsilon`** for infinitesimal Voigt strain, and avoid ambiguous bare `strain` without context.

## Weak-form assembly (required)

Every element must assemble tensors from a **weak form** discretized with **Gauss quadrature** on the reference map \(\xi \in [-1,1]\) with Jacobian **`|J| = dx/dξ`** (chord map unless the theory states otherwise). **Do not** use closed-form element matrices for stiffness, geometric tangent, mass, or distributed loads unless that matrix is **provably identical** to the Gauss sum for that operator (e.g. constant **`B`** with one-point exact rule) and the docstring **still states the integral** the code implements.

**Per quantity (state each that applies in a “Weak forms” block):**

1. **Stiffness:** \(K_e \mathrel{+}= B^\top D\, B\, w_g\, |J|\) at each Gauss point, or a **documented sum** of such terms for **selective / multi-rule** integration (different point sets or extracted sub-blocks are still Gauss weak-form sums).
2. **Distributed loads:** \(F_{\mathrm{dist}} = \int_{-1}^{1} N^\top q\, |J|\, d\xi \approx \sum_g w_g\, N(\xi_g)^\top q(\xi_g)\, |J|\).
3. **Point loads:** Work-equivalent \(F_{\mathrm{point}} = N(\xi_p)^\top P\) (not a domain integral; still document explicitly under **Loads**).
4. **Consistent mass:** \(M_e \mathrel{+}= \tfrac{1}{2}(\mu_i+\mu_j)\, (N_i\cdot N_j)\, w_g |J|\) (or the theory’s equivalent) summed over Gauss points; define \(\mu\) per DOF.
5. **Nonlinear internal force:** \(F_{\mathrm{int}} \mathrel{+}= B^\top S\, w_g |J|\) with the **`B`** used in code (e.g. linearized `B_lin` for TL).
6. **Material tangent (TL):** \(K_{\mathrm{mat}} \mathrel{+}= (B_{\mathrm{lin}}+B_{\mathrm{nl}})^\top D\, (B_{\mathrm{lin}}+B_{\mathrm{nl}})\, w_g |J|\) when that is the implemented operator.
7. **Geometric tangent (TL):** \(K_\sigma \mathrel{+}= \sum_g w_g |J|\) times **stress-weighted** products of shape-function derivatives (e.g. \((N + M_z/L)\) times \(\partial h/\partial x\) outer products for EB bending planes — see nonlinear `geometric_stiffness` implementation). **No** undocumented standalone beam-column templates.

**Plain-text lines for Python docstrings (same rules, copy-adapt):**

1. Stiffness: ``K_e += B.T @ D @ B * w_g * detJ`` per Gauss point (or documented split for selective integration).
2. Distributed loads: ``F_dist = sum_g w_g * N.T @ q * detJ`` (with `q` at Gauss `xi_g`).
3. Point loads: ``F_point = N.T @ P`` evaluated at the load station (work-equivalent).
4. Mass: consistent mass as sum over `g` of ``0.5 * (mu_i + mu_j) * dot(N_i, N_j) * w_g * detJ`` (or theory equivalent); define `mu` per DOF in the docstring.
5. Nonlinear internal force: ``F_int += B.T @ S * w_g * detJ`` with the **`B`** the code uses (e.g. `B_lin`).
6. Material tangent (TL): ``K_mat += B_tot.T @ D @ B_tot * w_g * detJ`` with ``B_tot = B_lin + B_nl`` when that matches implementation.
7. Geometric tangent (TL): ``K_sigma += (...) * w_g * detJ`` as stress-weighted products of shape derivatives — state which section forces and which displacement derivatives enter **your** file; no undocumented beam-column templates.

---

## Element class or module docstring

1. **Identity:** Number of nodes, DOF per node, ordering of `U_e` (match [`Element1DBase`](../../pre_processing/element_library/element_1D_base.py) / job convention).
2. **Weak forms:** Bulleted list of the integrals **this file implements** (copy from the list above; omit methods the element does not expose).
3. **Shapes:** `K_e`, `F_e`; per Gauss point `B`, `D`, strain vector `ε` (or `E`), stress resultant packing.
4. **Kinematics:** Strain definitions and **which axis / frame** (straight local `x`, curved `s`, etc.).
5. **Constitutive:** What `D` contains (EA, EI, GA, GJ, κ, …); **integration** (full / reduced / selective).
6. **Limits:** Reductions (e.g. κ₀→0); links to [`docs/proofs/`](../../docs/proofs) or [`docs/element_library/`](../../docs/element_library).
7. **Public methods:** Return types (`ElementObject`, `MassObject`, …) and what `element_stiffness_matrix`, `element_force_vector`, `element_mass_matrix`, `tangent_stiffness_matrix` compute.

## `B_matrix` / `D_matrix` / loads / shapes utilities

- Input/output **array shapes** and **index meaning** (row *i* of `B` = which strain measure).
- Natural vs physical form (`B̃` vs `B`) and **Jacobian** `dx/dξ`, `detJ`.
- For loads: how `q` maps to `F_e` (`∫ Nᵀ q detJ dξ`).

## Nonlinear Total Lagrangian — composition and utility docstrings

Nonlinear elements assemble **internal force** and **tangents** in the parent `*_3D` module by **summing** per–Gauss-point contributions with `w_g` and `detJ`. TL **utilities are not standalone physics packages**: each implements tensors the element loops over.

**Typical pipeline (Euler–Bernoulli / Timoshenko TL in this repo):**

1. **`green_lagrange_strain`** — Green–Lagrange strain `E` (or stacked form used in code) and strain–displacement operators `B_lin`, `B_nl` at a Gauss point; relationship to nodal `U_e`.
2. **`stress_resultant`** — `S` from `D` and strain (`S = D @ E` or equivalent); which entries feed the **material** tangent vs **geometric** stiffness; row order must match linear `D` / Voigt convention.
3. **`geometric_stiffness`** — `K_sigma` as a **Gauss sum** of stress-weighted products of shape-function derivatives (not an undocumented closed-form beam-column matrix); which axial / bending resultants and which slopes enter the implementation.

**Parent element module** must state, in **Weak forms**, at least: ``F_int += B.T @ S * w_g * detJ`` (with identified `B`); material and geometric tangent lines matching code; same load and mass rules as linear if exposed.

**Per-utility docstrings (TL):** Use NumPy sections. **Notes** should include **implementation assumptions** (e.g. moderate rotation, which curvature measure, values held constant vs evaluated at Gauss points). **See Also** [total_lagrangian_beam_formulation.md](../element_library/total_lagrangian_beam_formulation.md) or proofs when relevant. If long-form theory docs lag the code, the docstring must still describe **what the code does**. Optional **Limitations** bullet in **Notes** for provisional behaviour; avoid public `TODO` unless the project explicitly wants them.

**GEBT shear** and other families may use a different utility set; apply the same principles: list operators, Gauss composition, assumptions.

**Linear modules reused by TL:** Still document which linear `D_matrix`, `B_matrix` / shapes, and load operators are composed.

---

## Tensor and Voigt contract (reference for all beam elements)

Use this table as the **single source of truth** for straight 12-DOF beam elements unless the theory explicitly extends it (e.g. warping, bar/truss reductions). Element and `B_matrix` / `D_matrix` docstrings should say “Voigt order per FORMULATION_DOCSTRING_STANDARDS” or repeat this table once per family.

### Outer sizes as comparison frame

When multiple theories share the same **outer** tensor shapes (e.g. `B` `(6, 12)`, `U_e` `(12,)` for EB, Timoshenko, Levinson), use **Contract** / **Diff** notes in utility or class docstrings so reviewers see **row activity, order, zeros, and extensions** at a glance:

- **Baseline:** simpler models often have **zero rows** in `B` or `D` (e.g. EB shear strains constitutively zero).
- **Progression:** higher-order or shear-deformable theories **fill** those rows with kinematic terms and matching stiffness (e.g. `kappa*G*A` vs `G*A`); warping adds DOFs and an extra strain row; curved or TL formulations add parameters (`kappa0`, `B_nl`) while keeping the packing explicit.

Extensions such as `(14,) U_e` and `(7, 14) B` should state how the **first 12 DOFs and 6 strain rows** embed the linear baseline.

**Local frame:** Straight elements use local **`x`** along the chord from node 1 → node 2 (same isoparametric map \(x(\xi)\) on the chord). Reference curvature for initially curved geometry is supplied via **`precurvature.txt`** and straight beam implementations — state conventions in the element docstring when relevant.

**Nodal displacement vector `U_e`:** Length **`2 × dof_per_node`**, **node-major**: all DOFs at node 1, then node 2. For **`dof_per_node = 6`**:

| Index block | Physical DOFs (in order) |
|-------------|---------------------------|
| 0–5 (node 1) | \(u_x, u_y, u_z, \theta_x, \theta_y, \theta_z\) |
| 6–11 (node 2) | same order |

Match [`Element1DBase`](../../pre_processing/element_library/element_1D_base.py) and job/global assembly conventions.

**Engineering strain vector `ε` (6,) — Voigt row order**

| Row | Symbol | Meaning |
|-----|--------|---------|
| 0 | ε_x | Axial strain |
| 1 | κ_y | Bending curvature about **y** (x–z plane); definition depends on theory (EB: from \(w\); Timoshenko: from \(\theta_y\)) |
| 2 | κ_z | Bending curvature about **z** (x–y plane) |
| 3 | γ_xy | Transverse shear xy |
| 4 | γ_xz | Transverse shear xz |
| 5 | φ_x | Twist rate \(d\theta_x/dx\) |

**Stress resultants `S` (6,) paired with `ε`:** \(S = [N, M_y, M_z, V_y, V_z, T]^T\) with **`S = D @ ε`** in linear problems. (EB: rows 3–4 of `D` are zero so \(V_y=V_z=0\) from constitutive law; use equilibrium for shear if needed.)

**Matrices per Gauss point (12-DOF beam):**

- **`B`:** `(6, 12)` with **`ε = B @ U_e`** in physical **`x`** where the assembly uses `physical_coordinate_form`.
- **`D`:** `(6, 6)` symmetric (Timoshenko: shear diagonal uses κGA; Levinson/Reddy: GA without κ in those rows).
- **`N` (shape functions):** From `ShapeFunctionOperator.natural_coordinate_form(ξ)`: for one point, **`N`** shape **`(12, 6)`** — row `a` is DOF `a`, column `c` is displacement component **`c ∈ {0..5}`** = \((u_x, u_y, u_z, \theta_x, \theta_y, \theta_z)\). Batched over Gauss points: **`(n_gp, 12, 6)`** (see [API_STANDARDS.md](API_STANDARDS.md)).

**Weak form (stiffness):** For each Gauss point \(\xi_g\) with weight \(w_g\),

\[
K_e \mathrel{+}= B^T D\, B\, w_g\, |J|, \quad |J| = dx/d\xi = L/2
\]

**Distributed loads:** Equivalent nodal vector contribution

\[
F_{\mathrm{dist}} = \int_{-1}^{1} N^T q(\xi)\, |J|\, d\xi
\]

with **`q`** six components per station matching load columns (force and moment per API_STANDARDS).

**Consistent mass (when implemented):** \(M_e \mathrel{+}= \tfrac{1}{2}(\mu_i+\mu_j)\, (N_i\cdot N_j)\, w_g |J|\) with per-DOF weights \(\mu\) from \(\rho A\), \(\rho J_t\), \(\rho I_y\), \(\rho I_z\) on translation / torsion / bending-rotation DOFs (see linear Timoshenko implementation).

**Reddy vs Levinson:** Prefer **“same Voigt layout as Levinson; differences: …”** in Reddy docstrings to avoid duplicated tables that drift.

---

## Module + class docstring template (copy headings)

**Module one-liner (top of file):** Theory name, nodes, `U_e` length, `K_e`/`F_e` shapes, strain length, one integration/Jacobian fact.

**Module extended summary:** After the one-liner, use the bullet headings **Tensors**, **Weak forms**, **Kinematics**, **Constitutive**, **Quadrature**, **Public API** (plain-text weak-form lines per *Plain-text math* above).

**Class docstring — use these headings (as bullets or short paragraphs, aligned with numpydoc):**

1. **Identity** — 2 nodes; `dof_per_node`; `U_e` ordering (reference table above if 6/node).
2. **Weak forms** — Same integrals as module **Weak forms** block (may be one line “as module”).
3. **Tensors** — `K_e`, `F_e`, `B`, `D`, `ε` (and `S`) shapes; warping/nonlinear extensions if any.
4. **Kinematics** — Definitions of non-zero strain rows; frame (`x`, `s`, κ₀).
5. **Constitutive** — What enters `D`; selective/reduced quadrature if not full Gauss on all terms.
6. **Quadrature** — `|J| = L/2`; orders or multi-rule summary (Timoshenko).
7. **Public API** — One line each: `element_stiffness_matrix` → `ElementObject`; `element_force_vector` → `ForceObject`; `element_mass_matrix` → `MassObject` (if implemented); nonlinear: `tangent_stiffness_matrix` / internal force.
8. **Limits / links** — Reductions (κ₀→0, EB limit); optional `docs/element_library/`, `docs/proofs/`.

**Out of scope for this checklist:** Parsers, `run_job`, and non-element modules unless they expose the same tensors.
