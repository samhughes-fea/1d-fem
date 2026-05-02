# Euler–Bernoulli beam with Vlasov warping (7 DOF / node)

This document is a formulation reference for the 2-node 3D **Euler–Bernoulli** beam extended with **Vlasov-type** non-uniform torsion: two extra nodal degrees of freedom carry a **warping intensity** \(\chi\). The public implementation is [`LinearEulerBernoulliBeamElement3D`](../euler_bernoulli/linear_euler_bernoulli_3D.py) with `[warping]` in `element.txt` (this folder documents operators shared with that path). Rows **0–5** of strain/stress and columns **0–11** of **B** match the linear EB element documented in [linear Euler–Bernoulli (3D)](../euler_bernoulli/README.md). Voigt ordering for those rows follows [`docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`](../../../../../docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md).

---

## 1. Element identity

| Item | Value |
|------|--------|
| Class | `LinearEulerBernoulliBeamElement3D` (+ `[warping]`) |
| Nodes | 2 |
| DOFs | **7** per node → **14** element DOFs |
| Kinematics | ZOSDT Euler–Bernoulli + **non-uniform torsion** (warping strain row) |
| Stiffness | Geometrically linear |

**Embedding.** The first **12** displacement components and the first **six** strain rows are the same operators as in [standard EB](../euler_bernoulli/README.md). Indices **12–13** are warping intensities \(\chi^1, \chi^2\).

---

## 2. Reference configuration and coordinate mapping

Identical to linear EB: straight chord of length \(L\), isoparametric \(\xi \in [-1,1]\), \(\det J = L/2\), chain rules \(\partial\xi/\partial x = 2/L\), \(\partial^2\xi/\partial x^2 = 4/L^2\). See §2 of the [EB README](../euler_bernoulli/README.md).

---

## 3. Degrees of freedom and displacement vector

**Stiffness column order** (matches **B** columns):

\[
\mathbf{U}_e = [u_x^1,\, u_y^1,\, u_z^1,\, \theta_x^1,\, \theta_y^1,\, \theta_z^1,\,
u_x^2,\, u_y^2,\, u_z^2,\, \theta_x^2,\, \theta_y^2,\, \theta_z^2,\, \chi^1,\, \chi^2]^\top \in \mathbb{R}^{14}.
\]

| Indices | Content |
|---------|---------|
| 0–11 | Standard 12-DOF EB packing (same as [EB README](../euler_bernoulli/README.md) §3) |
| 12–13 | Nodal warping intensities \(\chi^1, \chi^2\) (sectorial / warping amplitude; implementation uses them as the discrete warping DOFs paired with \(E\Gamma\) stiffness) |

---

## 4. Vlasov warping theory — physical context

For **open** thin-walled sections (I, C, Z, channel, …), **uniform** St. Venant torsion alone does not capture **non-uniform** twist: longitudinal fibers stretch when warping is **restrained**, producing **bimoments** and **normal stresses** along the beam axis in addition to the shear flow from pure torque.

In this 1D discrete model:

| Mechanism | Strain row | Typical stiffness | Dominant modulus |
|-----------|------------|---------------------|------------------|
| St. Venant (uniform twist rate) | Row 5: \(\phi_x = \partial\theta_x/\partial x\) | \(G J_t\) on \(D_{55}\) | **G** (shear-dominated torsion) |
| Warping / non-uniform torsion | Row 6: \(\phi_x' = \partial\theta_x/\partial x + \partial\chi/\partial x\) | \(E \Gamma\) on \(D_{66}\) | **E** (longitudinal **normal** stresses from restrained warping — an axial-stretch-type mode in the warping sense) |

\(\Gamma\) is the **warping constant** (sectorial second moment), typically **m\(^6\)**. The implementation reads \(\Gamma\) from `section_array[9]` when `len(section_array) >= 10`; otherwise \(\Gamma = 0\) (no warping stiffness, but the element remains 14×14).

---

## 5. Shape functions

- **Rows 0–11 (stiffness path):** identical registry tensors as EB: \(\mathbf{N}\), \(\partial\mathbf{N}/\partial\xi\), \(\partial^2\mathbf{N}/\partial\xi^2\) with batch shape **(n_gp, 12, 6)** from `ShapeFunctionOperator` registered as for `LinearEulerBernoulliBeamElement3D`.

- **Mass path only:** `extend_natural_shape_to_warping` builds an extended \(\mathbf{N}\) of shape **(14, 6)** per Gauss point: rows 0–11 copy EB; rows 12–13 place linear Lagrange \(L_1, L_2\) in the \(\theta_x\) **component column** so consistent mass can couple \(\chi\) DOFs.

- **Stiffness row 6 of B** does **not** use an extended \(\partial\mathbf{N}/\partial\xi\) for \(\chi\). Instead, `WarpingStrainDisplacementOperator` forms row 6 from the existing **EB** \(\mathbf{B}_{6\times12}\): twist part from **row 5**, \(\partial\chi/\partial x\) from the **same** linear slopes as axial strain (row 0, columns 0 and 6), mapped to columns 12–13. See §10.

**Code:** [../euler_bernoulli/utilities/shape_functions.py](../euler_bernoulli/utilities/shape_functions.py), [utilities/shape_functions.py](utilities/shape_functions.py).

---

## 6. Kinematic assumptions

All assumptions in [EB README §5](../euler_bernoulli/README.md) apply (no transverse shear, small displacements, straight element, uniform section along the element, diagonal **D** without shear-centre coupling in the implemented matrices).

**Additional / modified for warping:**

1. **Vlasov thin-wall idealisation:** cross-section is rigid in its plane; out-of-plane warping is described by a scalar intensity \(\chi\) per node in this 1D element (not full sectorial field interpolation on the 2D section).
2. **Linear \(\chi\) along the element:** \(\chi\) uses the same two-node linear Lagrange pattern as \(u_x\); \(\partial\chi/\partial x\) is **piecewise constant** on each element (adequate for slender meshes; refine near strong warping gradients).
3. **Uncoupled warping in D:** \(D_{6i} = 0\) for \(i \neq 6\) — no bimoment–bending or bimoment–shear coupling terms in this implementation.
4. **Strain rows 5 and 6:** both involve \(\partial\theta_x/\partial x\); **constitutive** split (\(G J_t\) vs \(E\Gamma\)) separates St. Venant torsion from warping stiffness rather than merging them into a single strain measure.

---

## 7. Strain vector (7-component Voigt)

\[
\boldsymbol{\varepsilon} = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x,\, \varepsilon_w]^\top \in \mathbb{R}^7 .
\]

| Row | Symbol | Definition | Role |
|-----|--------|------------|------|
| 0–4 | — | Same as EB | Axial, bending, zero shear |
| 5 | \(\phi_x\) | \(\partial\theta_x/\partial x\) | St. Venant twist rate |
| 6 | \(\varepsilon_w\) (warping strain row) | \(\phi_x' = \partial\theta_x/\partial x + \partial\chi/\partial x\) | Non-uniform torsion / bimoment-type rate |

**Note on row 5 vs 6.** Row 6 **includes** \(\partial\theta_x/\partial x\) **and** \(\partial\chi/\partial x\). Row 5 isolates \(\partial\theta_x/\partial x\) for \(G J_t\). The **seventh** stress resultant \(S_w\) is work-conjugate to row 6 via \(D_{66} = E\Gamma\).

---

## 8. Stress resultants (7-component)

\[
\mathbf{S} = [N,\, M_y,\, M_z,\, V_y,\, V_z,\, T,\, S_w]^\top .
\]

- \(S_w\) denotes the **warping** resultant paired with \(\varepsilon_w\) (bimoment-type, notation as in code/docstrings).
- \(V_y, V_z\) remain zero from \(\mathbf{D}\boldsymbol{\varepsilon}\) as in EB; shear from equilibrium if needed.

---

## 9. Material matrix **D** (7×7)

Block structure:

\[
\mathbf{D} =
\begin{bmatrix}
\mathbf{D}_{\mathrm{EB}}^{(6\times6)} & \mathbf{0} \\
\mathbf{0}^\top & E\Gamma
\end{bmatrix},
\quad
\mathbf{D}_{\mathrm{EB}} = \mathrm{diag}(EA,\, EI_y,\, EI_z,\, 0,\, 0,\, GJ_t).
\]

```text
D[0:6, 0:6]  = same as MaterialStiffnessOperator (linear EB)
D[6, 6]      = E * Gamma
D[6, i] = D[i, 6] = 0  for i = 0..5
```

**Code:** `WarpingMaterialStiffnessOperator.assembly_form` in [utilities/D_matrix.py](utilities/D_matrix.py).

---

## 10. Strain–displacement matrix **B** (7×14)

- Batch shape **(n_gp, 7, 14)**.

**Pipeline**

1. \(\mathbf{B}_{6\times12} = \texttt{StrainDisplacementOperator.physical\_coordinate\_form}(\partial\mathbf{N}/\partial\xi, \partial^2\mathbf{N}/\partial\xi^2)\).
2. \(\mathbf{B}_{:,0:6,\,0:12} = \mathbf{B}_{6\times12}\).
3. Row 6:
   - \(\mathbf{B}_{6,\,0:12} = \mathbf{B}_{6\times12}[5, :]\) (twist rate part),
   - \(\mathbf{B}_{6,\,12} = \mathbf{B}_{6\times12}[0, 0]\), \(\mathbf{B}_{6,\,13} = \mathbf{B}_{6\times12}[0, 6]\) (\(\partial\chi/\partial x\) using the same \(L_1, L_2\) slopes as \(\varepsilon_x\)).

For 2-node linear \(L_1, L_2\), row 6 is **constant** across Gauss points.

**Code:** `WarpingStrainDisplacementOperator.physical_coordinate_form` in [utilities/B_matrix.py](utilities/B_matrix.py).

---

## 11. Weak form and stiffness / force

\[
\mathbf{K}_e = \sum_g \mathbf{B}_g^\top \mathbf{D}\,\mathbf{B}_g\, w_g\, \det J,
\quad \mathbf{K}_e \in \mathbb{R}^{14\times 14}.
\]

**Force vector** \(\mathbf{F}_e \in \mathbb{R}^{14}\): distributed and point loads populate **only the first 12 entries** via the same EB shape tensor **(12, 6)** as linear EB. **No** direct external load is assembled on \(\chi\) DOFs (indices 12–13) in the current implementation; warping is driven by stiffness coupling and boundary conditions.

**Code:** `element_stiffness_matrix`, `element_force_vector` in [linear_euler_bernoulli_3D.py](../euler_bernoulli/linear_euler_bernoulli_3D.py).

---

## 12. Consistent mass matrix (14×14)

Same Gauss loop pattern as EB, with per-DOF weights extended by **\(\rho\Gamma\)** on DOFs 12 and 13. The shape-function matrix \(\mathbf{N}_g\) is **(14, 6)** from `extend_natural_shape_to_warping`. Pairwise averaging \(m_{ij} = \tfrac{1}{2}(\mu_i + \mu_j)\) as in EB.

**Code:** `element_mass_matrix` in [linear_euler_bernoulli_3D.py](../euler_bernoulli/linear_euler_bernoulli_3D.py).

---

## 13. Capabilities

- All capabilities of [linear EB](../euler_bernoulli/README.md) on the first 12 DOFs (axial, bending, St. Venant torsion, consistent mass for those DOFs).
- Additional **warping stiffness** \(E\Gamma\) and **warping mass** contribution when \(\Gamma > 0\) and density is set.
- **Graceful degradation:** if \(\Gamma = 0\), \(D_{66} = 0\); the model has no warping stiffness but remains a 14-DOF element.

---

## 14. Assumptions and limitations

| Topic | Limitation |
|-------|------------|
| EB baseline | Same as [EB §13](../euler_bernoulli/README.md) (no shear deformation, geometric linearity, straight element, …). |
| Vlasov thin-wall | Intended for **open** thin-walled sections; **closed** tubes or thick-walled distortion are not represented by this 1D warping DOF. |
| Linear \(\chi\) | Piecewise constant \(\partial\chi/\partial x\) per element; mesh refinement near warping restraints may be needed. |
| **D** structure | No coupling between warping and bending/torsion off-diagonals; no shear-centre offset in **D**. |
| Loads | No distributed **warping** load on \(\chi\) DOFs; boundary-driven warping only via prescribed DOFs or indirectly through **K**. |
| Section data | `section_array` must have length in **\{5, 7, 9, 10\}**; \(\Gamma\) uses index **9** only if length \(\geq 10\). |
| Material | `material_array` size **4** (same convention as EB for \(\rho\), \(E\), \(G\), etc.). |

---

## 15. Operator / code reference

| Concept | Class / method | File |
|---------|----------------|------|
| Shape functions **N** (12×6) | `ShapeFunctionOperator.natural_coordinate_form` | [../euler_bernoulli/utilities/shape_functions.py](../euler_bernoulli/utilities/shape_functions.py) |
| Extended **N** (14×6) for mass | `extend_natural_shape_to_warping` | [utilities/shape_functions.py](utilities/shape_functions.py) |
| Strain–displacement **B** (7×14) | `WarpingStrainDisplacementOperator.physical_coordinate_form` | [utilities/B_matrix.py](utilities/B_matrix.py) |
| Material **D** (7×7) | `WarpingMaterialStiffnessOperator.assembly_form` | [utilities/D_matrix.py](utilities/D_matrix.py) |
| \(\mathbf{K}_e\) (14×14) | `element_stiffness_matrix` | [linear_euler_bernoulli_3D.py](../euler_bernoulli/linear_euler_bernoulli_3D.py) |
| \(\mathbf{F}_e\) (14,) | `element_force_vector` | [linear_euler_bernoulli_3D.py](../euler_bernoulli/linear_euler_bernoulli_3D.py) |
| \(\mathbf{M}_e\) (14×14) | `element_mass_matrix` | [linear_euler_bernoulli_3D.py](../euler_bernoulli/linear_euler_bernoulli_3D.py) |
