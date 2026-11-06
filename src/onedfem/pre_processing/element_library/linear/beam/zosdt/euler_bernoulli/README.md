# Linear Euler–Bernoulli beam element (3D)

This document is a formulation reference for the 2-node 3D **Euler–Bernoulli** beam in **zero-order shear deformation theory** (ZOSDT): transverse shear strains are set to zero kinematically. The implementation is [`LinearEulerBernoulliBeamElement3D`](linear_euler_bernoulli_3D.py). Voigt ordering follows [`docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`](../../../../../docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md).

**Relationship to other elements.** The 12-DOF displacement block and the first six strain rows of **B** / **D** are the baseline embedded unchanged in the warping extension: see [Euler–Bernoulli with Vlasov warping](../euler_bernoulli_with_warp/README.md).

---

## 1. Element identity

| Item | Value |
|------|--------|
| Class | `LinearEulerBernoulliBeamElement3D` |
| Nodes | 2 |
| DOFs | 6 per node → **12** element DOFs |
| Kinematics | Euler–Bernoulli (no shear deformation) |
| Stiffness | Geometrically linear |

---

## 2. Reference configuration and coordinate mapping

- **Geometry.** Chord length \(L\); local **x** runs along the undeformed straight element from node 1 to node 2.
- **Isoparametric map** \(\xi \in [-1,1]\):
  \[
  x(\xi) = \frac{1-\xi}{2} x_1 + \frac{1+\xi}{2} x_2, \quad
  \frac{\mathrm{d}x}{\mathrm{d}\xi} = \frac{L}{2}, \quad
  \frac{\partial\xi}{\partial x} = \frac{2}{L}, \quad
  \frac{\partial^2\xi}{\partial x^2} = \frac{4}{L^2}.
  \]
- **Jacobian** (1D): \(\det J = L/2\) (used as `detJ` in code). Curvature uses \(\partial^2/\partial x^2 = (4/L^2)\,\partial^2/\partial\xi^2\).

---

## 3. Degrees of freedom and displacement vector

**Node-major** element vector \(\mathbf{U}_e \in \mathbb{R}^{12}\):

\[
\mathbf{U}_e = [u_x^1,\, u_y^1,\, u_z^1,\, \theta_x^1,\, \theta_y^1,\, \theta_z^1,\,
u_x^2,\, u_y^2,\, u_z^2,\, \theta_x^2,\, \theta_y^2,\, \theta_z^2]^\top .
\]

| Index | DOF | Physical meaning |
|-------|-----|------------------|
| 0–2 | \(u_x^1, u_y^1, u_z^1\) | Translations at node 1 |
| 3–5 | \(\theta_x^1, \theta_y^1, \theta_z^1\) | Rotations at node 1 |
| 6–8 | \(u_x^2, u_y^2, u_z^2\) | Translations at node 2 |
| 9–11 | \(\theta_x^2, \theta_y^2, \theta_z^2\) | Rotations at node 2 |

---

## 4. Shape functions

- **Axial** (\(u_x\)) and **torsion** (\(\theta_x\)) use **linear Lagrange** polynomials \(L_1(\xi), L_2(\xi)\) on two nodes.
- **Bending** in each plane uses **Hermite cubic** shape functions (standard \(C^1\) beam pair) so that transverse displacement and section rotation are coupled along the element.

The registry operator returns, for Gauss points \(\xi_g\), tensors of shape **(n_gp, 12, 6)**:

- \(\mathbf{N}\): shape functions  
- \(\partial\mathbf{N}/\partial\xi\), \(\partial^2\mathbf{N}/\partial\xi^2\)

Rows index **global DOF** (0…11); columns index **displacement components** \((u_x, u_y, u_z, \theta_x, \theta_y, \theta_z)\). At Gauss point \(g\), \(\mathbf{u}_g = \mathbf{N}_g \mathbf{U}_e\) with \(\mathbf{u}_g \in \mathbb{R}^6\).

**Code:** `ShapeFunctionOperator.natural_coordinate_form` in [utilities/shape_functions.py](utilities/shape_functions.py).

---

## 5. Kinematic assumptions (Euler–Bernoulli)

These are the assumptions that define the discrete operators:

1. **Plane sections remain plane** and **normal** to the deformed mid-axis (Bernoulli hypothesis).
2. **No transverse shear deformation:** \(\gamma_{xy} = \gamma_{xz} = 0\) **kinematically** (not merely “small”).
3. **Kinematic rotations:** \(\theta_y = -\partial u_z/\partial x\), \(\theta_z = \partial u_y/\partial x\) (sign convention as in implementation).
4. **Small displacements and rotations** (linearized geometry).
5. **Linear elastic**, homogeneous, **isotropic** material in the 1D sense (axial/bending use \(E\); St. Venant torsion uses \(G\) and \(J_t\)).
6. **Uniform** cross-section properties along the element.
7. **Uncoupled bending planes** in **D**: no shear-centre or centroid offset coupling in the implemented diagonal **D**.
8. **St. Venant torsion only** in this element: **no** warping DOF, **no** bimoment. For open thin-walled sections with non-uniform torsion, use [warping EB](../euler_bernoulli_with_warp/README.md).

---

## 6. Strain vector (6-component Voigt)

\[
\boldsymbol{\varepsilon} = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top \in \mathbb{R}^6 .
\]

| Row | Component | Definition | EB |
|-----|-----------|------------|-----|
| 0 | \(\varepsilon_x\) | \(\partial u_x/\partial x\) | Active |
| 1 | \(\kappa_y\) | Curvature about \(y\) (from \(u_z\), \(\theta_y\)) | Active |
| 2 | \(\kappa_z\) | Curvature about \(z\) (from \(u_y\), \(\theta_z\)) | Active |
| 3 | \(\gamma_{xy}\) | Shear strain | **0** |
| 4 | \(\gamma_{xz}\) | Shear strain | **0** |
| 5 | \(\phi_x\) | \(\partial \theta_x/\partial x\) (twist rate) | Active |

**Shear forces** \(V_y, V_z\) are **not** obtained from \(\mathbf{D}\boldsymbol{\varepsilon}\) (rows 3–4 of **D** are zero by construction). Where needed, they follow from **equilibrium**, e.g. \(V = \mathrm{d}M/\mathrm{d}x\).

---

## 7. Stress resultants

Work-conjugate to \(\boldsymbol{\varepsilon}\) in Voigt order:

\[
\mathbf{S} = [N,\, M_y,\, M_z,\, V_y,\, V_z,\, T]^\top .
\]

- \(N, M_y, M_z, T\) follow from \(\mathbf{D}\boldsymbol{\varepsilon}\) on the active diagonal entries.
- \(V_y, V_z\) from \(\mathbf{D}\boldsymbol{\varepsilon}\) are **zero**; recover from equilibrium if required for post-processing.

---

## 8. Material matrix **D**

- **Shape:** \((6, 6)\), symmetric diagonal in this implementation.
- Voigt order as in §6:

\[
\mathbf{D} = \mathrm{diag}(EA,\, EI_y,\, EI_z,\, 0,\, 0,\, GJ_t),
\]

with \(EA = E A\), \(EI_y = E I_y\), \(EI_z = E I_z\), \(GJ_t = G J_t\).

**Rows 3 and 4 (shear)** are **exactly zero** because \(\gamma_{xy} = \gamma_{xz} = 0\) is a **kinematic** constraint of the theory, not a numerical approximation.

**Code:** `MaterialStiffnessOperator.assembly_form` in [utilities/D_matrix.py](utilities/D_matrix.py).

---

## 9. Strain–displacement matrix **B**

- **Per Gauss point:** \(\boldsymbol{\varepsilon} = \mathbf{B}\,\mathbf{U}_e\) with \(\mathbf{B} \in \mathbb{R}^{6\times 12}\); batch shape **(n_gp, 6, 12)**.

**Assembly pipeline**

1. `ShapeFunctionOperator.natural_coordinate_form` → \(\partial\mathbf{N}/\partial\xi\), \(\partial^2\mathbf{N}/\partial\xi^2\).
2. `StrainDisplacementOperator.physical_coordinate_form` applies \(\partial\xi/\partial x = 2/L\) and \(\partial^2\xi/\partial x^2 = 4/L^2\) to build **physical** \(\mathbf{B}\).

**Row content (summary)**

| Strain row | Primary fields | Derivative order in \(x\) |
|------------|----------------|---------------------------|
| 0 \(\varepsilon_x\) | \(u_x\) | First (\(\sim \partial N/\partial x\)) |
| 1–2 \(\kappa_y, \kappa_z\) | \(u_z, u_y, \theta_y, \theta_z\) | Second (\(\sim \partial^2 N/\partial x^2\)) |
| 3–4 \(\gamma_{xy}, \gamma_{xz}\) | — | Rows left **zero** |
| 5 \(\phi_x\) | \(\theta_x\) | First |

A **natural-coordinate** form \(\tilde{\mathbf{B}}\) (without \(2/L\) and \(4/L^2\)) exists in `natural_coordinate_form`; the element uses **physical** \(\mathbf{B}\) for \(\mathbf{K}_e\).

**Code:** [utilities/B_matrix.py](utilities/B_matrix.py).

---

## 10. Weak form and element stiffness / force

**Stiffness (Gauss–Legendre):**

\[
\mathbf{K}_e = \sum_g \mathbf{B}_g^\top \mathbf{D}\,\mathbf{B}_g\, w_g\, \det J ,
\quad \mathbf{K}_e \in \mathbb{R}^{12\times 12}.
\]

**Distributed loads** (conceptually): contributions of the form \(\int \mathbf{N}^\top \mathbf{q}\,\mathrm{d}x\) with consistent quadrature; **point loads** use \(\mathbf{N}\) evaluated at the load station.

**Code:** `element_stiffness_matrix`, `element_force_vector` in [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py). Distributed loads use `LoadInterpolationOperator` ([utilities/interpolate_loads.py](utilities/interpolate_loads.py)).

---

## 11. Consistent mass matrix

**Shape:** \(\mathbf{M}_e \in \mathbb{R}^{12\times 12}\), symmetric.

Lumped **per-DOF inertia weights** \(\mu_i\) (density \(\rho\) from material array):

| DOF indices | Weight |
|-------------|--------|
| 0,1,2,6,7,8 (translations) | \(\rho A\) |
| 3,9 (\(\theta_x\)) | \(\rho J_t\) |
| 4,10 (\(\theta_y\)) | \(\rho I_y\) |
| 5,11 (\(\theta_z\)) | \(\rho I_z\) |

Assembly uses shape functions \(\mathbf{N}_g\) and a **pairwise average** \(m_{ij} = \tfrac{1}{2}(\mu_i + \mu_j)\) for the coupling between DOFs \(i\) and \(j\), integrated with Gauss weights and \(\det J\).

**Code:** `element_mass_matrix` → `MassObject` in [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py).

---

## 12. Capabilities

- 3D **axial**, **biaxial bending**, **St. Venant torsion**, **consistent mass**.
- **Gauss–Legendre** quadrature; order configurable (element constructor / `element_array` convention).
- **Distributed** and **point** loads on the standard 12 DOFs.

---

## 13. Assumptions and limitations

| Topic | Consequence |
|-------|-------------|
| No shear strain \(\gamma_{xy}, \gamma_{xz}\) | Beam is **kinematically** stiffer than Timoshenko; for **short/deep** beams (low span-to-depth), a shear-deformable theory is often more accurate. |
| St. Venant torsion only | **No** warping restraint / **no** bimoment; open thin-walled sections may need [warping EB](../euler_bernoulli_with_warp/README.md). |
| Straight reference geometry | **Curved** geometry is not modelled by this element; see the curved-beam element family if applicable. |
| Geometric linearity | **No** large deflection or buckling in **K**; no follower loads in the linearization. |
| Diagonal **D**, no shear-centre coupling | No coupled bending–torsion from offset; specialised sections may need richer section models. |
| Uniform section along the element | **Taper** requires mesh refinement or a different formulation. |
| Data contracts | `material_array` size **4**; `section_array` size **5** or **7** (validation in element). |

---

## 14. Operator / code reference

| Concept | Class / method | File |
|---------|----------------|------|
| Shape functions \(\mathbf{N}\) | `ShapeFunctionOperator.natural_coordinate_form` | [utilities/shape_functions.py](utilities/shape_functions.py) |
| Strain–displacement \(\mathbf{B}\) | `StrainDisplacementOperator.physical_coordinate_form` | [utilities/B_matrix.py](utilities/B_matrix.py) |
| Material \(\mathbf{D}\) | `MaterialStiffnessOperator.assembly_form` | [utilities/D_matrix.py](utilities/D_matrix.py) |
| Load interpolation | `LoadInterpolationOperator.interpolate` | [utilities/interpolate_loads.py](utilities/interpolate_loads.py) |
| \(\mathbf{K}_e\) (12×12) | `LinearEulerBernoulliBeamElement3D.element_stiffness_matrix` | [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py) |
| \(\mathbf{F}_e\) (12,) | `element_force_vector` | [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py) |
| \(\mathbf{M}_e\) (12×12) | `element_mass_matrix` | [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py) |
