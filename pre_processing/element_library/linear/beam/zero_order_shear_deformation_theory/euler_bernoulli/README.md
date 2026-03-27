# Linear Euler–Bernoulli beam element (3D)

2-node 3D Euler–Bernoulli beam element: axial, bending (about y and z), and torsion; no shear deformation. The element uses the **generalised 3D** layout: **K_e** and **F_e** are always **(12×12)** and **(12×1)** (6 DOF per node × 2 nodes).

---

## Reference configuration and coordinate mapping

- **Reference**: Initial geometry; element length \(L\); local x along the element.
- **Coordinate mapping**:
  \[
  x(\xi) = \frac{1-\xi}{2} x_1 + \frac{1+\xi}{2} x_2, \quad
  \frac{\mathrm{d}x}{\mathrm{d}\xi} = \frac{L}{2}, \quad
  \frac{\partial\xi}{\partial x} = \frac{2}{L}, \quad
  \frac{\partial^2\xi}{\partial x^2} = \frac{4}{L^2}.
  \]
- Jacobian \(|J| = L/2\). Curvature terms use \(\partial^2/\partial x^2\) with \(\partial^2\xi/\partial x^2 = 4/L^2\).

---

## Strain vector (full 6-component view)

The **complete** strain vector is
\[
\varepsilon = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top \quad \text{shape } (6,).
\]

Definitions and Euler–Bernoulli usage:

| Component | Definition | Euler–Bernoulli |
|-----------|------------|-----------------|
| \(\varepsilon_x\) | Axial strain \(\partial u_x/\partial x\) | Active |
| \(\kappa_y\) | Curvature about y: \(\partial^2 w/\partial x^2\) (\(w = u_z\)) | Active |
| \(\kappa_z\) | Curvature about z: \(\partial^2 v/\partial x^2\) (\(v = u_y\)) | Active |
| \(\gamma_{xy}\) | Shear strain (xy) | **0** (no shear deformation) |
| \(\gamma_{xz}\) | Shear strain (xz) | **0** |
| \(\phi_x\) | Torsion \(\partial\theta_x/\partial x\) | Active |

Shear force in Euler–Bernoulli is not from the constitutive law (\(V_y = V_z = 0\) from \(\mathbf{D}\,\varepsilon\)); it is obtained from equilibrium (\(V = \mathrm{d}M/\mathrm{d}x\)) when needed.

---

## Stress resultants (complete set)

The **complete** set is **N** (axial), **M_y**, **M_z** (bending), **V_y**, **V_z** (shear), **T** (torsion). Constitutive relation \(\mathbf{N} = \mathbf{D}\,\varepsilon\) gives:

- **N**, **M_y**, **M_z**, **T**: from \(\mathbf{D}\,\varepsilon\) (axial, bending, torsion).
- **V_y**, **V_z**: **zero** from \(\mathbf{D}\,\varepsilon\) (shear rows of **D** are zero). Obtain shear from equilibrium if required.

---

## Material matrix D

- **Shape**: **(6, 6)**.
- **Entries** (Voigt order \(\varepsilon = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top\)):
  \[
  \mathbf{D} = \mathrm{diag}(EA,\, EI_y,\, EI_z,\, 0,\, 0,\, GJ_t),
  \]
  with \(EA = E\cdot A\), \(EI_y = E\cdot I_y\), \(EI_z = E\cdot I_z\), \(GJ_t = G\cdot J_t\). Rows 4 and 5 (shear) are zero.

**Code**: `MaterialStiffnessOperator` in [utilities/D_matrix.py](utilities/D_matrix.py).

---

## Strain–displacement matrix B

- **Shape**: at each Gauss point, **B** has shape **(6, 12)**.
- Relation: \(\varepsilon = \mathbf{B}\,\mathbf{u}_e\), \(\mathbf{u}_e\) shape **(12,)**.
- **Rows**: row 0 → \(\varepsilon_x\) (from \(\partial N/\partial x\) for \(u_x\)); rows 1–2 → \(\kappa_y\), \(\kappa_z\) (from \(\partial^2 N/\partial x^2\) for \(u_z\), \(u_y\) and \(\theta_y\), \(\theta_z\)); rows 3–4 → \(\gamma_{xy}\), \(\gamma_{xz}\) (zero, not populated); row 5 → \(\phi_x\) (from \(\partial N/\partial x\) for \(\theta_x\)).
- Physical coordinates: \(\partial/\partial x = (\partial\xi/\partial x)\,\partial/\partial\xi = (2/L)\,\partial/\partial\xi\); \(\partial^2/\partial x^2 = (4/L^2)\,\partial^2/\partial\xi^2\).

**Code**: `StrainDisplacementOperator.physical_coordinate_form` in [utilities/B_matrix.py](utilities/B_matrix.py). Returns **B** of shape **(n_gauss, 6, 12)**.

---

## Element stiffness and force

- **Stiffness**: \(\mathbf{K}_e = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,\mathrm{d}x = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,|J|\,\mathrm{d}\xi\). **Shape**: **(12, 12)**.
- **Force**: \(\mathbf{F}_e\) shape **(12,)** (consistent nodal forces from distributed and point loads).
- **Quadrature**: Gauss–Legendre; typically 2 or 3 points for stiffness (order set in element).
- Jacobian \(|J| = L/2\).

**Code**: [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py) (`element_stiffness_matrix`, `element_force_vector`).

---

## Operator / code reference

| Concept | Class / method |
|--------|-----------------|
| Strain–displacement **B** | `StrainDisplacementOperator.physical_coordinate_form`, [utilities/B_matrix.py](utilities/B_matrix.py) |
| Material matrix **D** | `MaterialStiffnessOperator.assembly_form`, [utilities/D_matrix.py](utilities/D_matrix.py) |
| Element stiffness **K_e** | `LinearEulerBernoulliBeamElement3D.element_stiffness_matrix`, [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py) |
| Element force **F_e** | `LinearEulerBernoulliBeamElement3D.element_force_vector`, [linear_euler_bernoulli_3D.py](linear_euler_bernoulli_3D.py) |
