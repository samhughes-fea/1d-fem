# Linear Truss element (3D, axial + transverse shear + torsion)

2-node 3D Truss element with axial, transverse shear, and torsional deformation. The element uses the **generalised 3D** layout: **K_e** and **F_e** are always **(12×12)** and **(12×1)** (6 DOF per node × 2 nodes); inactive DOF appear as zero entries.

---

## Reference configuration and coordinate mapping

- **Reference**: Initial geometry; element length \(L\); local axial and transverse unit vectors from the grid.
- **Coordinate mapping**:
  \[
  x(\xi) = \frac{1-\xi}{2} x_1 + \frac{1+\xi}{2} x_2, \quad
  \frac{\mathrm{d}x}{\mathrm{d}\xi} = \frac{L}{2}, \quad
  \frac{\partial\xi}{\partial x} = \frac{2}{L}.
  \]
- Jacobian \(|J| = \mathrm{d}x/\mathrm{d}\xi = L/2\).

---

## Strain vector (full 6-component view)

The **complete** strain vector in the program’s generalised 3D view is
\[
\varepsilon = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top \quad \text{shape } (6,).
\]

For the **Truss** element three components are active; the rest are zero:

| Component | Definition | Truss |
|-----------|------------|-------|
| \(\varepsilon_x\) | Axial strain \(\partial u_{\mathrm{axial}}/\partial x\) | Active |
| \(\kappa_y\) | Curvature about y | 0 |
| \(\kappa_z\) | Curvature about z | 0 |
| \(\gamma_{\mathrm{transverse}}\) | Transverse shear (in local transverse direction) | Active |
| \(\phi_x\) | Torsion \(\partial\theta_x/\partial x\) | Active |

**Reduced form used in code**: \(\varepsilon_{\mathrm{red}} = [\varepsilon_{\mathrm{axial}},\, \gamma_{\mathrm{transverse}},\, \phi_x]^\top\), shape **(3,)**. This maps to the full 6-component view with \(\varepsilon_x\), one shear component (e.g. \(\gamma_{xy}\) or \(\gamma_{xz}\) in local convention), and \(\phi_x\) active; \(\kappa_y\), \(\kappa_z\), and the other shear are zero.

---

## Stress resultants (complete set)

The **complete** set of stress resultants is **N** (axial force), **M_y**, **M_z** (bending moments), **V_y**, **V_z** (shear forces), **T** (torsional moment). For the Truss element:

- **N**: from axial strain (constitutive).
- **V** (transverse): from transverse shear strain (constitutive; one component in local frame).
- **T**: from torsion (constitutive).
- **M_y**, **M_z**: zero (no bending stiffness).

Constitutive relation in reduced form: stress resultants from \(\mathbf{N}_{\mathrm{red}} = \mathbf{D}\,\varepsilon_{\mathrm{red}}\).

---

## Material matrix D

- **Reduced form (used in assembly)**: shape **(3, 3)**,
  \[
  \mathbf{D} = \mathrm{diag}(EA,\, \kappa GA,\, GJ_t), \quad
  EA = E\cdot A, \quad \kappa GA = \kappa\cdot G\cdot A, \quad GJ_t = G\cdot J_t.
  \]
  \(\kappa\) is the shear correction factor (default 5/6).
- In the full 6-component stress-resultant view, **D** would be 6×6 with only three diagonal entries non-zero (axial, transverse shear, torsion); the implementation uses the 3×3 reduced **D** and assembles a 12×12 **K_e**.

**Code**: `MaterialStiffnessOperator` in [utilities/D_matrix.py](utilities/D_matrix.py).

---

## Strain–displacement matrix B

- **Shape**: at each Gauss point, **B** has shape **(3, 12)** in the reduced form (3 strain components, 12 element DOF).
- Relation: \(\varepsilon_{\mathrm{red}} = \mathbf{B}\,\mathbf{u}_e\), with \(\mathbf{u}_e\) of shape **(12,)**.
- **Rows**: row 0 → \(\varepsilon_{\mathrm{axial}}\); row 1 → \(\gamma_{\mathrm{transverse}}\) (transverse displacement derivative along local transverse direction); row 2 → \(\phi_x\) (torsion).
- **B** is constant in \(\xi\). Derivatives w.r.t. \(x\) use \(\partial\xi/\partial x = 2/L\); axial and transverse strains are \((u_2 - u_1)/L\) in the respective directions, torsion is \((\theta_{x,2} - \theta_{x,1})/L\).

**Code**: `StrainDisplacementOperator.natural_coordinate_form` in [utilities/B_matrix.py](utilities/B_matrix.py). Returns **B** of shape **(n_gauss, 3, 12)**.

---

## Element stiffness and force

- **Stiffness**: \(\mathbf{K}_e = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,\mathrm{d}x = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,|J|\,\mathrm{d}\xi\). **Shape**: **(12, 12)**.
- **Force**: \(\mathbf{F}_e\) of shape **(12,)**.
- **Quadrature**: Gauss–Legendre; **one** point suffices (constant **B**).
- Jacobian \(|J| = L/2\).

**Code**: [linear_truss_3D.py](linear_truss_3D.py) (`element_stiffness_matrix`, `element_force_vector`).

---

## Operator / code reference

| Concept | Class / method |
|--------|-----------------|
| Strain–displacement **B** | `StrainDisplacementOperator.natural_coordinate_form`, [utilities/B_matrix.py](utilities/B_matrix.py) |
| Material matrix **D** | `MaterialStiffnessOperator.assembly_form`, [utilities/D_matrix.py](utilities/D_matrix.py) |
| Element stiffness **K_e** | `LinearTrussElement3D.element_stiffness_matrix`, [linear_truss_3D.py](linear_truss_3D.py) |
| Element force **F_e** | `LinearTrussElement3D.element_force_vector`, [linear_truss_3D.py](linear_truss_3D.py) |
