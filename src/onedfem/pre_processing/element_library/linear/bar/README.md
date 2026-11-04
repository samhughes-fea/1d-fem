# Linear Bar element (3D, axial + torsion)

2-node 3D Bar element with axial and torsional deformation only. The element uses the **generalised 3D** layout: **K_e** and **F_e** are always **(12×12)** and **(12×1)** (6 DOF per node × 2 nodes); inactive DOF appear as zero entries.

---

## Reference configuration and coordinate mapping

- **Reference**: Initial geometry; element length \(L\); local axis from the grid (axial unit vector).
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

For the **Bar** element only two components are active; the rest are zero:

| Component | Definition | Bar |
|-----------|------------|-----|
| \(\varepsilon_x\) | Axial strain \(\partial u_{\mathrm{axial}}/\partial x\) | Active |
| \(\kappa_y\) | Curvature about y | 0 |
| \(\kappa_z\) | Curvature about z | 0 |
| \(\gamma_{xy}\) | Shear (xy) | 0 |
| \(\gamma_{xz}\) | Shear (xz) | 0 |
| \(\phi_x\) | Torsion \(\partial\theta_x/\partial x\) | Active |

**Reduced form used in code**: \(\varepsilon_{\mathrm{red}} = [\varepsilon_{\mathrm{axial}},\, \phi_x]^\top\), shape **(2,)**. This maps to the full 6-component view with \(\varepsilon_x = \varepsilon_{\mathrm{axial}}\), \(\phi_x\) as above, and the other four entries zero.

---

## Stress resultants (complete set)

The **complete** set of stress resultants is **N** (axial force), **M_y**, **M_z** (bending moments), **V_y**, **V_z** (shear forces), **T** (torsional moment). For the Bar element:

- **N**: from axial strain (constitutive).
- **T**: from torsion (constitutive).
- **M_y**, **M_z**, **V_y**, **V_z**: zero (no bending or shear stiffness).

Constitutive relation in reduced form: \([N_{\mathrm{axial}},\, M_{\mathrm{torsion}}]^\top = \mathbf{D}\,\varepsilon_{\mathrm{red}}\).

---

## Material matrix D

- **Reduced form (used in assembly)**: shape **(2, 2)**,
  \[
  \mathbf{D} = \mathrm{diag}(EA,\, GJ_t), \quad
  EA = E\cdot A, \quad GJ_t = G\cdot J_t.
  \]
- In the full 6-component stress-resultant view, **D** would be 6×6 with only the (1,1) and (6,6) entries non-zero (axial and torsion); the implementation uses the 2×2 reduced **D** and assembles a 12×12 **K_e** with zeros elsewhere.

**Code**: `MaterialStiffnessOperator` in [utilities/D_matrix.py](utilities/D_matrix.py).

---

## Strain–displacement matrix B

- **Shape**: at each Gauss point, **B** has shape **(2, 12)** in the reduced form (2 strain components, 12 element DOF).
- Relation: \(\varepsilon_{\mathrm{red}} = \mathbf{B}\,\mathbf{u}_e\), with \(\mathbf{u}_e\) of shape **(12,)**.
- **Rows**: row 0 → \(\varepsilon_{\mathrm{axial}}\) (axial strain from translational DOF along the bar axis); row 1 → \(\phi_x\) (torsion from \(\theta_x\) at the two nodes).
- **B** is constant in \(\xi\) (linear shape functions in \(x\)). Derivatives w.r.t. \(x\) use \(\partial\xi/\partial x = 2/L\); for the bar, axial strain is \((u_{\mathrm{axial},2} - u_{\mathrm{axial},1})/L\), torsion is \((\theta_{x,2} - \theta_{x,1})/L\).

**Code**: `StrainDisplacementOperator.natural_coordinate_form` in [utilities/B_matrix.py](utilities/B_matrix.py). Returns **B** of shape **(n_gauss, 2, 12)**.

---

## Element stiffness and force

- **Stiffness**: \(\mathbf{K}_e = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,\mathrm{d}x = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,|J|\,\mathrm{d}\xi\). **Shape**: **(12, 12)**.
- **Force**: \(\mathbf{F}_e\) of shape **(12,)** (consistent nodal forces from distributed/point loads).
- **Quadrature**: Gauss–Legendre; **one** point suffices for exact integration (constant **B**).
- Jacobian \(|J| = L/2\).

**Code**: [linear_bar_3D.py](linear_bar_3D.py) (`element_stiffness_matrix`, `element_force_vector`).

---

## Operator / code reference

| Concept | Class / method |
|--------|-----------------|
| Strain–displacement **B** | `StrainDisplacementOperator.natural_coordinate_form`, [utilities/B_matrix.py](utilities/B_matrix.py) |
| Material matrix **D** | `MaterialStiffnessOperator.assembly_form`, [utilities/D_matrix.py](utilities/D_matrix.py) |
| Element stiffness **K_e** | `LinearBarElement3D.element_stiffness_matrix`, [linear_bar_3D.py](linear_bar_3D.py) |
| Element force **F_e** | `LinearBarElement3D.element_force_vector`, [linear_bar_3D.py](linear_bar_3D.py) |
