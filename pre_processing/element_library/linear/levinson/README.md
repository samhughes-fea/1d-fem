# Linear Levinson beam element (3D)

2-node 3D Levinson beam element: axial, bending (rotation-based curvature), **higher-order shear** (with \(\alpha\) terms), and torsion. The element uses the **generalised 3D** layout: **K_e** and **F_e** are always **(12×12)** and **(12×1)** (6 DOF per node × 2 nodes). Unlike Timoshenko, Levinson uses no shear correction factor \(\kappa\) in the constitutive law (shear stiffness is \(GA\)); shear deformation is refined by higher-order terms in the strain–displacement relation.

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
- Jacobian \(|J| = L/2\). The higher-order shear terms use \(\partial^2\theta/\partial x^2\) with \(\partial^2\xi/\partial x^2 = 4/L^2\).

---

## Strain vector (full 6-component view)

The **complete** strain vector is
\[
\varepsilon = [\varepsilon_x,\, \kappa_z,\, \kappa_y,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top \quad \text{shape } (6,).
\]
(Note: in this implementation the curvature components are ordered \(\kappa_z\), \(\kappa_y\) in the strain vector; the full stress-resultant set remains N, M_y, M_z, V_y, V_z, T.)

Definitions for **Levinson** (all six components active, with higher-order shear):

| Component | Definition | Levinson |
|-----------|------------|----------|
| \(\varepsilon_x\) | Axial strain \(\partial u_x/\partial x\) | Active |
| \(\kappa_z\) | Curvature (x–y plane, bending about z): \(\partial\theta_z/\partial x\) | Active |
| \(\kappa_y\) | Curvature (x–z plane, bending about y): \(\partial\theta_y/\partial x\) | Active |
| \(\gamma_{xy}\) | Shear (xy): \(\partial u_y/\partial x - \theta_z + \alpha\,\partial^2\theta_z/\partial x^2\) | Active |
| \(\gamma_{xz}\) | Shear (xz): \(\partial u_z/\partial x - \theta_y + \alpha\,\partial^2\theta_y/\partial x^2\) | Active |
| \(\phi_x\) | Torsion \(\partial\theta_x/\partial x\) | Active |

\(\alpha\) is the higher-order shear coefficient (e.g. \(h^2/12\) for rectangular sections; \(h\) = depth). It may be set from section properties or default 0.

---

## Stress resultants (complete set)

The **complete** set is **N**, **M_y**, **M_z**, **V_y**, **V_z**, **T**. Constitutive relation \(\mathbf{N} = \mathbf{D}\,\varepsilon\) gives all six. Shear stiffness is **GA** (no \(\kappa\)): \(\mathbf{D}\) uses \(GA\) for the shear diagonal entries, unlike Timoshenko’s \(\kappa GA\).

---

## Material matrix D

- **Shape**: **(6, 6)**.
- **Entries** (same Voigt order as strain):
  \[
  \mathbf{D} = \mathrm{diag}(EA,\, EI_z,\, EI_y,\, GA,\, GA,\, GJ_t),
  \]
  with \(EA = E\cdot A\), \(EI_y = E\cdot I_y\), \(EI_z = E\cdot I_z\), \(GA = G\cdot A\) (no shear correction factor), \(GJ_t = G\cdot J_t\). All diagonal entries are non-zero.

**Code**: `MaterialStiffnessOperator` in [utilities/D_matrix.py](utilities/D_matrix.py).

---

## Strain–displacement matrix B

- **Shape**: at each Gauss point, **B** has shape **(6, 12)**.
- Relation: \(\varepsilon = \mathbf{B}\,\mathbf{u}_e\), \(\mathbf{u}_e\) shape **(12,)**.
- **Rows**: row 0 → \(\varepsilon_x\); row 1 → \(\kappa_z\) (\(\partial\theta_z/\partial x\), x–y plane); row 2 → \(\kappa_y\) (\(\partial\theta_y/\partial x\), x–z plane); row 3 → \(\gamma_{xy}\) (derivative + rotation + \(\alpha\,\partial^2\theta_z/\partial x^2\)); row 4 → \(\gamma_{xz}\) (derivative + rotation + \(\alpha\,\partial^2\theta_y/\partial x^2\)); row 5 → \(\phi_x\). The \(\partial^2/\partial x^2\) terms use \(\partial^2\xi/\partial x^2 = 4/L^2\).
- Physical coordinates: \(\partial/\partial x = (2/L)\,\partial/\partial\xi\); \(\partial^2/\partial x^2 = (4/L^2)\,\partial^2/\partial\xi^2\).

**Code**: `StrainDisplacementOperator.physical_coordinate_form` in [utilities/B_matrix.py](utilities/B_matrix.py). Returns **B** of shape **(n_gauss, 6, 12)**. Requires **N** and **d2N_dξ2** for the higher-order shear terms when \(\alpha \neq 0\).

---

## Element stiffness and force

- **Stiffness**: \(\mathbf{K}_e = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,\mathrm{d}x = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,|J|\,\mathrm{d}\xi\). **Shape**: **(12, 12)**.
- **Force**: \(\mathbf{F}_e\) shape **(12,)**.
- **Quadrature**: Gauss–Legendre; number of points set in element.
- Jacobian \(|J| = L/2\).

**Code**: [linear_levinson_3D.py](linear_levinson_3D.py) (`element_stiffness_matrix`, `element_force_vector`).

---

## Operator / code reference

| Concept | Class / method |
|--------|-----------------|
| Strain–displacement **B** | `StrainDisplacementOperator.physical_coordinate_form`, [utilities/B_matrix.py](utilities/B_matrix.py) |
| Material matrix **D** | `MaterialStiffnessOperator.assembly_form`, [utilities/D_matrix.py](utilities/D_matrix.py) |
| Element stiffness **K_e** | `LinearLevinsonBeamElement3D.element_stiffness_matrix`, [linear_levinson_3D.py](linear_levinson_3D.py) |
| Element force **F_e** | `LinearLevinsonBeamElement3D.element_force_vector`, [linear_levinson_3D.py](linear_levinson_3D.py) |
