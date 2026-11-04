# Linear Timoshenko beam element (3D)

2-node 3D Timoshenko beam element: axial, bending (rotation-based curvature), shear, and torsion. The element uses the **generalised 3D** layout: **K_e** and **F_e** are always **(12×12)** and **(12×1)** (6 DOF per node × 2 nodes).

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
- Jacobian \(|J| = L/2\). Timoshenko uses first derivatives only for curvature (rotation-based); \(\partial^2\xi/\partial x^2\) is available but not used for curvature in the strain definition.

---

## Strain vector (full 6-component view)

The **complete** strain vector is
\[
\varepsilon = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top \quad \text{shape } (6,).
\]

Definitions for **Timoshenko** (all six components active):

| Component | Definition | Timoshenko |
|-----------|------------|------------|
| \(\varepsilon_x\) | Axial strain \(\partial u_x/\partial x\) | Active |
| \(\kappa_y\) | Curvature about y: \(\partial\theta_y/\partial x\) (rotation-based) | Active |
| \(\kappa_z\) | Curvature about z: \(\partial\theta_z/\partial x\) (rotation-based) | Active |
| \(\gamma_{xy}\) | Shear (xy): \(\partial u_y/\partial x - \theta_z\) | Active |
| \(\gamma_{xz}\) | Shear (xz): \(\partial u_z/\partial x + \theta_y\) (or \(\partial u_z/\partial x - \theta_y\) per sign convention in code) | Active |
| \(\phi_x\) | Torsion \(\partial\theta_x/\partial x\) | Active |

---

## Stress resultants (complete set)

The **complete** set is **N**, **M_y**, **M_z**, **V_y**, **V_z**, **T**. Constitutive relation \(\mathbf{N} = \mathbf{D}\,\varepsilon\) gives all six: axial force, bending moments, shear forces, and torsional moment from the corresponding strain components.

---

## Material matrix D

- **Shape**: **(6, 6)**.
- **Entries** (Voigt order as above):
  \[
  \mathbf{D} = \mathrm{diag}(EA,\, EI_y,\, EI_z,\, \kappa GA,\, \kappa GA,\, GJ_t),
  \]
  with \(EA = E\cdot A\), \(EI_y = E\cdot I_y\), \(EI_z = E\cdot I_z\), \(\kappa GA = \kappa\cdot G\cdot A\) (shear correction factor \(\kappa\), default 5/6), \(GJ_t = G\cdot J_t\). All diagonal entries are non-zero.

**Code**: `MaterialStiffnessOperator` in [utilities/D_matrix.py](utilities/D_matrix.py).

---

## Strain–displacement matrix B

- **Shape**: at each Gauss point, **B** has shape **(6, 12)**.
- Relation: \(\varepsilon = \mathbf{B}\,\mathbf{u}_e\), \(\mathbf{u}_e\) shape **(12,)**.
- **Rows**: row 0 → \(\varepsilon_x\) (\(\partial u_x/\partial x\)); row 1 → \(\kappa_y\) (\(\partial\theta_y/\partial x\)); row 2 → \(\kappa_z\) (\(\partial\theta_z/\partial x\)); row 3 → \(\gamma_{xy}\) (\(\partial u_y/\partial x - \theta_z\)); row 4 → \(\gamma_{xz}\) (\(\partial u_z/\partial x - \theta_y\)); row 5 → \(\phi_x\) (\(\partial\theta_x/\partial x\)). Shear rows use both derivative and shape-function terms (no \(\mathrm{d}\xi/\mathrm{d}x\) on the rotation terms).
- Physical coordinates: \(\partial/\partial x = (2/L)\,\partial/\partial\xi\).

**Code**: `StrainDisplacementOperator.physical_coordinate_form` in [utilities/B_matrix.py](utilities/B_matrix.py). Returns **B** of shape **(n_gauss, 6, 12)**. Requires shape functions **N** for the shear terms.

---

## Element stiffness and force

- **Stiffness**: \(\mathbf{K}_e = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,\mathrm{d}x = \int \mathbf{B}^\top \mathbf{D}\,\mathbf{B}\,|J|\,\mathrm{d}\xi\). **Shape**: **(12, 12)**.
- **Force**: \(\mathbf{F}_e\) shape **(12,)**.
- **Quadrature**: Material stiffness uses ``assemble_timoshenko_K0`` with ``TimoshenkoQuadratureOrders`` from the mesh integration columns (``axial``, ``bending_y``, ``bending_z``, ``shear_y``, ``shear_z``, ``torsion``). The full-rule order is their maximum; bending and shear-stiffness rows use separate Gauss rules (default shear block = 1 point). Post-processing caches Gauss data on that full rule. Mass and distributed loads use ``loop_order`` (\(\max(\cdot), 2\)).
- Jacobian \(|J| = L/2\).

**Code**: [linear_timoshenko_3D.py](linear_timoshenko_3D.py); shared assembly [utilities/k0_timoshenko.py](utilities/k0_timoshenko.py).

---

## Operator / code reference

| Concept | Class / method |
|--------|-----------------|
| Strain–displacement **B** | `StrainDisplacementOperator.physical_coordinate_form`, [utilities/B_matrix.py](utilities/B_matrix.py) |
| Material matrix **D** | `MaterialStiffnessOperator.assembly_form`, [utilities/D_matrix.py](utilities/D_matrix.py) |
| Element stiffness **K_e** | `assemble_timoshenko_K0`, [utilities/k0_timoshenko.py](utilities/k0_timoshenko.py); `LinearTimoshenkoBeamElement3D.element_stiffness_matrix`, [linear_timoshenko_3D.py](linear_timoshenko_3D.py) |
| Element force **F_e** | `LinearTimoshenkoBeamElement3D.element_force_vector`, [linear_timoshenko_3D.py](linear_timoshenko_3D.py) |
