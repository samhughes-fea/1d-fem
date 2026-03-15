# Nonlinear Euler–Bernoulli beam element (3D, Total Lagrangian)

2-node 3D **geometrically nonlinear** Euler–Bernoulli beam: Green–Lagrange strain, no shear deformation (\(\gamma_{xy} = \gamma_{xz} = 0\)). The element uses the **generalised 3D** layout: tangent stiffness **K_T** and internal force **F_int** are **(12×12)** and **(12×1)**. Formulation is **Total Lagrangian** (all quantities referred to the initial configuration).

---

## Full formulation document

The tensor mathematics (strain, stress, tangent stiffness, Newton–Raphson) is described in:

**[docs/element_library/total_lagrangian_beam_formulation.md](../../../../docs/element_library/total_lagrangian_beam_formulation.md)**

The following summarises the **tensor overview** for the Euler–Bernoulli nonlinear element.

---

## Reference configuration and coordinate mapping

- **Reference**: Initial (undeformed) geometry; element length \(L\); same local frame as the linear formulation.
- **Coordinate mapping**: \(x(\xi) = \frac{1-\xi}{2} x_1 + \frac{1+\xi}{2} x_2\), \(\frac{\mathrm{d}x}{\mathrm{d}\xi} = L/2\), \(\frac{\partial\xi}{\partial x} = 2/L\).

---

## Strain vector (full 6-component view, Green–Lagrange)

\[
\mathbf{E} = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top \quad \text{shape } (6,).
\]

**Decomposition**: \(\mathbf{E} = \mathbf{E}_{\mathrm{lin}}(\mathbf{u}) + \mathbf{E}_{\mathrm{nl}}(\mathbf{u})\).

For **nonlinear Euler–Bernoulli**:

- **Axial (Green–Lagrange)**: \(\varepsilon_x = \frac{\partial u_x}{\partial x} + \frac{1}{2}\bigl(\frac{\partial u_x}{\partial x}\bigr)^2 + \frac{1}{2}\bigl(\frac{\partial u_y}{\partial x}\bigr)^2 + \frac{1}{2}\bigl(\frac{\partial u_z}{\partial x}\bigr)^2\).
- **Bending** \(\kappa_y\), \(\kappa_z\): linear part as in the linear beam; nonlinear corrections (axial–curvature coupling) as in the strain operator.
- **Shear**: \(\gamma_{xy} = \gamma_{xz} = 0\) (Euler–Bernoulli).
- **Torsion**: \(\phi_x = \frac{\partial\theta_x}{\partial x}\) (linear).

---

## Stress resultants (complete set)

Full set **N**, **M_y**, **M_z**, **V_y**, **V_z**, **T**. Second Piola–Kirchhoff stress resultants from \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) (same **D** as linear). Section forces from strain: **N**, **M_y**, **M_z** (and **V_y** = **V_z** = 0 from constitutive) drive the geometric stiffness **K_σ**.

---

## Tangent stiffness

\[
\mathbf{K}_T = \mathbf{K}_0 + \mathbf{K}_\sigma \quad \text{both shape } (12, 12).
\]

- **K_0 (material stiffness)**: From linearization of \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) w.r.t. \(\mathbf{u}\) using the **linear** strain–displacement operator (same as the linear Euler–Bernoulli element stiffness).
- **K_σ (geometric stiffness)**: Depends on current section forces **N**, **M_y**, **M_z** and shape-function derivatives; assembled at Gauss points.

---

## Operator / code reference

| Concept | Class / method |
|--------|-----------------|
| Green–Lagrange strain | `GreenLagrangeStrainOperator`, [utilities/green_lagrange_strain.py](utilities/green_lagrange_strain.py) |
| Linear / nonlinear strain | `strain_linear_part`, `strain_nonlinear_part` |
| Strain–displacement (linearized) | `linearized_strain_displacement` |
| Section forces from strain | `StressResultantOperator.section_forces_from_strain`, [utilities/stress_resultant.py](utilities/stress_resultant.py) |
| Geometric stiffness **K_σ** | `GeometricStiffnessOperator.assemble_K_sigma`, [utilities/geometric_stiffness.py](utilities/geometric_stiffness.py) |

For residual, Newton–Raphson, and logging details, see the [Total Lagrangian document](../../../../docs/element_library/total_lagrangian_beam_formulation.md).
