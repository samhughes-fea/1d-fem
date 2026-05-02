# Nonlinear Timoshenko beam element (3D, Total Lagrangian)

2-node 3D **geometrically nonlinear** Timoshenko beam: Green–Lagrange strain with shear deformation. The element uses the **generalised 3D** layout: tangent stiffness **K_T** and internal force **F_int** are **(12×12)** and **(12×1)** when no warping DOFs are used; with **`[warping]`** and a mesh that allocates χ per node, **(14×14)** and **(14×1)** (same Vlasov row as linear Timoshenko + warping). Formulation is **Total Lagrangian** (all quantities referred to the initial configuration).

---

## Full formulation document

The tensor mathematics (strain, stress, tangent stiffness, Newton–Raphson) is described in:

**[docs/element_library/total_lagrangian_beam_formulation.md](../../../../docs/element_library/total_lagrangian_beam_formulation.md)**

The following summarises the **tensor overview** for the Timoshenko nonlinear element.

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

For **nonlinear Timoshenko**:

- **Axial (Green–Lagrange)**: same as in the Total Lagrangian document (nonlinear in \(\partial u/\partial x\)).
- **Bending** \(\kappa_y\), \(\kappa_z\): linear part as in the linear Timoshenko beam; nonlinear corrections as in the strain operator.
- **Shear** \(\gamma_{xy}\), \(\gamma_{xz}\): same linear relation as in the linear Timoshenko formulation (optional nonlinear terms can be added in the operator).
- **Torsion**: \(\phi_x = \frac{\partial\theta_x}{\partial x}\) (linear).

---

## Stress resultants (complete set)

Full set **N**, **M_y**, **M_z**, **V_y**, **V_z**, **T**. Second Piola–Kirchhoff stress resultants from \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) (same **D** as linear Timoshenko). Section forces from strain drive the geometric stiffness **K_σ**; all six resultants are in the complete data-structure view.

---

## Tangent stiffness

\[
\mathbf{K}_T = \mathbf{K}_0 + \mathbf{K}_\sigma \quad \text{both shape } (12, 12).
\]

- **K_0 (material stiffness)**: Same ``assemble_timoshenko_K0`` / ``TimoshenkoQuadratureOrders`` as linear Timoshenko ``K_e`` (block-wise axial / bending / shear / torsion quadrature).
- **K_σ (geometric stiffness)**: Depends on current section forces **N**, **M_y**, **M_z** (and **V_y**, **V_z**, **T** in the full view) and shape-function derivatives; assembled at Gauss points using one rule of order **loop_order** (default: maximum mesh integration order, at least 2). Internal force **F_int** uses the same loop points.

**Migration:** meshes that used ``GEBTShearBeamElement3D`` should use ``NonlinearTimoshenkoBeamElement3D`` with the same ``integration_orders``; see [docs/element_library/gebt_shear_formulation.md](../../../../docs/element_library/gebt_shear_formulation.md).

---

## Operator / code reference

| Concept | Class / method |
|--------|-----------------|
| Green–Lagrange strain | `GreenLagrangeStrainOperator`, [utilities/green_lagrange_strain.py](utilities/green_lagrange_strain.py) |
| Linear / nonlinear strain | `strain_linear_part`, `strain_nonlinear_part` |
| Strain–displacement (linearized) | `linearized_strain_displacement` |
| Section forces from strain | `StressResultantOperator.section_forces_from_strain`, [utilities/stress_resultant.py](utilities/stress_resultant.py) |
| Geometric stiffness **K_σ** | Re-exported from EB utility — same `assemble_K_sigma` Gauss API, [utilities/geometric_stiffness.py](utilities/geometric_stiffness.py) |

For residual, Newton–Raphson, and logging details, see the [Total Lagrangian document](../../../../docs/element_library/total_lagrangian_beam_formulation.md).
