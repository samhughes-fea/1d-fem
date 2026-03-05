# Total Lagrangian beam formulation (geometric nonlinearity)

This document describes the geometrically nonlinear 2-node 3D beam formulation used for **Euler-Bernoulli-3D-Nonlinear** and **Timoshenko-3D-Nonlinear** elements. The formulation is **Total Lagrangian**: all quantities are referred to the **initial (undeformed) configuration**.

## Reference configuration

- **Reference**: The **initial** geometry (element length \(L\), local axis from the grid) is fixed. No moving frame; the same local frame as in the linear formulation is used.
- **Coordinate mapping**: \(x(\xi) = \frac{1-\xi}{2} x_1 + \frac{1+\xi}{2} x_2\), \(\frac{\mathrm{d}x}{\mathrm{d}\xi} = L/2\), \(\frac{\partial\xi}{\partial x} = 2/L\).

## Strain: Green–Lagrange

Strain vector \(\mathbf{E} = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top\) with **nonlinear** terms in displacement:

- **Axial (Green–Lagrange)**  
  \(\varepsilon_x = \frac{\partial u_x}{\partial x} + \frac{1}{2}\left(\frac{\partial u_x}{\partial x}\right)^2 + \frac{1}{2}\left(\frac{\partial u_y}{\partial x}\right)^2 + \frac{1}{2}\left(\frac{\partial u_z}{\partial x}\right)^2\).

- **Bending** \(\kappa_y\), \(\kappa_z\): linear part as in the linear beam; nonlinear corrections as implemented in the strain operator.

- **Shear** \(\gamma_{xy}\), \(\gamma_{xz}\): zero for Euler–Bernoulli; for Timoshenko, same linear relation as in the linear formulation (optional nonlinear terms can be added).

- **Torsion** \(\phi_x = \frac{\partial\theta_x}{\partial x}\) (linear).

So \(\mathbf{E} = \mathbf{E}_{\mathrm{lin}}(\mathbf{u}) + \mathbf{E}_{\mathrm{nl}}(\mathbf{u})\).

## Stress and section forces

- **Stress**: Second Piola–Kirchhoff stress in the reference configuration; constitutive \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) (same \(\mathbf{D}\) as in linear theory).
- **Section forces**: From \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) at each Gauss point we obtain **axial force \(N\)** and **moments \(M_y\), \(M_z\)** (and torsion/shear as in the linear 6-component beam). These are the stress resultants that drive the geometric stiffness.

**Operator**: `StressResultantOperator.section_forces_from_strain(E, D)` → \((N,\, M_y,\, M_z)\).

## Internal force and tangent stiffness

- **Internal force**: From the weak form, \(\mathbf{F}_{\mathrm{int}} = \int \mathbf{B}(\mathbf{u})^\top \mathbf{S}\,\mathrm{d}V\) (beam analogue with section forces). \(\mathbf{B}\) includes the linearized strain–displacement contribution (and nonlinear contribution in the residual).
- **Tangent stiffness**: Linearization of the residual w.r.t. \(\mathbf{U}\) gives  
  \(\mathbf{K}_T = \mathbf{K}_0 + \mathbf{K}_\sigma\):
  - **\(\mathbf{K}_0\) (material stiffness)**: From linearizing \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) w.r.t. \(\mathbf{u}\) using the **linear** part of the strain–displacement operator (same as the linear element stiffness).
  - **\(\mathbf{K}_\sigma\) (geometric stiffness)**: Depends on the **current** section forces \(N,\, M_y,\, M_z\) and shape-function derivatives. Assembled at Gauss points from current \(N,\, M\) and \(\mathrm{d}N/\mathrm{d}x\).

**Operators**:
- `GreenLagrangeStrainOperator`: `strain_linear_part`, `strain_nonlinear_part`, `linearized_strain_displacement` (definitions above).
- `GeometricStiffnessOperator.assemble_K_sigma(N, M_y, M_z, xi, weights, dN_dx, jacobian)` → 12×12 \(\mathbf{K}_\sigma\).

## Residual and Newton–Raphson

- **Residual**: \(\mathbf{R} = \mathbf{F}_{\mathrm{ext}} - \mathbf{F}_{\mathrm{int}}(\mathbf{U})\).
- **Newton step**: Solve \(\mathbf{K}_T(\mathbf{U})\,\Delta\mathbf{U} = \mathbf{R}\); update \(\mathbf{U} \leftarrow \mathbf{U} + \Delta\mathbf{U}\) until \(\|\mathbf{R}\|\) (and optionally \(\|\Delta\mathbf{U}\|\)) is below tolerance.

## Operator names in code

| Concept | Class / method |
|--------|------------------|
| Green–Lagrange strain | `GreenLagrangeStrainOperator` (utilities/total_lagrangian_beam.py) |
| Linear / nonlinear strain | `strain_linear_part`, `strain_nonlinear_part` |
| Strain–displacement (linearized) | `linearized_strain_displacement` |
| Section forces from strain | `StressResultantOperator.section_forces_from_strain` |
| Geometric stiffness | `GeometricStiffnessOperator.assemble_K_sigma` |

## Nonlinear static results (secondary/tertiary)

For **nonlinear static** runs, the formulation cache holds the **initial** (U = 0) element evaluation. Secondary and tertiary results use the **converged** \(\mathbf{U}_{\mathrm{global}}\) (and disassembled \(\mathbf{U}_e\)) with this cache: strain and stress are computed as \(\boldsymbol{\varepsilon} = \mathbf{B}_{\mathrm{lin}}\,\mathbf{U}_e\) and \(\boldsymbol{\sigma} = \mathbf{D}\,\boldsymbol{\varepsilon}\). Thus they reflect a **linearized** strain/stress at the converged displacement; the Gauss point data in the cache do not represent the full Green–Lagrange strain at the last configuration unless the cache is updated from the last converged element evaluation (future improvement).

## Element formulation logging

When `logger_operator` is set (e.g. when a job is run with `job_results_dir` and the standard formulation log directories exist), **Euler–Bernoulli-3D-Nonlinear** and **Timoshenko-3D-Nonlinear** elements write the same level of tensor detail as the linear beam elements. Stiffness logs include element length \(L\), material matrix \(\mathbf{D}\), and per–Gauss-point blocks: \(\xi\), \(x\), weight \(w\), shape-function derivatives \(\mathrm{d}N/\mathrm{d}\xi\) and \(\mathrm{d}^2N/\mathrm{d}\xi^2\), strain–displacement matrix \(\mathbf{B}\), and contribution \(\mathbf{B}^\top\mathbf{D}\mathbf{B}\). Force logs include distributed-load blocks (shape functions and load at each GP, total \(\mathbf{F}_e\)) and point-load blocks (force, shape functions, translation/rotation contributions), then the final force vector. Log files are written under `element_stiffness_matrices` and `element_force_vectors` in the job results directory.

## References

Standard Total Lagrangian beam formulations (e.g. Bathe, Crisfield) for the form of \(\mathbf{E}_{\mathrm{GL}}\) and \(\mathbf{K}_\sigma\).
