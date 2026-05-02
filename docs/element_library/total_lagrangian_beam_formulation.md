# Total Lagrangian beam formulation (geometric nonlinearity)

This document describes the geometrically nonlinear 2-node 3D beam formulation used for **Euler-Bernoulli-3D-Nonlinear** and **Timoshenko-3D-Nonlinear** elements. The formulation is **Total Lagrangian**: all quantities are referred to the **initial (undeformed) configuration**.

**Operators in code:** Green–Lagrange strain and \\(\\mathbf{B}_\\mathrm{lin}\\)/\\(\\mathbf{B}_\\mathrm{nl}\\) are implemented in ``pre_processing/element_library/nonlinear/{euler_bernoulli,timoshenko}/utilities/green_lagrange_strain.py``; section resultants \\(\\mathbf{S}=\\mathbf{D}\\mathbf{E}\\) (2PK-type, StVK beam reduction) in ``stress_resultant.py``; geometric tangent \\(\\mathbf{K}_\\sigma\\) in ``geometric_stiffness.py``. Linear infinitesimal theory uses ``linear/.../utilities/B_matrix.py`` and ``D_matrix.py`` instead.

## Scope

Implementations are **1D beam elements in 3D space** (the library’s line mesh). There are **no** shell elements or separate 2D finite-element meshes of the cross-section.

## Frozen kinematic summary

| Item | Euler–Bernoulli NL | Timoshenko NL |
|------|---------------------|----------------|
| Linear strain \(\mathbf{E}_{\mathrm{lin}}\) | Hermite-based \(\varepsilon_x\), curvature from transverse \(u_y,u_z\) / rotations | \(\kappa_y=\partial\theta_y/\partial x\), \(\kappa_z=\partial\theta_z/\partial x\), shear \(\gamma_{xy},\gamma_{xz}\), \(\phi_x\) — same **linear** `B` as linear Timoshenko |
| \(\mathbf{E}_{\mathrm{nl}}\) | Axial Green–Lagrange quadratics; \(\kappa\) axial–curvature coupling (`GreenLagrangeStrainOperator`, EB utility) | Same axial row; **rotation-based** \(\kappa\) coupling \(\kappa_{y,\mathrm{nl}} \approx u_x'\,\partial^2\theta_y/\partial x^2\) (and \(\theta_z\) analogue); centroid nonlinear shear remainders when shear NL is on |
| Equilibrium \(\mathbf{F}_{\mathrm{int}}\) | \(\sum_g \mathbf{B}_{\mathrm{tot}}^\top \mathbf{S}\, w_g \det J\), \(\mathbf{B}_{\mathrm{tot}}=\mathbf{B}_{\mathrm{lin}}+\mathbf{B}_{\mathrm{nl}}\) | Same structure using Timoshenko \(\mathbf{B}_{\mathrm{lin}}\) from `StrainDisplacementOperator` |

## Reference configuration

- **Reference**: The **initial** geometry (element length \(L\), local axis from the grid) is fixed. No moving frame; the same local frame as in the linear formulation is used.
- **Coordinate mapping**: \(x(\xi) = \frac{1-\xi}{2} x_1 + \frac{1+\xi}{2} x_2\), \(\frac{\mathrm{d}x}{\mathrm{d}\xi} = L/2\), \(\frac{\partial\xi}{\partial x} = 2/L\).

## Strain: Green–Lagrange

Strain vector \(\mathbf{E} = [\varepsilon_x,\, \kappa_y,\, \kappa_z,\, \gamma_{xy},\, \gamma_{xz},\, \phi_x]^\top\) with **nonlinear** terms in displacement:

- **Axial (Green–Lagrange)**  
  \(\varepsilon_x = \frac{\partial u_x}{\partial x} + \frac{1}{2}\left(\frac{\partial u_x}{\partial x}\right)^2 + \frac{1}{2}\left(\frac{\partial u_y}{\partial x}\right)^2 + \frac{1}{2}\left(\frac{\partial u_z}{\partial x}\right)^2\).

- **Bending** \(\kappa_y\), \(\kappa_z\): linear part as in the **corresponding linear** beam; nonlinear supplements from `GreenLagrangeStrainOperator` (EB vs Timoshenko utilities differ — see table above).

- **Shear** \(\gamma_{xy}\), \(\gamma_{xz}\): zero for Euler–Bernoulli; for Timoshenko, linear relation as in linear Timoshenko plus optional centroid nonlinear remainder in `strain_nonlinear_part` when `include_shear=True`.

- **Torsion** \(\phi_x = \frac{\partial\theta_x}{\partial x}\) (linear).

So \(\mathbf{E} = \mathbf{E}_{\mathrm{lin}}(\mathbf{u}) + \mathbf{E}_{\mathrm{nl}}(\mathbf{u})\).

## Stress and section forces

- **Stress**: Second Piola–Kirchhoff stress in the reference configuration; constitutive \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) (same \(\mathbf{D}\) as in linear theory).
- **Section forces**: From \(\mathbf{S} = \mathbf{D}\,\mathbf{E}\) at each Gauss point we obtain **axial force \(N\)** and **moments \(M_y\), \(M_z\)** (and torsion/shear as in the linear 6-component beam). These are the stress resultants that drive the geometric stiffness.

**Operator**: `StressResultantOperator.section_forces_from_strain(E, D)` → \((N,\, M_y,\, M_z)\).

## Internal force and tangent stiffness

### Euler–Bernoulli NL

- **Internal force**: \(\mathbf{F}_{\mathrm{int}} = \sum_g \mathbf{B}_{\mathrm{tot}}^\top \mathbf{S}\, w_g \det J\) with \(\mathbf{S}=\mathbf{D}\mathbf{E}\), \(\mathbf{B}_{\mathrm{tot}}=\mathbf{B}_{\mathrm{lin}}+\mathbf{B}_{\mathrm{nl}}\) from `GreenLagrangeStrainOperator`.
- **Tangent**: \(\mathbf{K}_T = \mathbf{K}_{\mathrm{mat}} + \mathbf{K}_\sigma\) with \(\mathbf{K}_{\mathrm{mat}} = \sum_g \mathbf{B}_{\mathrm{tot}}^\top \mathbf{D}\,\mathbf{B}_{\mathrm{tot}}\, w_g \det J\) and \(\mathbf{K}_\sigma\) from `GeometricStiffnessOperator`.

### Timoshenko NL

- **Internal force**: Same weak-form structure with \(\mathbf{B}_{\mathrm{lin}}\) the **physical** linear Timoshenko `B` matrix and \(\mathbf{B}_{\mathrm{nl}}\) from `nonlinear.timoshenko.utilities.green_lagrange_strain`.
- **Tangent**: \(\mathbf{K}_T = \mathbf{K}_0 + \mathbf{K}_{\delta} + \mathbf{K}_\sigma\):
  - \(\mathbf{K}_0\): selective assembly `assemble_timoshenko_K0` — **same** linear Timoshenko stiffness as the linear element.
  - \(\mathbf{K}_{\delta}\): \(\sum_g \bigl(\mathbf{B}_{\mathrm{tot}}^\top \mathbf{D}\,\mathbf{B}_{\mathrm{tot}} - \mathbf{B}_{\mathrm{lin}}^\top \mathbf{D}\,\mathbf{B}_{\mathrm{lin}}\bigr)\, w_g \det J\) on the TL Gauss loop so \(\mathbf{K}_{\delta}=\mathbf{0}\) at \(\mathbf{U}=\mathbf{0}\) and \(\mathbf{K}_T\) matches \(\mathbf{K}_0\) initially (with \(\mathbf{K}_\sigma=\mathbf{0}\)).
  - \(\mathbf{K}_\sigma\): geometric stiffness from section forces as for EB.

**Operators**:
- Canonical strain utilities: `nonlinear.euler_bernoulli.utilities.green_lagrange_strain`, `nonlinear.timoshenko.utilities.green_lagrange_strain` (prefer these over `utilities/total_lagrangian_beam.py` stubs).
- `GeometricStiffnessOperator.assemble_K_sigma(N_gp, M_y_gp, M_z_gp, weights, dN_dx, jacobian)` → 12×12 \(\mathbf{K}_\sigma\).

## Residual and Newton–Raphson

- **Residual**: \(\mathbf{R} = \mathbf{F}_{\mathrm{ext}} - \mathbf{F}_{\mathrm{int}}(\mathbf{U})\).
- **Newton step**: Solve \(\mathbf{K}_T(\mathbf{U})\,\Delta\mathbf{U} = \mathbf{R}\); update \(\mathbf{U} \leftarrow \mathbf{U} + \Delta\mathbf{U}\) until \(\|\mathbf{R}\|\) (and optionally \(\|\Delta\mathbf{U}\|\)) is below tolerance.

## Operator names in code

| Concept | Class / method |
|--------|------------------|
| Green–Lagrange strain | `GreenLagrangeStrainOperator` in EB or Timoshenko `utilities/` |
| Linear / nonlinear strain | `strain_linear_part`, `strain_nonlinear_part` |
| Strain–displacement (linearized / NL gradient) | `linearized_strain_displacement`, `nonlinear_strain_displacement_gradient` |
| Section forces from strain | `StressResultantOperator.section_forces_from_strain` |
| Geometric stiffness | `GeometricStiffnessOperator.assemble_K_sigma` |

## Nonlinear static results (secondary/tertiary)

For **nonlinear static** runs, the formulation cache may hold the **initial** (U = 0) element evaluation. Secondary and tertiary results use the **converged** \(\mathbf{U}_{\mathrm{global}}\) with this cache: strain and stress are sometimes computed as \(\boldsymbol{\varepsilon} = \mathbf{B}_{\mathrm{lin}}\,\mathbf{U}_e\) and \(\boldsymbol{\sigma} = \mathbf{D}\,\boldsymbol{\varepsilon}\). Thus they may reflect a **linearized** strain/stress at the converged displacement unless the cache is updated from the last converged element evaluation with full \(\mathbf{E}\).

## Element formulation logging

When `logger_operator` is set (e.g. when a job is run with `job_results_dir` and the standard formulation log directories exist), **Euler–Bernoulli-3D-Nonlinear** and **Timoshenko-3D-Nonlinear** elements write the same level of tensor detail as the linear beam elements. Stiffness logs include element length \(L\), material matrix \(\mathbf{D}\), and per–Gauss-point blocks: \(\xi\), \(x\), weight \(w\), shape-function derivatives \(\mathrm{d}N/\mathrm{d}\xi\) and \(\mathrm{d}^2N/\mathrm{d}\xi^2\), strain–displacement matrix \(\mathbf{B}\), and contribution \(\mathbf{B}^\top\mathbf{D}\mathbf{B}\). Force logs include distributed-load blocks (shape functions and load at each GP, total \(\mathbf{F}_e\)) and point-load blocks (force, shape functions, translation/rotation contributions), then the final force vector. Log files are written under `element_stiffness_matrices` and `element_force_vectors` in the job results directory.

## References

Standard Total Lagrangian beam formulations (e.g. Bathe, Crisfield) for the form of \(\mathbf{E}_{\mathrm{GL}}\) and \(\mathbf{K}_\sigma\).
