# GESDB weak form (locked reference for `GeometricallyExactShearDeformableBeam3D`)

This appendix fixes notation for the **12-DOF** shear-deformable beam element registered as
`GeometricallyExactShearDeformableBeam3D`. It aligns Voigt ordering with
`docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`.

## Strain vector and stress resultants

At each Gauss point along the straight reference chord (natural coordinate \(\xi \in [-1,1]\), physical \(x\)):

- **Voigt strain** \(\mathbf{E} \in \mathbb{R}^6\): axial \(\varepsilon_x\), shear \(\gamma_{xy}\), \(\gamma_{xz}\), torsion \(\kappa_x\), bending curvatures \(\kappa_y\), \(\kappa_z\) — same layout as Total Lagrangian Timoshenko in this repository.
- **Constitutive:** \(\mathbf{S} = \mathbf{D}\,(\mathbf{E} - \mathbf{E}_0)\) with \(\mathbf{D}\) from `MaterialStiffnessOperator` (diagonal Timoshenko pattern).
- **Section forces for geometric stiffness:** axial force \(N\) and bending moments \(M_y, M_z\) extracted from \(\mathbf{S}\) via `StressResultantOperator.section_forces_from_strain`.

## Discrete weak form (current production kernel)

For Total Lagrangian kinematics on the **fixed chord frame**,

\[
\mathbf{E}(\mathbf{U}_e) = \mathbf{E}_{\mathrm{lin}} + \mathbf{E}_{\mathrm{nl}}(\mathbf{U}_e),
\quad
\mathbf{E}_{\mathrm{lin}} = \mathbf{B}_{\mathrm{lin}}\mathbf{U}_e,
\]

with \(\mathbf{B}_{\mathrm{lin}}\) from `StrainDisplacementOperator.physical_coordinate_form` and
\(\mathbf{E}_{\mathrm{nl}}\) from `GreenLagrangeStrainOperator.strain_nonlinear_part` (same operators as
`NonlinearTimoshenkoBeamElement3D`).

**Internal force**

\[
\mathbf{F}_{\mathrm{int}} = \sum_g w_g\, |J|_g\, \mathbf{B}_{\mathrm{tot},g}^\top \mathbf{S}_g,
\quad
\mathbf{B}_{\mathrm{tot}} = \mathbf{B}_{\mathrm{lin}} + \mathbf{B}_{\mathrm{nl}}.
\]

**Material tangent increment**

\[
\mathbf{K}_\Delta = \sum_g w_g\, |J|_g \left(
\mathbf{B}_{\mathrm{tot},g}^\top \mathbf{D} \mathbf{B}_{\mathrm{tot},g}
- \mathbf{B}_{\mathrm{lin},g}^\top \mathbf{D} \mathbf{B}_{\mathrm{lin},g}
\right).
\]

**Geometric tangent**

\[
\mathbf{K}_\sigma = \text{assemble\_}K_\sigma(N, M_y, M_z, \ldots)
\]

via `GeometricStiffnessOperator` on the same Gauss loop.

**Total tangent:** \(\mathbf{K}_T = \mathbf{K}_0 + \mathbf{K}_\Delta + \mathbf{K}_\sigma\) with \(\mathbf{K}_0\) from selective Timoshenko assembly (`assemble_timoshenko_K0`).

## Director pathway (`gesdb_kinematics`)

The class overrides `_tl_voigt_strain_at_gauss` so **strain sampling** goes through
`gesdb_director_voigt_strain_timoshenko_12`. For the present shape functions and chord frame, that
function returns the **same** \(\mathbf{E}\) as the chord-frame Green–Lagrange reduction above — i.e.
the Simo–Vu-Quoc-style spatial strain vector collapses to this Voigt form for the chosen \(C^0\)
kinematics (see Simo & Vu-Quoc, CMAME, 1986+; Jelenić & Crisfield, 1999 for general moving-frame
forms). When a distinct director \(\mathbf{B}\) matrix is introduced, replace only the strain map;
the Gauss structure for \(\mathbf{F}_{\mathrm{int}}\) and \(\mathbf{K}_T\) remains as stated.

## TL fallback

`simulation_settings['nonlinear']['gesdb_tl_fallback'] = true` forces the parent
`_tl_voigt_strain_at_gauss` implementation (pure TL delegate) for side-by-side regression.

## Native kernel (`gesdb_kernel = native`)

When ``simulation_settings['nonlinear']['gesdb_kernel']`` is ``native`` (constructor mirrored by the factory), strains use:

- **Axial:** engineering chord stretch ``ε_ax = L_\mathrm{def}/L_\mathrm{ref} - 1`` with ``L_\mathrm{def} = \|(\mathbf{X}_2+\mathbf{u}_2)-(\mathbf{X}_1+\mathbf{u}_1)\|`` and ``L_\mathrm{ref}`` the reference chord length ``L``.
- **Remaining rows:** same linear Timoshenko rows as ``B_\mathrm{lin}\mathbf{U}`` at the Gauss point (``κ_y``, ``κ_z``, ``γ_xy``, ``γ_xz``, ``φ_x``).

The strain Jacobian ``\mathbf{B}_\mathrm{eng} = \partial \mathbf{E}/\partial \mathbf{U}_e`` has row 0 from ``∂ε_ax/∂\mathbf{U}_e`` (nonzero on translation DOFs only) and rows 1–5 equal to ``\mathbf{B}_\mathrm{lin}[1{:}6,:]``.

**Internal force:** ``\mathbf{F}_\mathrm{int} = \sum_g w_g |J| \mathbf{B}_\mathrm{eng}^\top \mathbf{S}`` with ``\mathbf{S}=\mathbf{D}(\mathbf{E}-\mathbf{E}_0)``.

**Tangent:** ``\mathbf{K}_T \approx \sum_g w_g |J| \mathbf{B}_\mathrm{eng}^\top \mathbf{D} \mathbf{B}_\mathrm{eng} + \mathbf{K}_\sigma`` using section forces from ``\mathbf{S}`` and the existing Gauss assembly for ``\mathbf{K}_\sigma``. This omits the TL ``\mathbf{K}_\delta`` split; at finite rotation native axial strain differs from chord-frame Green–Lagrange, so regressions vs TL are documented with rtol bands at moderate ``\|\mathbf{U}\|``.
