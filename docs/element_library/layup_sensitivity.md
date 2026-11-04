# Sensitivities to ply angles and thicknesses

With `layup_sensitivity.txt` in a static or eigen job, 1d-fem writes the derivatives of
chosen responses (nodal displacements, the compliance, natural frequencies) with respect to
chosen ply angles and ply thicknesses of the thin-walled sections (for Timoshenko, Benscoter
and GBT beams), and ply angles of sections meshed in 2D. The derivatives are exact
for the discrete model: the section derivatives are found by differentiating the section
analysis, not by finite differences, and the beam derivatives follow by the adjoint method
(static) or the eigenvalue derivative (eigen). One solve gives every derivative.

## Section

### Laminate

A ply angle enters the laminate only through the rotated stiffness of that ply,
\(\bar Q = T^T Q\,T\), so

\[
\frac{\partial \bar Q}{\partial \theta} = \frac{\partial T}{\partial \theta}^T Q\,T + T^T Q\,\frac{\partial T}{\partial \theta},
\]

with \(\partial T / \partial \theta\) in closed form and the result per degree. A thickness
change of ply \(k\) moves every interface \(z_j\) of the laminate by
\(\partial z_j / \partial t_k = [j > k] - \tfrac12\) (the laminate stays centred on its wall
midline), so

\[
\frac{\partial}{\partial t_k}\,\frac1p \sum_i \bar Q_i \left(z_{i+1}^p - z_i^p\right)
= \sum_i \bar Q_i \left(z_{i+1}^{p-1}\,\frac{\partial z_{i+1}}{\partial t_k} - z_i^{p-1}\,\frac{\partial z_i}{\partial t_k}\right)
\]

for \(A\), \(B\) and \(D\) (\(p = 1, 2, 3\)) and the transverse shear stiffness, and the
same for the mass moments \(\int \rho\, z^n\,dz\). The wall law is the laminate condensed to
the wall's strains, \(K_c = P^T \mathrm{ABD}\,P\) with \(P = [I;\ -K_{rr}^{-1}K_{rk}]\), so
\(\partial K_c = P^T\,\partial \mathrm{ABD}\,P\) (the terms in \(\partial P\) cancel because
the condensed resultants are zero). With `hoop = restrained` the same holds on the smaller
set of condensed strains.

### Section stiffness and mass

The section analysis (`ThinWalledSectionAnalysis`) solves linear systems whose matrices are
linear in the wall laws, so the derivative of its flexibility \(F\) is found by
differentiating those systems directly:

\[
\partial F = \sum_e w_e\left(\partial\Phi_e^T K_e \Phi_e + \Phi_e^T K_e\,\partial\Phi_e + \Phi_e^T\,\partial K_e\,\Phi_e\right),
\qquad \partial S = -S\,\partial F\,S,
\]

where \(\Phi\) are the warping solutions and \(\partial\Phi\) solves the same system with the
right-hand side from \(\partial K\). The warping constraint multipliers drop out. The mass
matrix is linear in the mass moments, so its derivative is exact by construction.

`section_derivative(analysis, LayupParameter(laminate_id, ply, kind))`
(`pre_processing/section_integration/thin_walled/layup_sensitivity.py`) returns
`SectionDerivative` with `stiffness`, `flexibility`, `mass` and `mass_matrix`. Against
central differences the stiffness derivative agrees to about \(10^{-9}\), with the hoop
free or restrained.

### Benscoter stiffness and warping inertia

With `benscoter=True` the derivative also carries the 8 x 8 Benscoter stiffness and the
warping inertia. The torsional warping \(\omega\) is one more solve with the same wall-law
matrices (\(E\,\omega = -R_\phi\)), and its rigid part is removed by a 3 x 3 solve in the walls'
axial stiffness; both are differentiated. Then

\[
\partial K_B = \sum_e w_e\left(\partial G_e^T K_e G_e + G_e^T K_e\,\partial G_e + G_e^T\,\partial K_e\,G_e\right),
\qquad
\partial I_\omega = \sum_e w_e\left(\partial m_e\,\omega_e^2 + 2 m_e\,\omega_e\,\partial\omega_e\right),
\]

with \(G\) the wall strains of the beam strains, \(\omega_s\,d\) and \(\omega\,\psi'\). Both
agree with central differences to about \(10^{-9}\).

### Second-order terms

With `geometric=True` it carries the derivative of `geometric_terms` \(H\), the wall
resultants of the Saint-Venant stresses against the section geometry. \(H\) is linear in the
wall law and in the warping fields, so \(\partial H\) takes \(\partial K_e\) and \(\partial\Phi_e\).

### Sections meshed in 2D

A 2D mesh made by `SectionMesh.from_thin_walled` tags each element with its laminate ply
(`SectionMesh.plies`, `[Section_Element_Ply]` in `section_mesh_2d.txt`). An `angle`
parameter turns the tagged elements together. The element stiffness is
\(C = T C_m T^T\) with \(T\) the bond matrix of the material axes \(R\); \(T\) is bilinear in
\(R\), so \(\partial T = T(\partial R, R) + T(R, \partial R)\) and
\(\partial C = \partial T\,C_m T^T + T C_m\,\partial T^T\). The Giavotto saddle systems of
`SolidSectionAnalysis` are linear in \(C\) and are differentiated as above. With
`benscoter=True` so are the torsional warping solve, its shear-centre shift in the axial
modulus \(E_x = C_{xx} - C_{xi} C_{ii}^{-1} C_{ix}\), and the relaxed-gradient solve behind
\(K_{77}\); the second-order terms and axial weights are linear in the stresses and \(E_x\).
The mass does not depend on ply angles.

`solid_section_derivative(analysis, LayupParameter(laminate_id, ply, "angle"))`
(`pre_processing/section_integration/solid/layup_sensitivity.py`) returns the same
`SectionDerivative`, and jobs on `section_mesh_2d.txt` sections use it. A thickness
parameter moves the mesh nodes and is refused.

### GBT sections

`LinearGBTBeamElement3D` is linear in the GBT section matrices: the stiffness
\(\sum w\,H^T K H\), the mass and the geometric tensors, with \(H\) the wall strains of the
beam strains, the warping modes and the distortional modes, and \(K\) the GBT wall law
(\(N_{ss}\) released). The modes change with the layup as well, and holding them fixed is
not good enough: on the box cantilever with two distortional modes it puts the derivative of
the first (distortional) frequency with respect to a ply thickness 18 % off. So
`gbt_section_derivative` differentiates them:

- the distortional modes \(x_k\) are eigenvectors of \(B x = \lambda C x\) on the admissible
  frame space, which the other eigenvectors and the rigid modes complete. With
  \(r_k = (\partial B - \lambda_k\,\partial C) x_k\), the component of \(\partial x_k\) on another mode
  \(x_j\) is \(-x_j^T r_k / ((\lambda_j - \lambda_k)\,x_j^T C x_j)\), and its rigid part solves
  \(R^T C R\,c = -R^T\,\partial C\,x_k\);
- the Saint-Venant warping modes span the warping fields \(W_0, W_1\) with the rigid axial
  fields removed in the axial stiffness \(A\), which differentiates with \(\partial W_0\),
  \(\partial W_1\) and \(\partial A\);
- the shear-lag modes solve \(E x = \mu A x\) A-orthogonal to the rigid and Saint-Venant
  fields, and are differentiated like the distortional ones, with the multipliers of that
  constraint.

Components of a derivative within the span of the modes kept are left out: they only
change the basis of a block, which a GBT beam does not see because its amplitudes are free.
That holds for a beam whose GBT elements share one section; on beams with several sections
the bases are aligned between neighbours (`gbt_alignment`), the alignment is not
differentiated, and the job is refused.

## Beam

The Timoshenko element is linear in the section's 6 x 6 stiffness and mass, and the
Benscoter element in its 8 x 8 stiffness, mass and warping inertia, so the element
derivatives \(\partial K_e\) and \(\partial M_e\) are the element's own matrices formed with
the derivative section in place of its own. From the solved job:

- **static:** \(K U = F\) with \(F\) independent of the layup. For a displacement
  \(J = g^T U\), \(\partial J = -\lambda^T \partial K\,U\) with \(K\lambda = g\) on the free
  DOFs. The compliance \(C = F^T U\) is self-adjoint: \(\partial C = -U^T \partial K\,U\).
- **eigen:** \(\partial(\omega^2) = \varphi^T(\partial K - \omega^2 \partial M)\varphi\,/\,\varphi^T M \varphi\)
  for each (simple) mode, and \(\partial f = \partial(\omega^2) / (8\pi^2 f)\).
- **linear buckling:** \((K + \lambda K_g(U_0))\varphi = 0\) about the prestress \(K U_0 = F\), so
  \(\partial\lambda = -\varphi^T(\partial K + \lambda\,\partial K_g)\varphi\,/\,\varphi^T K_g\varphi\).
  \(K_g\) changes with the section at fixed \(U_0\): its section forces \(D e\) and its
  second-order terms \(H\) enter linearly, so that part is the element's \(K_g\) formed with
  \(\partial D\) plus the one formed with \(\partial H\). It also changes through the prestress:
  \(K_g\) is linear in \(U_0\), so \(\varphi^T K_g(\partial U_0)\varphi = g^T \partial U_0\) with
  \(g_i = \varphi^T K_g(e_i)\varphi\) (element by element), and one adjoint solve \(K\mu = g\) gives
  \(g^T\partial U_0 = -\mu^T \partial K\,U_0\).

**Pretwisted, curved and tapered elements.** The pretwist coupling (`with_pretwist`) is
linear in the walls' axial stiffness weights, so a pretwisted element formed with the
derivative section (with derivative weights) gives its derivative. A curved element is a
chain of straight sub-elements condensed by their static shapes \(\Phi\):
\(K = \Phi^T K_c \Phi\) is stationary in \(\Phi\), so \(\partial K = \Phi^T \partial K_c \Phi\), while
its mass takes the change of the shapes, \(\partial\Phi_i = -K_{ii}^{-1}(\partial K_c\Phi)_i\). The
taper coupling is \(\sum_l A_l\,G_l\) over the laminates, with \(A_l\) the axial wall stiffness
and \(G_l\) a geometric factor that `section_taper` now writes per laminate
(`[Section_Taper_Laminate]` in `section_taper.txt`, written by `blade_beam_job`), so its
derivative is \(\partial A_l\,G_l\).

The products are formed element by element (\(\lambda_e^T \partial K_e U_e\)), so no
global derivative matrix is assembled.

## In a job

`layup_sensitivity.txt` in the job directory:

```
[Layup_Parameter]
[laminate_id] [ply] [kind]
1 5 angle
1 5 thickness
[Response]
[type] [node_id] [dof] [mode]
displacement 12 UZ -
displacement 12 RX -
compliance - - -
frequency - - 1
buckling - - 1
```

- `ply` counts from 1 at the bottom of the laminate, as `[ply]` in `laminate.txt`; `kind` is
  `angle` (derivative per degree) or `thickness` (per metre).
- `displacement` and `compliance` responses are for static jobs, `frequency` (mode number
  from 1, in Hz) for eigen jobs, `buckling` (the load factor of a mode, from 1) for linear
  buckling jobs with a linear static prestress.

The job writes `secondary_results/sensitivity/layup_sensitivity.csv` with columns
`parameter, laminate_id, ply, kind, value, response, response_value, derivative,
relative_derivative`. The parameter label is `L{laminate_id}.p{ply}.{kind}`. For thickness
parameters `relative_derivative` is the logarithmic derivative \((t / J)\,\partial J/\partial t\).

## Checks

`tests/test_layup_sensitivity.py` checks the section stiffness, mass, Benscoter and
second-order derivatives, and the job derivatives against central differences of re-solved
jobs: static displacements, twist and compliance and eigen frequencies with Timoshenko and
Benscoter elements, and the first two buckling load factors of a box cantilever under
compression, shear, torque and bending, and static and eigen responses of a pretwisted,
prebent box and of a tapered, pretwisted one, and static (Timoshenko), eigen (Benscoter) and
buckling responses of the box meshed in 2D, for a ply angle (all to \(2\times10^{-5}\)), and
static, eigen and buckling responses of GBT cantilevers with two distortional modes, with and
without two shear-lag modes (to \(3\times10^{-6}\), and \(2\times10^{-5}\) with the shear-lag modes). On the IEA 15 MW
outboard segment (tapered, pretwisted, 6 elements) the cap angle derivatives of tip
deflection, twist, compliance and frequencies agree with re-built jobs to \(10^{-5}\). The web
core thickness derivatives agree to \(4\times10^{-4}\): the derivative holds the wall
mid-lines, while a re-cut loft moves the walls with the thickness.

## Example

[examples/thin_walled/layup_sensitivity](../../examples/thin_walled/layup_sensitivity/README.md):
a bend-twist coupled CFRP box spar. Every derivative of the tip deflection, tip twist,
compliance and first four frequencies with respect to nine ply angles and thicknesses,
checked against central differences, and a Newton loop on the off-axis ply angle that
reaches a target twist per unit flap deflection.

## Limits

- Linear static, eigen and linear buckling (linear static prestress) jobs with
  `LinearTimoshenkoBeamElement3D`, `LinearBenscoterBeamElement3D` (not buckling) or
  `LinearGBTBeamElement3D` elements. Timoshenko and Benscoter elements take thin-walled
  sections or, for ply angles only, sections meshed in 2D with ply tags; pretwisted, curved
  (not in buckling) and tapered elements are included (a tapered one needs
  `[Section_Taper_Laminate]`). GBT elements need one thin-walled section along the beam
  and straight elements without precurvature, pretwist or wall loads. Other element types
  (third-order, corotational), ply thicknesses of 2D sections, 2D meshes without ply tags
  (blade cuts from `section_mesh_from_cut`), nonlinear prestress, second-order static and
  rotating (centrifugal) jobs are refused.
- Wall mid-lines are held fixed: a thickness change grows the stack about its mid-surface.
- Eigenvalue derivatives assume simple eigenvalues. For repeated frequencies (a doubly
  symmetric section) the derivative of each of them depends on the chosen modes.
- Loads are taken as independent of the layup, so gravity and inertial loads, whose
  magnitude depends on the ply thicknesses, are not differentiated.
