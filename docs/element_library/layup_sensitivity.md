# Sensitivities to ply angles and thicknesses

With `layup_sensitivity.txt` in a static or eigen job, 1d-fem writes the derivatives of
chosen responses (nodal displacements, the compliance, natural frequencies) with respect to
chosen ply angles and ply thicknesses of the thin-walled sections. The derivatives are exact
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

## Beam

The Timoshenko element is linear in the section's 6 x 6 stiffness and mass, so the element
derivatives \(\partial K_e\) and \(\partial M_e\) are the element's own matrices formed with
the derivative section in place of its own. From the solved job:

- **static:** \(K U = F\) with \(F\) independent of the layup. For a displacement
  \(J = g^T U\), \(\partial J = -\lambda^T \partial K\,U\) with \(K\lambda = g\) on the free
  DOFs. The compliance \(C = F^T U\) is self-adjoint: \(\partial C = -U^T \partial K\,U\).
- **eigen:** \(\partial(\omega^2) = \varphi^T(\partial K - \omega^2 \partial M)\varphi\,/\,\varphi^T M \varphi\)
  for each (simple) mode, and \(\partial f = \partial(\omega^2) / (8\pi^2 f)\).

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
```

- `ply` counts from 1 at the bottom of the laminate, as `[ply]` in `laminate.txt`; `kind` is
  `angle` (derivative per degree) or `thickness` (per metre).
- `displacement` and `compliance` responses are for static jobs, `frequency` (mode number
  from 1, in Hz) for eigen jobs.

The job writes `secondary_results/sensitivity/layup_sensitivity.csv` with columns
`parameter, laminate_id, ply, kind, value, response, response_value, derivative,
relative_derivative`. The parameter label is `L{laminate_id}.p{ply}.{kind}`. For thickness
parameters `relative_derivative` is the logarithmic derivative \((t / J)\,\partial J/\partial t\).

## Checks

`tests/test_layup_sensitivity.py` checks the section stiffness and mass derivatives and the
job derivatives (static displacements, twist, compliance; eigen frequencies) against
central differences of re-solved jobs.

## Example

[examples/thin_walled/layup_sensitivity](../../examples/thin_walled/layup_sensitivity/README.md):
a bend-twist coupled CFRP box spar. Every derivative of the tip deflection, tip twist,
compliance and first four frequencies with respect to nine ply angles and thicknesses,
checked against central differences, and a Newton loop on the off-axis ply angle that
reaches a target twist per unit flap deflection.

## Limits

- Linear static and eigen jobs with `LinearTimoshenkoBeamElement3D` elements on thin-walled
  sections. Other element types, 2D sections, pretwisted and tapered elements, second-order
  and rotating (centrifugal) jobs are refused.
- Eigenvalue derivatives assume simple eigenvalues. For repeated frequencies (a doubly
  symmetric section) the derivative of each of them depends on the chosen modes.
- Loads are taken as independent of the layup, so gravity and inertial loads, whose
  magnitude depends on the ply thicknesses, are not differentiated.
