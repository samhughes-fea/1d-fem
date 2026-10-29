# Wall and ply stress recovery for thin-walled sections

Every beam element with a thin-walled section (`element_section.txt`) can give
its wall resultants, ply stresses and interlaminar stresses at any point of the
section. Code:

- `pre_processing/section_integration/thin_walled/stress_recovery.py`: the
  wall-level steps and `SaintVenantRecovery` (section forces in, stresses out).
- `pre_processing/section_integration/thin_walled/element_recovery.py`: the
  per-element wall strains (`element_wall_field`, `element_wall_stress`).
- `post_processing/thin_walled/wall_stresses.py`: the CSV written by
  `process_job`.

## 1. Wall strains a beam theory carries

| Element | Kept wall strains | Their rate along the beam |
|---|---|---|
| Euler-Bernoulli | Saint-Venant strains of `[N, M_y, M_z, V_y, V_z, T]`, with `V_z = M_y'` and `V_y = -M_z'` from the Hermite curvature | beam equilibrium, `F' = T_RATE F` |
| Timoshenko | Saint-Venant strains of `F = K psi` | beam equilibrium |
| Benscoter | Saint-Venant strains plus `eps_xx += omega psi'`, `gamma_xs += omega_s (psi - phi)` | equilibrium for the Saint-Venant part, the element interpolation for the warping part |
| GBT | plane-section, warping-mode and distortional-mode strains, including the hoop curvature `kappa_ss` | the element interpolation |

The kept strains are `[eps_xx, gamma_xs, kappa_xx, kappa_xs]`, plus `kappa_ss`
for GBT. With no distributed load the Saint-Venant rate is exact: the wall
strains are linear in the section forces and `T_RATE` is nilpotent, so
`d/dx` of the wall strains is the wall strains of `T_RATE F`.

## 2. Full mid-surface strains

The wall law is `[N; M] = ABD [eps0; kappa]` in the order
`[eps_xx, eps_ss, gamma_xs, kappa_xx, kappa_ss, kappa_xs]`. The section is free
to deform in its own plane, so `N_ss = 0`. For the hoop bending row:

- `hoop="free"` (Saint-Venant, Benscoter by default): `M_ss = 0`.
- `hoop="restrained"`: `kappa_ss = 0`, or the GBT value. This adds the plate
  factor `1 / (1 - nu^2)` to wall bending.

Solving the released rows gives the six strains, then the resultants
`[N_xx, N_ss, N_xs, M_xx, M_ss, M_xs] = ABD e`.

## 3. Ply stresses

Through the thickness `eps(n) = eps0 + n kappa`, with `n` along the wall normal
`n = x × s` (on a wall running in `-y`, `n` points in `-z`). Each ply gives
`[sigma_xx, sigma_ss, tau_xs] = Q_bar eps` in wall axes and
`[sigma_11, sigma_22, tau_12] = Q T_eps eps` in ply axes.

## 4. Interlaminar stresses

Integrating 3D equilibrium up from the traction-free bottom face `n = -h/2`:

```
tau_xn(n)   = -∫ (d sigma_xx/dx + d tau_xs/ds) dn
tau_sn(n)   = -∫ (d tau_xs/dx  + d sigma_ss/ds) dn
sigma_nn(n) = -∫ (d tau_xn/dx  + d tau_sn/ds) dn
```

Each ply's in-plane stresses are linear in `n`, so the integrals are exact
polynomials per ply (`numpy.polynomial`). The `s` derivatives come from a
three-point stencil inside the wall element that holds `s`. At the top face
the values are the equilibrium residuals: `-(N_xx,x + N_xs,s)` and
`-(N_xs,x + N_ss,s)` vanish for a Saint-Venant field, and
`∫ tau_xn dn = Q_x`, `∫ tau_sn dn = Q_s` are the plate transverse shears.

## 5. Output

When a linear static job has thin-walled sections, `process_job` writes
`secondary_results/wall_stresses/wall_stresses.csv` with one row per element,
wall point (`points_per_wall = 5` along each wall), ply and ply face, at each
element centre:

```
element_id, x, wall_id, s, y, z,
N_xx, N_ss, N_xs, M_xx, M_ss, M_xs,
ply, n, sigma_xx, sigma_ss, tau_xs, sigma_11, sigma_22, tau_12,
tau_xn, tau_sn, sigma_nn
```

Elements without wall recovery are skipped.

## 6. Python use

```python
from pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

rec = SaintVenantRecovery(analysis)                 # a ThinWalledSectionAnalysis
p = rec.at(wall_id=3, s=0.3, forces=[N, M_y, M_z, V_y, V_z, T], points_per_ply=5)
p.resultants           # [N_xx, N_ss, N_xs, M_xx, M_ss, M_xs]
p.plies[0].stress_ply  # [sigma_11, sigma_22, tau_12] at the bottom of ply 1
p.plies[0].tau_xn
```

For a solved element use `element_wall_stress(element, U_e, xi, analysis, wall_id, s)`.

## 7. Checks

`tests/test_thin_walled_stress_recovery.py`:

- box bending `sigma = M z / I`, Bredt torsion `q = T / 2A` and the web shear flow;
- released resultants are zero and the ply stresses integrate to the resultants;
- interlaminar top-face residual about `1e-10`, `∫ tau_xn = Q_x` and the
  parabolic `1.5 Q / t` peak of a single-ply wall;
- every element type against the exact cantilever forces, including
  Benscoter and GBT with free warping matching Saint-Venant;
- `process_job` writes the CSV and it matches the direct recovery.

The worked examples in `examples/thin_walled/` compare against Saint-Venant tip
deflections and multi-cell Bredt-Batho shear flows. A comparison against a
shell model (the plan's 3 % check) needs a shell solver and is not in this
repository.
