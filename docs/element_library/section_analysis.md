# SectionAnalysis: section stiffness for every element

`pre_processing/section_integration/section_analysis.py` defines `SectionAnalysis`. Every element can read it as `self.section_analysis`.

It holds:

- the 6 × 6 section stiffness, in the element Voigt order `[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]`;
- the shear-rigid 4 × 4 `classical_stiffness` for `[eps_x, kappa_y, kappa_z, phi_x]`;
- \(E\Gamma\), the shear and tension centres, and the mass properties.

## Sources

| Source | When | Stiffness |
|---|---|---|
| `"thin_walled"` | element listed in `element_section.txt` | `ThinWalledSectionAnalysis`, mapped to element signs |
| `"isotropic"` | otherwise | `diag(EA, EI_y, EI_z, kappa GA, kappa GA, GJ_t)` from `section.txt` and `material.txt` (`kappa` defaults to 5/6) |

The isotropic form equals the linear Timoshenko element's `D` exactly. A test checks this.

## Element sign convention

The section analysis works with right-handed strains \(\boldsymbol\psi\). Probing the linear Euler-Bernoulli and Timoshenko `B` operators with right-hand-rule nodal rotations shows how the element strains relate to them:

| Row | Element strain | In right-handed terms |
|---|---|---|
| 0 | \(\varepsilon_x\) | \(u'\) |
| 1 | \(\kappa_y\) | \(-\theta_y'\) (\(= w''\) for Euler-Bernoulli) |
| 2 | \(\kappa_z\) | \(\theta_z'\) |
| 3 | \(\gamma_{xy}\) | \(v' - \theta_z\) |
| 4 | \(\gamma_{xz}\) | \(w' + \theta_y\) |
| 5 | \(\phi_x\) | \(\theta_x'\) |

So the element strains are \(\boldsymbol\varepsilon_e = \mathbf S\boldsymbol\psi\) with \(\mathbf S = \mathrm{diag}(1, -1, 1, 1, 1, 1)\). The element stiffness is \(\mathbf D_e = \mathbf S\mathbf K\mathbf S\), and the element moment \(M_y\) is \(-\int z\sigma\).

Couplings change sign accordingly. For example, a section whose tension centre is at \(z_0\) has \(D_e[0,1] = -EA\,z_0\), because \(\varepsilon_{xx} = \varepsilon - z\kappa_y\) in element terms.

Tests: `tests/test_section_analysis_element_wiring.py`.

## Elements

An element type sets `supports_section_analysis = True` once it builds its `D` from `self.section_analysis`. If any other element type is assigned a thin-walled section, it raises `NotImplementedError`.

| Element | `D` from a thin-walled section | Mass |
|---|---|---|
| `LinearEulerBernoulliBeamElement3D` | `euler_bernoulli_D()`: the shear-rigid `classical_stiffness` on rows `[0, 1, 2, 5]`. This is the inverse of the force-flexibility block, not a sub-block of the 6 × 6. | rigid-section mass with offsets, using the Hermite field interpolation |
| `LinearTimoshenkoBeamElement3D` | the full 6 × 6 `stiffness`, with selective integration that keeps every coupling ([proof](../proofs/timoshenko/coupled_selective_integration.md)) | rigid-section mass with offsets, using linear interpolation |

**Mass.** `SectionAnalysis.mass_matrix()` is the rigid-section mass per length. It includes the coupling terms from a mass centre that is off the reference line.

**Not yet supported** (both elements raise `NotImplementedError`):
- Thin-walled sections with the warping DOF (`[warping] = 1`). This comes with the Benscoter element in the next PR.

**Approximations:**
- The geometric stiffness used for buckling takes \(N\), \(M_y\) and \(M_z\) from the coupled `D`. It does not yet include the Wagner or shear-centre terms (planned for second order).

**Accuracy (Euler-Bernoulli).** Bending is cubic and twist is linear. Constant moments and torque are exact on any mesh. A shear load on a section with bending-twist coupling converges at \(O(h^2)\), because the twist rate then varies along the element.
