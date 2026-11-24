# SectionAnalysis: section stiffness for every element

`anisotropic_section_model/section_analysis.py` defines `SectionAnalysis`. Every element can read it as `self.section_analysis`.

It holds:

- the 6 × 6 section stiffness, in the element Voigt order `[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x]`;
- the shear-rigid 4 × 4 `classical_stiffness` for `[eps_x, kappa_y, kappa_z, phi_x]`;
- \(E\Gamma\), the shear and tension centres, and the mass properties.

## Sources

| Source | When | Stiffness |
|---|---|---|
| `"thin_walled"` | element listed in `element_section.txt` | `ThinWalledSectionAnalysis`, mapped to element signs |
| `"isotropic"` | otherwise | the provider's `SectionAnalysis.from_isotropic`: `diag(EA, EI_y, EI_z, kappa GA, kappa GA, GJ_t)` about the shear centre, from `section.txt` and `material.txt` (`kappa` defaults to 5/6) |

1D-FEM computes no isotropic section properties of its own. `element_library/isotropic_section.py` makes the one `from_isotropic` call, and every isotropic element reads it: the bar and truss take their rows of its `D`, the Euler-Bernoulli elements its `euler_bernoulli_D()`, the Timoshenko, nonlinear and corotational beams its full `D`, and all of them their consistent mass from its mass per length. The twist DOFs carry the polar moment \(\rho I_p\): `I_x` from `section.txt`, or \(I_y + I_z\) when `I_x` is `0` or absent. A shear centre `y_sc`, `z_sc` off the centroid couples shear and twist about it: `D[3,5] = -z_sc κGA`, `D[4,5] = y_sc κGA`, `D[5,5] = G J_t + κGA (y_sc² + z_sc²)`, so a shear force through the shear centre twists nothing. Axial and bending stay uncoupled, the shear-rigid `euler_bernoulli_D()` keeps `G J_t`, and the bar and truss take no offset. `tests/test_isotropic_sections_from_provider.py` checks each element against the provider and swaps in a provider that doubles every stiffness to show that no element keeps its own formulas.

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
