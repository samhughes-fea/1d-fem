# Thin-walled section analysis (Saint-Venant level)

`ThinWalledSectionAnalysis` turns a thin-walled section (walls of laminates, see [thin_walled_section_walls.md](thin_walled_section_walls.md)) into the fully coupled 6 × 6 Timoshenko section stiffness, plus the shear centre, tension centre, torsional warping, the Vlasov warping stiffness \(E\Gamma\) and mass properties.

Code: `pre_processing/section_integration/thin_walled/section_analysis.py`. Derivation: [../proofs/section/saint_venant_thin_walled.md](../proofs/section/saint_venant_thin_walled.md).

## What it solves

The Giavotto (1983) Saint-Venant problem, restricted to thin-wall kinematics:

- The section is rigid in its own plane.
- Each wall warps out of plane by \(w(s)\) and follows Kirchhoff plate kinematics through its thickness.
- Warping is solved for, not assumed. Open branches, closed cells and any number of cells all come out of the same equations, with no cuts, Bredt loops or shear-flow bookkeeping.
- All couplings from the wall laminates are kept: extension-twist, bending-twist, extension-shear, bending-shear.

## Conventions

Right-handed axes, rotations by the right-hand rule. The section lies in the \((y, z)\) plane and the beam axis is \(x\). Results are about the coordinate origin of the section geometry.

| Index | Strain | Definition | Force | Definition |
|---|---|---|---|---|
| 0 | \(\varepsilon\) | \(u'\) | \(N\) | \(\int\sigma\) |
| 1 | \(\kappa_y\) | \(\theta_y'\) | \(M_y\) | \(\int z\sigma\) |
| 2 | \(\kappa_z\) | \(\theta_z'\) | \(M_z\) | \(-\int y\sigma\) |
| 3 | \(\gamma_{xy}\) | \(v' - \theta_z\) | \(V_y\) | \(\int\tau_{xy}\) |
| 4 | \(\gamma_{xz}\) | \(w' + \theta_y\) | \(V_z\) | \(\int\tau_{xz}\) |
| 5 | \(\phi\) | \(\theta_x'\) | \(T\) | \(\int(y\tau_{xz} - z\tau_{xy})\) |

The axial displacement is \(u_x = u - y\theta_z + z\theta_y + \text{warping}\). The row order matches the element Voigt order \([\varepsilon_x, \kappa_y, \kappa_z, \gamma_{xy}, \gamma_{xz}, \phi_x]\). The signs of \(\theta_y\) and \(\gamma_{xz}\) are the right-hand ones. The next PR maps them onto each element's own convention.

## Wall strains

On a straight wall with tangent \(\mathbf t\), normal \(\mathbf n = (-t_z, t_y)\), \(h = y t_z - z t_y\) and \(r_n = y n_z - z n_y\):

\[
\begin{aligned}
\varepsilon_{xx} &= \varepsilon + z\kappa_y - y\kappa_z + \partial_x w \\
\gamma_{xs} &= t_y\gamma_{xy} + t_z\gamma_{xz} + h\phi + \partial_s w \\
\kappa_{xx} &= n_z\kappa_y - n_y\kappa_z - (n_y\gamma_{xy}' + n_z\gamma_{xz}' + r_n\phi') \\
\kappa_{xs} &= -2\phi
\end{aligned}
\]

- The bracket in \(\kappa_{xx}\) is the wall bending caused by transverse shear strain that varies along the beam. It is how the walls carry Kirchhoff transverse shear \(Q_x = \partial M_{xx}/\partial x\).
- \(\kappa_{xs} = -2\phi\) is the wall's own twisting. It gives the open-wall torsion term \(\sum G L t^3/3\).

## Wall law: hoop options

| `hoop` | Hoop force \(N_{ss}\) | Hoop curvature \(\kappa_{ss}\) | Use |
|---|---|---|---|
| `"free"` (default) | 0 | free (\(M_{ss}=0\)) | Saint-Venant: walls bend freely in \(s\); the textbook thin-wall values such as \(Et^3/12\) |
| `"restrained"` | 0 | 0 | the section keeps its shape, so walls do not bend in \(s\); plate stiffnesses such as \(Et^3/12(1-\nu^2)\). GBT sections always use it |

The default is `"free"` because it is the Saint-Venant answer: the 2D section solver
(`SolidSectionAnalysis`), which solves the in-plane deformation exactly, matches it
(0.2 % on a CFRP I-section, 0.8 % on a CAS box), while `"restrained"` reads `D26`
stiff and on a ±45° I-section gives a GJ 30 % higher. Until 2026-10 the default was
`"restrained"`, after the anisotropic-section-model comparison where a finite shell
segment of the IEA 15 MW blade, whose end constraints hold the section's shape, read
GJ 16 % above `"free"`. Choose `"restrained"` for a section held in shape by
diaphragms, ribs or close end constraints. GBT cross-section modes carry the hoop
bending themselves, so sections used by GBT elements are analysed restrained.

## Discretisation

- Cubic Lagrange elements for \(w\) along each wall, shared at section points, with 5-point Gauss integration.
- On straight walls the exact warping is at most cubic for the forces used here, so **one element per wall is exact**. The tests check that refining (`n_sub`) changes nothing.
- Curved walls are represented by more straight walls.

## Outputs

| Attribute | Meaning |
|---|---|
| `stiffness` | 6 × 6 Timoshenko stiffness \(\mathbf F = \mathbf K\boldsymbol\psi\) |
| `flexibility` | its inverse, built directly from complementary energy |
| `classical_stiffness` | 4 × 4 for \([\varepsilon, \kappa_y, \kappa_z, \phi]\) (shear-rigid: Euler-Bernoulli with Saint-Venant torsion) |
| `shear_centre` | point where a transverse force causes no twist |
| `tension_centre` | point where an axial force causes no bending |
| `warping_constant` | \(E\Gamma\): primary and secondary (through-thickness) warping, warping about the shear centre and orthogonal to extension and bending |
| `point_warping()` | torsional warping \(\omega\) at each section point |
| `mass_properties` | mass per length, mass centre, mass moments of inertia (through-thickness distribution exact) |
| `wall_strains(F)` | wall strains at the Gauss points for section forces \(\mathbf F\) |

`stiffness_about(K, y0, z0)` moves a stiffness to another reference point. `section_analyses_from_files` builds one analysis per section from `laminate.txt` and `section_geometry.txt`.

## Verification

All checks are in `tests/test_thin_walled_section_analysis.py`.

| Check | Reference | Agreement |
|---|---|---|
| Polygonal tube EA, EI, GJ | exact thin-wall line integrals, Bredt plus \(GLt^3/3\) | 1e-10 |
| Tube shear stiffness | \(G\pi R t\) (coefficient 1/2) | 2e-3 (polygon vs circle) |
| Box with two thicknesses: GJ, EI | Bredt with \(\oint ds/t\), line integrals | 1e-10 |
| I-beam GJ, \(E\Gamma\) | \(\sum Gbt^3/3\); \(t_fb^3h^2/24 + t_f^3b^3/72 + t_w^3h^3/144\) | 1e-10 |
| I-beam shear stiffness \(V_z\) | complementary energy of the classical shear flow | 1e-4 (wall-thickness terms) |
| Channel shear centre and centroid | \(3b^2/(h+6b)\), \(b^2/(h+2b)\) | 1e-6 |
| CUS tube extension-twist | \(EA = P A_{11}^*\), \(K_{\varepsilon\phi} = 2A_e A_{16}^*\), \(GJ = 4A_e^2A_{66}^*/P\) | 1e-10 |
| Translation and rotation of the geometry | rigid transformation of the strains | 1e-7 |
| Mirror-symmetric section and layup | no couplings | 1e-10 |
| Wall resultants | \(\int \mathbf Z^\top\mathbf K\mathbf e_w\,ds = \mathbf F\) | 1e-8 |

## Limits

- Thin walls only: walls are Kirchhoff plates, and the wall transverse shear stiffness \(\mathbf H\) is not used here. A wall's own transverse shear would add boundary layers about \(t/3\) wide at free edges and junctions, which would need graded meshes. Thick walls and through-thickness stresses are left to the 2D section solver (PRs 14–16).
- Straight walls. Curvature within a wall is represented by splitting it.
- Saint-Venant level: free warping. Restrained warping (Benscoter, Vlasov) and section distortion (GBT) come in later PRs and build on `torsion_warping`.
