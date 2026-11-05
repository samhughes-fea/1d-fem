# Cross-section modes and the GBT element

`CrossSectionModes` (`anisotropic_section_model/thin_walled/cross_section_modes.py`)
finds the conventional (Vlasov) deformation modes and the warping modes of any open,
closed or multi-cell thin-walled section. `LinearGBTBeamElement3D` is a beam with
plane-section kinematics, the warping modes (Saint-Venant warping plus optional
shear-lag modes) and the lowest `m` distortional modes.

## Modes from the section frame

The walls are a plane frame, split into the `n_sub` elements of the
`ThinWalledSectionAnalysis`. Each node has an in-plane displacement \((d_y, d_z)\),
a rotation \(\beta\) and a warping \(u\). On an element with tangent \(t\) and normal
\(n\): \(v\) (tangential) and \(u\) are linear, \(w\) (normal) is a cubic Hermite with
\(\mathrm dw/\mathrm ds = \beta\). The Vlasov hypotheses are linear constraints:

- no hoop membrane strain: \(v_{,s} = 0\);
- no membrane shear: \(u_{,s} + v = 0\).

Their null space is the admissible space. It contains the rigid-section modes: axial,
two translations with plane-section warping and, for open sections only, the rotation
with its sectorial warping. Those are carried by the beam DOFs. The distortional modes
are the admissible vectors \(C\)-orthogonal to them, from

\[
B\,x = \lambda\,C\,x,\qquad
C = \int \mathbf e_2^{\mathsf T}\mathbf K\,\mathbf e_2\,\mathrm ds,\quad
B = \int \mathbf e_0^{\mathsf T}\mathbf K\,\mathbf e_0\,\mathrm ds,
\]

lowest \(\lambda\) first. \(\mathbf K\) is the 5 × 5 GBT wall law
(`gbt_wall_stiffness`, \(N_{ss} = 0\)). Each mode is scaled to a unit largest in-plane
displacement. With `n_sub = 1` the modes are the "natural node" modes (for a box, the
distortion mode first). With `n_sub > 1` local plate-bending modes follow.

## Wall strains of a mode

A mode with amplitude \(a(x)\) displaces the walls by \(u_x = u(s)\,a'\) and
\((u_y, u_z) = (v t + w n)\,a\). In the order \([\varepsilon_{xx}, \gamma_{xs},
\kappa_{xx}, \kappa_{ss}, \kappa_{xs}]\) its wall strains are

| amplitude term | wall strains |
|---|---|
| \(a''\) | \([u, 0, -w, 0, 0]\) |
| \(a'\) | \([0, u_{,s} + v\ (=0), 0, 0, -2w_{,s}]\) |
| \(a\) | \([0, 0, 0, -w_{,ss}, 0]\) |

## Warping modes

In the element all out-of-plane warping is carried by warping modes \(\Omega_j(s)\) with
their own amplitudes \(b_j(x)\): \(u_x \mathrel{+}= \Omega_j b_j\), so
\(\varepsilon_{xx} \mathrel{+}= \Omega_j b_j'\) and \(\gamma_{xs} \mathrel{+}= \Omega_{j,s} b_j\).
They use the warping mesh of the `ThinWalledSectionAnalysis` (`n_sub` cubic elements per wall).

- The first \(n_{sv}\) modes span the Saint-Venant warping \(W_0, W_1\) of the six beam strains,
  less its rigid part (\(1, y, z\)). They are always carried. For a rectangular box
  \(n_{sv} = 3\) (shear warping for \(F_y\), \(F_z\), and torsional warping); for a square
  box the torsional warping vanishes and \(n_{sv} = 2\).
- Further shear-lag modes come from \(E\,x = \mu\,A\,x\) in the complement, with \(A\) the
  axial warping stiffness and \(E\) the membrane shear stiffness. `[shear_modes]` adds the
  first few.
- All modes are \(A\)-orthogonal to \(1, y, z\) and to each other.

Why plane-section kinematics: if the warping modes were added on top of the Saint-Venant
warping, a uniform or linearly varying force state would not couple to them, and
\(b = 0\) at a support would not restrain the Saint-Venant warping. Carrying all warping
in \(b\) means a support with \(b = 0\) keeps the section plane, which gives shear lag and
the warping-restraint end effect. With \(b\) free and tip loads, minimising over \(b\)
gives back the Saint-Venant section. The small term of strains from \(\psi'\) in the
Saint-Venant solution is left out, which is exact for tip loads.

## Section stiffness

`CrossSectionModes.stiffness(m, p_extra)` is \(\int \mathbf H^{\mathsf T}\mathbf K\,\mathbf H\,\mathrm ds\)
over the strains \([\boldsymbol\psi_6, b'_{1..p}, b_{1..p}, a''_{1..m}, a'_{1..m}, a_{1..m}]\),
\(p = n_{sv} + p_{extra}\):

\[
\begin{aligned}
\varepsilon_{xx} &= Z_0\,\boldsymbol\psi_6 + \Omega\,b' + u_k a_k'' \\
\gamma_{xs} &= Z_1\,\boldsymbol\psi_6 + \Omega_{,s}\,b \\
\kappa_{xx} &= Z_2\,\boldsymbol\psi_6 - w_k a_k'' \\
\kappa_{ss} &= -w_{k,ss}\,a_k \\
\kappa_{xs} &= -2\phi - 2 w_{k,s}\,a_k'
\end{aligned}
\]

\(Z\) is the plane-section map of the walls. \(\kappa_{ss}\) is only carried by the modes,
which is why the section must use `hoop="restrained"`. The mode-mode blocks are the
classical GBT matrices \(C\), \(D\) and \(B\). The mode-beam and mode-warping blocks couple
distortion to torsion, bending, extension and warping; for anisotropic walls they come out
of the same integral.

`CrossSectionModes.mass(m, p_extra)` is the section mass for
\([u, v, w, \theta_x, \theta_y, \theta_z, b, a, a']\) from the wall displacement field.
`SectionAnalysis.from_thin_walled(..., gbt=True)` stores both, element-signed, and
`gbt_matrices(m, p_extra)` slices the first modes.

## Wall loads

A surface load \(\mathbf p = p_x\mathbf e_x + p_s\mathbf t + p_n\mathbf n\) on a wall does work
\(\int \mathbf p\cdot\mathbf u\,\mathrm ds\) per unit length. `CrossSectionModes.wall_load` gives
the generalised load on \([u, v, w, \theta_x, \theta_y, \theta_z, b, a, a']\), so an eccentric
flange load drives distortion as well as bending and torsion. Elements without modes use
`ThinWalledSectionAnalysis.wall_load_resultant`: \([N, F_y, F_z, M_x, M_y, M_z]\) and the
Benscoter \(\psi\) load \(\int \omega\,p_x\,\mathrm ds\). Loads come from `wall_load.txt`.

## Element

- `LinearGBTBeamElement3D` (a Timoshenko subclass). Nodal extras are
  \([b_1, \dots, b_p, a_1, a_1', \dots, a_m, a_m']\).
- \(b\) is linear in \(x\); the mode amplitudes are cubic Hermite.
- Integration: a full rule of at least 4 points. The shear block
  \(\{\gamma_{xy}, \gamma_{xz}, b\}\) is an assumed strain: its projection on polynomials
  of degree shear-block order − 1 (default order 1, the centre value), so the element does
  not lock. The energy is \(\int \bar B^T D \bar B\,dx\) on the full rule, which stays
  positive when the section couples the shear block with the cubic mode amplitudes. (Taking
  \(D\)'s shear block on one point and the rest on the full rule did not: on the IEA blade's
  2D sections it gave a spurious 5.2 Hz distortional mode.) Stress recovery and the
  corotational element use the same \(\bar B\).
- Input: element type `LinearGBTBeamElement3D`, and optional `[modes]` and `[shear_modes]`
  columns in `element_section.txt`.
- Supports: warping restraint is \(b_j = 0\); a diaphragm or rib is \(a_k = 0\); a fully
  fixed end has all extras zero.

## Validation

`tests/test_gbt_cross_section_modes.py`:

- **Square box distortion mode.** The corners move along the diagonals and the joints do not
  rotate. \(B = 192 D/a^3\) (four fixed-fixed walls, end offset 2) and
  \(C = E t a^3/3 + 68 D a/35\), both from hand frame analysis.
- **Open I-section.** Four rigid modes; the modes are \(B\)- and \(C\)-orthogonal and ordered.
- **Single-mode cantilever.** Checked against the exact solution of
  \(C a'''' - (D - 2K_{20})a'' + B a = 0\) (proof:
  [../proofs/gbt/single_mode_cantilever.md](../proofs/gbt/single_mode_cantilever.md)),
  converging as \(O(h^4)\).
- **Free warping is Saint-Venant.** With \(b\) free at the root, tip forces and torque give the
  Timoshenko beam with the Saint-Venant section to \(2\times10^{-5}\).
- **Restrained torsion.** An I-beam with \(b = 0\) at the root matches the Benscoter element to
  \(10^{-8}\), with or without extra shear-lag modes.
- **Shear lag.** A wide box (0.8 × 0.15 m, 1.5 m span) with a plane root section is about 1.7 %
  stiffer than with free warping. Each added set of shear-lag modes softens it, and it stays
  below the free case.
- **GBT against Benscoter.** With no distortional modes, a box cantilever (0.3 × 0.12 m)
  equals the Benscoter cantilever in torsion to round-off. In flap bending it is 0.11 %
  stiffer, because the plane root section also restrains the shear warping. Holding the
  flap shear-warping amplitude at every node makes the box 2.18 % stiffer in flap with the
  same tip rotation. The section then cannot warp under transverse shear, so its shear
  stiffness is the plane-section value, not the Saint-Venant one. This is the 2.2 %
  difference first seen when the curved GBT element was checked.
- **Element consistency.** With the modes fixed the element is the warping-only element, and
  freeing them only softens it.
- **End to end.** `process_job` runs with `[modes]`, `[shear_modes]`, `[Section_Mesh]`,
  `wall_load.txt` and the `Q` supports.

`tests/test_gbt_shear_modes_and_wall_loads.py` checks the warping modes (span and
orthogonality), the wall-load resultants by hand, and the parsers.

## GBT on 2D sections

`SolidCrossSectionModes` (`anisotropic_section_model/solid/cross_section_modes.py`)
gives the same generalised strains, stiffness and mass for a section meshed in 2D, so
`LinearGBTBeamElement3D` takes either kind of section. Each section point moves as

\[
u_x = u - y\theta_z + z\theta_y + \Omega_j b_j + u_k a_k', \qquad
(u_y, u_z) = (v - z\theta_x,\; w + y\theta_x) + S\,\psi + d_k a_k .
\]

- \(\Omega_j\): the axial part of the Saint-Venant warping of the six beam strains and their
  rates, less \(1, y, z\), made orthonormal.
- \(S\psi\): the in-plane (Poisson) Saint-Venant warping per unit beam strain, slaved to the
  beam strains (without it EA is held in plane strain).
- \(d_k\): the lowest in-plane eigenmodes of the section's plane-strain stiffness against its
  area, rigid in-plane motions removed. They include transverse wall extension and local
  plate modes of each wall, which the thin-walled basis cannot show.
- \(u_k\): the axial warping that minimises the out-of-plane shear energy of \(a_k'\), with
  zero mean, cleared of the rigid-motion parts (the translations paired with \(-y, -z\) and the
  rotation paired with the Saint-Venant torsion warping) so that \(a_k\) does not duplicate
  bending or torsion.
- The axial strain of \(a_k''\), \(\partial_x u_k\), comes with a slaved strain across the wall
  that makes \(\sigma_{nn} = 0\) (\(n\) normal to each element's wall direction
  \(\alpha\)), as in plate theory: the wall is free to thin. Without it the longitudinal
  bending of a wall is held in plane strain and `C` is
  \((1-\nu)/((1+\nu)(1-2\nu))\) times the plate value instead of \(1/(1-\nu^2)\), 22 % high
  for \(\nu = 0.3\), which put 2D local buckling 10 % high. The membrane part of
  \(\partial_x u_k\) keeps its hoop strain restrained, whereas the thin-walled modes release
  it (\(N_{ss} = 0\)), so for a membrane-dominated mode `C` can be up to \(1/(1-\nu^2)\)
  above the thin-walled value (none for a ply whose \(\nu_{21}\) is small).

The stiffness is \(\int H^T C H\,\mathrm dA\) with the 3D material law, so anisotropy and
Poisson coupling between distortion and the beam strains are kept. Extra shear-lag modes are
not offered for 2D sections.

`tests/test_gbt_2d_sections.py`, on a 0.3 × 0.12 m aluminium box with 2 mm walls (1.5 m
cantilever, 16 elements). The 2D mesh comes from `SectionMesh.from_thin_walled`, which
mitres the corners of a closed section like this one (see below):

- Statics: with 0 or 6 modes, the 2D section gives the thin-walled GBT tip response within
  3 %.
- Free vibration against a CalculiX `C3D20R` model of the same mitred box
  ([gbt_2d_box](../../validation/gbt_2d_box/make_reference.py)).
  CalculiX gives flap bending at 59.1 Hz, torsion with cross-section distortion at 76.9 Hz,
  local wall modes from 83.9 Hz, and edge bending at 128.9 Hz.

| Lowest eight modes (Hz) | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| CalculiX `C3D20R` | 59.1 | 76.9 | 83.9 | 87.5 | 94.7 | 105.8 | 106.0 | 113.9 |
| 2D GBT, 4 modes | 59.7 | 83.3 | 85.8 | 89.3 | 96.6 | 108.2 | 110.5 | 118.8 |
| Thin-walled GBT, 4 modes | 58.2 | 76.4 | 82.2 | 85.6 | 92.6 | 100.8 | 103.6 | 108.5 |

- Both GBT models find the local wall modes within 5 %. 2D GBT puts the torsion-distortion
  mode 8 % high; thin-walled GBT gets it within 1 %.
- Without distortional modes, edge bending (130.1 Hz) is within 1 % of CalculiX. Flap
  bending is then 11 % high (65.5 Hz), because the wide flanges' shear lag needs the modes.
- With more than about 8 modes the 2D GBT eigenproblem becomes ill-conditioned and gives
  spurious frequencies. Use up to 6.

**Junctions.** `SectionMesh.from_thin_walled(..., junctions=...)` joins walls in one of
three ways.

- **Mitre.** At a point where two walls meet, one ending and the next starting there, with
  the same ply thicknesses, both strips end on the line that bisects the corner. Each ply
  row meets the same row of the other wall there, and the two share nodes: the mesh is
  conforming.
- **Embedded tie.** At a T-junction, the joining wall's end edge lies inside the through
  wall's strip, and each of its nodes moves with that strip at its own position.
- **Butt joint.** At a two-wall corner that cannot be mitred (different thicknesses, or
  walls running against each other), the first wall runs on past the corner by half the
  other wall's thickness. The other wall's end edge is then tied inside it, as at a T.

`"auto"`, the default, uses all three, so no tied node is ever pinned. `"mitre"` refuses a
section with a junction it cannot mitre. `"tie"` is the old meshing: strips end on their
mid-line points, and the part of a corner's end edge that lies outside the master strip is
pinned to the strip's end, so the edge cannot rotate with the master wall. That stiffens
the corner:

- tied distortional modes come out too stiff (2D GBT found no distortional mode below
  124 Hz on the tied box);
- tied 3D extrusions are too stiff, as the table below shows.

| 3 m `C3D20R` prism, tip rotation / free-warping beam | `"tie"` | `"auto"` |
|---|---|---|
| 0.1 m CAS CFRP box: torsion, flap bending | 0.55, 0.65 | 1.01, 0.99 |
| Two-cell box with a middle web: torsion, flap, edge | 0.92, 0.77, 0.89 | 1.00, 0.97, 0.99 |
| L of 8 mm and 4 mm walls: torsion, bending y, bending z | 0.04, 0.90, 0.58 | 0.95, 0.99, 1.00 |

The rest of the `"auto"` difference is end restraint: the clamped root and the rigid tip
hold the warping that the free-warping beam lets go.

An earlier version of this check compared tied boxes. The 1.5 % agreement it first reported
came from boxes that had lost their ties, and the 2.4 % and 4.2 % that replaced it came from
tied boxes. Both are withdrawn.

## Tapered and pretwisted beams

Neighbouring GBT elements share the warping amplitudes \(b_j\) and the mode amplitudes
\(a_k\), so on a beam with its own section per element the \(j\)-th field of one section must
continue the \(j\)-th field of the next. Each section builds its bases on its own (an
eigen-decomposition each), and between two close sections a basis can come out turned within
its span, two modes can swap, or a mode can change sign.

`gbt_alignment.align_gbt_sections` walks the beam from its first element and turns each
section's warping block and mode block to the neighbour's already aligned, by the orthogonal
\(R = \operatorname{polar}(S_{\mathrm{ref}} S^T)\) (Procrustes, regularised towards the
identity). A row of \(S\) holds one field's shape moments over its section
(`SectionAnalysis.gbt_moments`: its displacements against the monomials up to degree 3 in the
section's bounding-box coordinates), which vary smoothly along a smoothly varying beam. A
turn within a block changes the basis, not the span, so every GBT matrix becomes
\(T^T K T\) with \(T\) the block-diagonal of \(R^T\) and each section's model is unchanged.
A match below 0.9 (the lowest cosine between turned and reference rows) is logged: refine the
mesh or change the number of modes. Thin-walled and 2D section jobs align their GBT sections
automatically.

The Resal taper coupling and the pretwist coupling (`with_taper`, `with_pretwist`) act on the
beam strains only, so they are added to the leading 6×6 block of `gbt_stiffness`, as for the
Benscoter section. On the IEA 15 MW outboard segment (20 elements, 2D sections, tapered and
pretwisted), GBT with the modes held gives the Timoshenko 2D tip deflections within 0.28 % in
flap and 0.09 % in edge; in torsion it is 3.7 % stiffer, from the root's warping restraint.

## Corotational GBT

`CorotationalGBTBeamElement3D` puts the GBT element in the corotational frame of
[corotational_section_beams.md](corotational_section_beams.md). The mode amplitudes pass
through as local DOFs. The second-order terms come from the section forces of the beam
strains, so the distortional modes carry no wall-level geometric stiffness. Checked at element
level (tangent, objectivity) and as a job: at small load it is the linear GBT element, and at
a tip deflection of 14 % of the span it stiffens as the corotational Timoshenko beam does.

### Brazier flattening

A mode moves the walls in the section plane by \(\mathbf d_k(s)\), so the beam curvature
\(\boldsymbol\kappa = (\kappa_y, \kappa_z)\) strains the moved fibre axially by

\[
e_B = \sum_k a_k\, \mathbf g_k \cdot \boldsymbol\kappa, \qquad \mathbf g_k = (d_{z,k},\ -d_{y,k}).
\]

Its work against the wall force of the GBT strains and its own energy are added to the
element's local energy at each Gauss point:

\[
\Pi_B = \sum_k a_k\, \boldsymbol\kappa^T \mathbf T_k \mathbf E
  + \tfrac12 \sum_{kl} a_k a_l\, \boldsymbol\kappa^T \mathbf Q_{kl} \boldsymbol\kappa,
\qquad
\mathbf T_k = \int \mathbf g_k N_{xx}\, ds,\quad \mathbf Q_{kl} = \int A_{11}\, \mathbf g_k \mathbf g_l^T ds,
\]

where \(N_{xx}\) is the wall force per unit GBT strain \(\mathbf E\) (on 2D sections, the
integrals are over the mesh with \(\sigma_{xx}\) and the axial modulus). The bent section's
compressed and stretched flanges are pushed towards the neutral axis, the modes take that
up, and the bending stiffness falls (Brazier 1927). The terms are of third and fourth order
in the local DOFs, like the other corotational terms. `SectionAnalysis.gbt_brazier_terms`
builds them (`CrossSectionModes.brazier_tensor`, `SolidCrossSectionModes.brazier_tensor`).
Linear GBT has none.

| Check | Reference | Agreement |
|---|---|---|
| Energy gradient and Hessian | central differences | 1e-7 |
| Isotropic tube (24-sided, R/t = 50) under a tip moment | Brazier's ovalisation \(\zeta = (1-\nu^2)\kappa^2R^4/t^2\) | 1 % at \(\zeta\) = 1.2 %, 4 % at 5.5 % |
| Same tube | bending stiffness falls by \(1.5\zeta\) | within 10 % of the drop |
| CAS composite box, flattening at 1 P (19 % span deflection) | CalculiX `C3D20R` `NLGEOM` solid | 86 % of the solid's |
| CAS box, tip deflection up to 2 P (35 % span) | same solid | within 0.3 % (rigid section: 1.3 %) |
| CAS box, tip deflection at 3 P (49 % span) | same solid | −3.5 % (rigid section: −5.3 %) |

The composite box is
[nonlinear_solid_box](../../validation/nonlinear_solid_box/README.md).
The beam's flattening agrees with a frame estimate of the flanges (rigid corners, Brazier
pressure \(N_{xx}\kappa\)): 0.28 mm at a quarter span at P with four modes, 0.34 mm with
eight, against 0.39 mm in the solid. Between 2 P and 3 P the solid's flattening grows faster
than the square of the curvature, to 3.3 mm, while the beam's reaches 2.2 mm. So at 3 P the
beam is still 3.5 % stiff, against a 2 % target. Modes beyond eight and finer wall meshes
change it by less than 0.1 %. Adding the modes' wall geometric stiffness changes it by less
than 0.1 %, and on a coarse mesh it loses convergence near 2.75 P. It is left out.

## Buckling

The element's geometric stiffness takes the wall stresses of the GBT strains, so the
distortional and local modes buckle: see [gbt_buckling.md](gbt_buckling.md).

## Not yet included

Transverse-extension modes on thin-walled sections (in-plane wall stretching,
\(N_{ss} \neq 0\)); use a 2D section, whose modes include them.
