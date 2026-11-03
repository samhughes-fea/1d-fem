# Cross-section modes and the GBT element

`CrossSectionModes` (`pre_processing/section_integration/thin_walled/cross_section_modes.py`)
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
  \(\{\gamma_{xy}, \gamma_{xz}, b\}\) uses the shear-block order (default 1), so the element
  does not lock.
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
- **Element consistency.** With the modes fixed the element is the warping-only element, and
  freeing them only softens it.
- **End to end.** `process_job` runs with `[modes]`, `[shear_modes]`, `[Section_Mesh]`,
  `wall_load.txt` and the `Q` supports.

`tests/test_gbt_shear_modes_and_wall_loads.py` checks the warping modes (span and
orthogonality), the wall-load resultants by hand, and the parsers.

## GBT on 2D sections

`SolidCrossSectionModes` (`pre_processing/section_integration/solid/cross_section_modes.py`)
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

The stiffness is \(\int H^T C H\,\mathrm dA\) with the 3D material law, so anisotropy and
Poisson coupling between distortion and the beam strains are kept. Extra shear-lag modes are
not offered for 2D sections.

`tests/test_gbt_2d_sections.py`, on a 0.3 × 0.12 m aluminium box with 2 mm walls (1.5 m
cantilever, 16 elements):

- Statics: the 2D section with 0 or 6 modes gives the thin-walled GBT tip response within
  0.8 % (largest on the torsion-coupled lateral rotation).
- Free vibration against a CalculiX `C3D20R` model of the same box
  ([gbt_2d_box](../../post_processing/validation_visualisers/gbt_2d_box/make_reference.py)):
  the lowest eight 2D GBT frequencies (6 modes) are each within 1.5 % of a CalculiX mode,
  the first flange mode within 0.4 % (44.50 vs 44.33 Hz). The thin-walled GBT with 6 modes
  starts at 58.1 Hz and misses the flange modes.

## Corotational GBT

`CorotationalGBTBeamElement3D` puts the GBT element in the corotational frame of
[corotational_section_beams.md](corotational_section_beams.md). The mode amplitudes pass
through as local DOFs. The second-order terms come from the section forces of the beam
strains, so the distortional modes carry no wall-level geometric stiffness. Checked at element
level (tangent, objectivity) and as a job: at small load it is the linear GBT element, and at
a tip deflection of 14 % of the span it stiffens as the corotational Timoshenko beam does.

## Not yet included

Transverse-extension modes on thin-walled sections (in-plane wall stretching,
\(N_{ss} \neq 0\)); use a 2D section, whose modes include them. Wall-level geometric
stiffness of the distortional modes (local buckling).
