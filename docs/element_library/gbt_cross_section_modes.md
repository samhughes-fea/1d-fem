# Cross-section modes and the GBT element

`CrossSectionModes` (`pre_processing/section_integration/thin_walled/cross_section_modes.py`)
finds the conventional (Vlasov) deformation modes of any open, closed or multi-cell
thin-walled section, and `LinearGBTBeamElement3D` adds the lowest `m` of them to the
Benscoter beam.

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

## Wall strains and section stiffness

A mode with amplitude \(a(x)\) displaces the walls by \(u_x = u(s)\,a'\) and
\((u_y, u_z) = (v t + w n)\,a\). In the order \([\varepsilon_{xx}, \gamma_{xs},
\kappa_{xx}, \kappa_{ss}, \kappa_{xs}]\) its wall strains are

| amplitude term | wall strains |
|---|---|
| \(a''\) | \([u, 0, -w, 0, 0]\) |
| \(a'\) | \([0, u_{,s} + v\ (=0), 0, 0, -2w_{,s}]\) |
| \(a\) | \([0, 0, 0, -w_{,ss}, 0]\) |

The section stiffness `CrossSectionModes.stiffness(m)` is
\(\int \mathbf H^{\mathsf T}\mathbf K\,\mathbf H\,\mathrm ds\) over the strains
\([\boldsymbol\varepsilon_8, a''_{1..m}, a'_{1..m}, a_{1..m}]\). Here \(\mathbf H\) stacks the
Benscoter wall strains (with \(\kappa_{ss} = 0\), which is why the section must use
`hoop="restrained"`) and the mode columns above. The leading 8 × 8 block is
`benscoter_stiffness`. The mode-mode blocks are the classical GBT matrices \(C\), \(D\)
and \(B\). The mode-beam blocks couple distortion to torsion, bending and extension; for
anisotropic walls they come out of the same integral.

`CrossSectionModes.mass(m)` is the section mass for
\([u, v, w, \theta_x, \theta_y, \theta_z, \psi, a, a']\) from the wall displacement field.
The element uses it with linear rigid fields and Hermite \(a\).

## Element

- `LinearGBTBeamElement3D`: Benscoter element plus \(m\) modes. Nodal extras are
  \([\psi, a_1, a_1', \dots, a_m, a_m']\).
- Amplitudes are cubic Hermite in \(x\).
- Integration: a full rule of at least 4 points, with the \(\{\gamma_{xy}, \gamma_{xz}, d\}\)
  block reduced as in the Benscoter element.
- Input: `[modes]` in `element_section.txt`, `[warping] = 1` in `element.txt`.
- Supports: a diaphragm or rib is \(a_k = 0\); a fully fixed end also has \(a_k' = 0\).

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
- **Element consistency.** With the modes fixed the element equals the Benscoter element, and
  freeing them only softens it.
- **End to end.** `process_job` runs with `[modes]` and the `Q` supports.

## Not yet included

Wall surface loads on the modes, shear (non-Vlasov) modes for closed-section shear
lag, and transverse-extension modes. These are the next step of the plan.
