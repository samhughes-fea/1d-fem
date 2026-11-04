# Local, distortional and global buckling with GBT elements

`LinearGBTBeamElement3D.linear_geometric_stiffness_matrix(U_e)` gives the geometric
stiffness of the GBT element. With it, linear buckling (`[Type] buckling`) and P-delta
work with GBT elements as with the other beams ([second_order.md](second_order.md)). The
modes the section carries ([gbt_cross_section_modes.md](gbt_cross_section_modes.md)) then
buckle too, so one model covers three kinds of buckling:

- local plate buckling of the walls (modes that leave the folds in place);
- distortional buckling (modes that move the folds, such as a flange and lip rotating
  about the web junction);
- global flexural, torsional and flexural-torsional buckling (the beam DOFs).

The interaction between them is in the same eigenproblem, so the lowest mode can mix them.

## Formulation

The section moves in its own plane by

\[
\mathbf d(s, x) = \Pi(s)\,\mathbf q(x), \qquad \mathbf q = [v, w, \theta_x, a_1, \dots, a_m],
\]

where the rigid columns of \(\Pi\) are \((1, 0)\), \((0, 1)\) and \((-z, y)\), and mode
\(k\) adds its in-plane field \(\mathbf d_k(s) = v_k \mathbf t + w_k \mathbf n\). The
quadratic parts of the wall strains are

\[
\varepsilon_{xx}^{(2)} = \tfrac12 |\mathbf d_{,x}|^2, \qquad
\gamma_{xs}^{(2)} = \mathbf d_{,x} \cdot \mathbf d_{,s},
\]

and their work against the prebuckling membrane resultants \(N_{xx}\) and \(N_{xs}\) is, per
unit length,

\[
\tfrac12\,\mathbf q'^T X\,\mathbf q' + \mathbf q'^T Y\,\mathbf q, \qquad
X = \int N_{xx}\,\Pi^T\Pi\,ds, \qquad Y = \int N_{xs}\,\Pi^T\Pi_{,s}\,ds .
\]

- For the rigid rows these are the section terms of the beams (`N`, `S_y`, `S_z`, the
  Wagner term `W`, `V_y`, `V_z`, `Q_r`). The difference is that here they are integrals of
  the GBT element's own stresses, which include shear lag, warping restraint and the
  distortion of the prebuckling state.
- The mode rows are the classical GBT geometric matrices (\(X_{jk} = \int N_{xx}\,
  \mathbf d_j \cdot \mathbf d_k\,ds\) for uniform compression), plus their coupling with the
  rigid motion and the shear terms.
- The axial displacement and the through-thickness variation of the in-plane motion are
  left out, as in the beams' own second-order terms.

`SectionAnalysis.gbt_geometric_terms(m, p)` stores \(X\) and \(Y\) per unit GBT strain
(`CrossSectionModes.geometric_tensor` on thin-walled sections). At each Gauss point the
element evaluates them on its strains \(\bar B U_e\) (the assumed-strain B̄ of the stiffness)
and integrates

\[
K_\sigma = \int G^T \begin{bmatrix} X & Y \\ Y^T & 0 \end{bmatrix} G\,dx,
\qquad G\,U_e = [\mathbf q', \mathbf q],
\]

on the full rule. The slopes \(v'\) and \(w'\) are the linked slopes of the Timoshenko
family, and the mode amplitudes are the element's cubic Hermite fields.

On 2D sections (`SolidCrossSectionModes.geometric_tensor`) the same terms come from the 3D
stresses: \(X = \int \sigma_{xx}\,\Pi^T\Pi\,dA\) and
\(Y = \int (\tau_{xy}\,\Pi^T\Pi_{,y} + \tau_{xz}\,\Pi^T\Pi_{,z})\,dA\).

## Supports

A simply supported end (fork plus diaphragm) is `UY UZ RX` and `a_k = 0` (`Q(p+1)`,
`Q(p+3)`, …, with `p` the warping DOFs), with the warping `b` free. `Q*` fixes every extra
DOF, which also restrains the warping and the mode slopes.

## Mesh

A local mode has a half-wavelength of about the plate width, and the cubic mode amplitudes
need about six elements per half-wave for 1 %. A long column with local buckling therefore
needs short elements along all of its length, even though the global mode would not. A
coarse mesh can only raise the local and distortional loads, never invent a lower one.

`n_sub` sets the local modes. With `n_sub = 1` a wall has no internal node, so the section
has distortional modes only. With `n_sub = 2` each wall gets one local plate mode, which is
within 1 % of plate theory for the first local mode.

## Checks

- **Square tube, local buckling:** 300 × 300 × 2 mm, 1.2 m, `n_sub = 4`, 8 modes, 16
  elements: 27 001 N against the plate value \(k = 4\), 26 994 N (+0.03 %).
- **Lipped channel column:** 100 × 60 × 10 × 2 mm, `n_sub = 2`, 10 modes, against a finite
  strip model of the same column with full plate theory, lowest over the number of
  half-waves. From 40 mm to 3 m GBT is within 1.3 %: 0.3 to 1.3 % stiff in local and
  distortional buckling (the Vlasov hypotheses), within 0.1 % in global buckling.
- **Global:** with no modes, the Euler load of a square tube is the Timoshenko element's to
  2e-5.
- **2D sections:** on a 300 × 120 × 2 mm box meshed in 2D, the Euler loads are the 2D
  Timoshenko beam's to 4e-4. The geometric term of the first local mode
  (\(X_{kk}/N\)) equals the thin-walled section's within 0.4 %, and its warping stiffness
  `C` equals it to 0.03 %. The local buckling stress of a 0.3 m half-wave converges on the
  thin-walled GBT and finite strip value (15.70 MPa) from above as the 2D mesh is refined:

  | Elements along each wall, rows through it | 2D local buckling stress | Against thin-walled |
  |---|---|---|
  | 12, 1 (`box_mesh`) | 16.67 MPa | +6.2 % |
  | 24, 2 | 16.00 MPa | +1.9 % |
  | 48, 2 | 15.84 MPa | +0.9 % |

  What is left is the transverse bending `B` of the coarse mesh (+9 % at 12 elements along,
  +1.2 % at 48). It was 10 % high (17.3 MPa) before the modes' walls were free to thin:
  the axial strain of \(a_k''\) held the wall in plane strain, so `C` was
  \((1-\nu)/((1+\nu)(1-2\nu))\) times the plate value instead of \(1/(1-\nu^2)\), 22 % high
  (`gbt_cross_section_modes.md`).

## Example

[examples/thin_walled/gbt_buckling](../../examples/thin_walled/gbt_buckling/README.md)
computes the lowest buckling stress of the lipped channel from 40 mm to 3 m.

## Linear buckling solver

`solve_linear_buckling_eigenpairs` uses Lanczos (`eigsh` on \(-K_\sigma x = \mu K x\), with
\(\lambda = 1/\mu\) and `K` factorised once) when the free problem has more than 600 DOFs
and `K` is symmetric. A GBT column with ten modes
has 29 DOFs per node. It falls back to the dense solver for a nonsymmetric `K`
(Levinson elements) or when Lanczos does not return enough positive eigenvalues.

## Limits

- The section modes follow the Vlasov hypotheses (no membrane shear, no transverse
  extension), so GBT is slightly stiffer than the finite strip model, which has neither
  restriction.
- Transverse-extension modes are not offered on thin-walled sections.
- `CorotationalGBTBeamElement3D` keeps the beam-only geometric stiffness of the section
  forces, so its distortional modes still do not buckle.
