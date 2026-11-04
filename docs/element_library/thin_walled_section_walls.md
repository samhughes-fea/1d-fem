# Thin-walled section walls

This page covers the first step of anisotropic thin-walled section analysis: the stiffness of each wall. Later steps (Saint-Venant warping, shear flows, section stiffness, cross-section modes) build on it.

Code: `pre_processing/section_integration/thin_walled/`.

## Axes

A section lies in the \((y, z)\) plane and the beam axis is \(x\). Each wall is a straight segment from a start point \(a\) to an end point \(b\).

| Axis | Direction |
|------|-----------|
| \(x\) | beam axis |
| \(s\) | along the wall, \(\mathbf{e}_s = (\mathbf{b}-\mathbf{a})/L\) |
| \(n\) | through the thickness, \(\mathbf{e}_n = \mathbf{e}_x \times \mathbf{e}_s = (-e_{s,z},\ e_{s,y})\) |

Ply angle \(\theta\) is measured from \(x\) toward \(s\). Plies are listed bottom first, from \(n=-h/2\) to \(n=+h/2\), and the reference surface is the wall mid-surface. On a counter-clockwise closed cell, \(\mathbf{e}_n\) points into the cell.

## Ply stiffness

Plane-stress reduced stiffness in ply axes, with \(\nu_{21} = \nu_{12}E_2/E_1\):

\[
\mathbf{Q} = \frac{1}{1-\nu_{12}\nu_{21}}
\begin{bmatrix} E_1 & \nu_{12}E_2 & 0 \\ \nu_{12}E_2 & E_2 & 0 \\ 0 & 0 & (1-\nu_{12}\nu_{21})G_{12} \end{bmatrix}.
\]

With engineering strains \(\boldsymbol\varepsilon_{12} = \mathbf{T}_\varepsilon\,\boldsymbol\varepsilon_{xs}\),

\[
\mathbf{T}_\varepsilon = \begin{bmatrix} c^2 & s^2 & cs \\ s^2 & c^2 & -cs \\ -2cs & 2cs & c^2-s^2 \end{bmatrix},
\qquad c=\cos\theta,\ s=\sin\theta,
\]

equal strain energy gives \(\bar{\mathbf{Q}} = \mathbf{T}_\varepsilon^\top \mathbf{Q}\,\mathbf{T}_\varepsilon\). The tests check this against the explicit textbook formulas for \(\bar Q_{11}\) to \(\bar Q_{66}\).

Transverse shear rotates as a vector: \([\gamma_{13}, \gamma_{23}] = \mathbf{R}[\gamma_{xn}, \gamma_{sn}]\) with \(\mathbf{R} = [[c, s], [-s, c]]\), so

\[
\bar{\mathbf{Q}}_s = \begin{bmatrix} G_{13}c^2 + G_{23}s^2 & (G_{13}-G_{23})cs \\ (G_{13}-G_{23})cs & G_{23}c^2 + G_{13}s^2 \end{bmatrix}.
\]

## Wall laminate

For ply \(k\) between \(n_{k-1}\) and \(n_k\):

\[
(\mathbf{A},\mathbf{B},\mathbf{D}) = \sum_k \bar{\mathbf{Q}}_k \left(n_k - n_{k-1},\ \tfrac12(n_k^2-n_{k-1}^2),\ \tfrac13(n_k^3-n_{k-1}^3)\right),
\qquad
\mathbf{H} = k_s \sum_k \bar{\mathbf{Q}}_{s,k}(n_k-n_{k-1}),
\]

with \(k_s = 5/6\) by default. The wall law is

\[
\begin{bmatrix} \mathbf{N} \\ \mathbf{M} \end{bmatrix} = \begin{bmatrix} \mathbf{A} & \mathbf{B} \\ \mathbf{B} & \mathbf{D} \end{bmatrix}
\begin{bmatrix} \boldsymbol\varepsilon \\ \boldsymbol\kappa \end{bmatrix},
\qquad
\begin{bmatrix} Q_x \\ Q_s \end{bmatrix} = \mathbf{H}\begin{bmatrix} \gamma_{xn} \\ \gamma_{sn} \end{bmatrix},
\]

in the order \([N_{xx}, N_{ss}, N_{xs}, M_{xx}, M_{ss}, M_{xs}]\) against \([\varepsilon_{xx}, \varepsilon_{ss}, \gamma_{xs}, \kappa_{xx}, \kappa_{ss}, \kappa_{xs}]\).

## Condensed wall law

A beam theory does not model every wall resultant. A resultant it does not model is released (set to zero) and its strain is left free. Static condensation gives the stiffness of the kept strains:

\[
\mathbf{K}_c = \mathbf{K}_{kk} - \mathbf{K}_{kr}\mathbf{K}_{rr}^{-1}\mathbf{K}_{rk} = \left(\mathbf{C}_{kk}\right)^{-1},
\qquad \mathbf{C} = \mathbf{K}^{-1}.
\]

The proof is in [../proofs/section/wall_condensation.md](../proofs/section/wall_condensation.md).

| Function | Released | Kept strains | Size | Used by |
|----------|----------|--------------|------|---------|
| `beam_wall_stiffness` | \(N_{ss}, M_{ss}\) | \(\varepsilon_{xx}, \gamma_{xs}, \kappa_{xx}, \kappa_{xs}\) | 4 × 4 | thin-walled beam theory (rigid section in its own plane) |
| `gbt_wall_stiffness` | \(N_{ss}\) | \(\varepsilon_{xx}, \gamma_{xs}, \kappa_{xx}, \kappa_{ss}, \kappa_{xs}\) | 5 × 5 | cross-section distortion modes |
| `condensed_wall_stiffness` | any set | the rest | any | studies of the release assumption |

For an isotropic wall of thickness \(t\), `beam_wall_stiffness` is \(\mathrm{diag}(Et,\ Gt,\ Et^3/12,\ Gt^3/12)\). The plane-strain factor \(1/(1-\nu^2)\) drops out because \(N_{ss}\) and \(M_{ss}\) are free.

### Caveat for unsymmetric walls in closed cells

Releasing \(M_{ss}\) lets \(\kappa_{ss}\) take whatever value makes \(M_{ss}=0\). For an unsymmetric wall (\(\mathbf{B}\neq 0\)) that curvature is coupled to the membrane strains, so releasing it softens \(\varepsilon_{xx}\) and \(\gamma_{xs}\). In a closed cell the neighbouring walls partly restrain \(\kappa_{ss}\), so the fully released law can under-predict torsional stiffness. The anisotropic-section-model package measured GJ about 16% low on the IEA 15 MW blade with the fully released law. The release set is a parameter, so the section analysis in the next PR can compare both choices against 2D section results.
