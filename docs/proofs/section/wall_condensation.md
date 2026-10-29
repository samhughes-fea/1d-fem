# Condensed wall stiffness

**Claim.** Let \(\mathbf{K}\) be a symmetric positive definite wall stiffness, split into kept (\(k\)) and released (\(r\)) components. If the released forces are zero, the kept forces and strains are related by

\[
\mathbf{f}_k = \mathbf{K}_c\,\mathbf{e}_k,\qquad
\mathbf{K}_c = \mathbf{K}_{kk} - \mathbf{K}_{kr}\mathbf{K}_{rr}^{-1}\mathbf{K}_{rk},
\]

and \(\mathbf{K}_c\) equals \((\mathbf{C}_{kk})^{-1}\), the inverse of the kept block of the compliance \(\mathbf{C} = \mathbf{K}^{-1}\). \(\mathbf{K}_c\) is symmetric positive definite and never stiffer than \(\mathbf{K}_{kk}\).

## Elimination

Write the wall law in blocks:

\[
\begin{bmatrix} \mathbf{f}_k \\ \mathbf{0} \end{bmatrix}
= \begin{bmatrix} \mathbf{K}_{kk} & \mathbf{K}_{kr} \\ \mathbf{K}_{rk} & \mathbf{K}_{rr} \end{bmatrix}
\begin{bmatrix} \mathbf{e}_k \\ \mathbf{e}_r \end{bmatrix}.
\]

\(\mathbf{K}_{rr}\) is a principal sub-matrix of a positive definite matrix, so it is positive definite and invertible. The second row gives \(\mathbf{e}_r = -\mathbf{K}_{rr}^{-1}\mathbf{K}_{rk}\mathbf{e}_k\). Putting this into the first row gives \(\mathbf{f}_k = \mathbf{K}_c\mathbf{e}_k\).

## Equal to the inverse of the kept compliance

The block-inverse formula for a symmetric matrix gives the kept block of \(\mathbf{K}^{-1}\) as the inverse of the Schur complement of \(\mathbf{K}_{rr}\):

\[
(\mathbf{K}^{-1})_{kk} = \left(\mathbf{K}_{kk} - \mathbf{K}_{kr}\mathbf{K}_{rr}^{-1}\mathbf{K}_{rk}\right)^{-1} = \mathbf{K}_c^{-1}.
\]

So \(\mathbf{K}_c = (\mathbf{C}_{kk})^{-1}\). Physically: with the released forces zero, \(\mathbf{e}_k = \mathbf{C}_{kk}\mathbf{f}_k\), which is the compliance statement of the same law.

## Positive definite, and never stiffer

\(\mathbf{C}\) is positive definite, so \(\mathbf{C}_{kk}\) is too, and so is its inverse \(\mathbf{K}_c\). Symmetry follows from the symmetry of \(\mathbf{K}\).

\(\mathbf{K}_{kk} - \mathbf{K}_c = \mathbf{K}_{kr}\mathbf{K}_{rr}^{-1}\mathbf{K}_{rk}\) is positive semi-definite, because \(\mathbf{x}^\top\mathbf{K}_{kr}\mathbf{K}_{rr}^{-1}\mathbf{K}_{rk}\mathbf{x} = \mathbf{y}^\top\mathbf{K}_{rr}^{-1}\mathbf{y} \ge 0\) with \(\mathbf{y} = \mathbf{K}_{rk}\mathbf{x}\). Releasing a resultant can only soften the kept law. It leaves it unchanged when the kept and released strains are uncoupled (\(\mathbf{K}_{kr}=0\)).

## Isotropic limit

For an isotropic wall of thickness \(t\), \(\mathbf{B}=0\) and

\[
\mathbf{A} = \frac{Et}{1-\nu^2}\begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \tfrac{1-\nu}{2} \end{bmatrix},\qquad
\mathbf{D} = \frac{t^2}{12}\mathbf{A}.
\]

Release \(N_{ss}\) and \(M_{ss}\). Membrane and bending are uncoupled, so each condenses alone. For the membrane block, \(A_{11} - A_{12}^2/A_{22} = \frac{Et}{1-\nu^2}(1-\nu^2) = Et\), \(A_{16}=A_{26}=0\), and the shear term stays at \(A_{66} = Et/(2(1+\nu)) = Gt\). The bending block gives \(Et^3/12\) and \(Gt^3/12\) in the same way. So

\[
\mathbf{K}_c = \mathrm{diag}\left(Et,\ Gt,\ \frac{Et^3}{12},\ \frac{Gt^3}{12}\right)
\]

in the order \([\varepsilon_{xx}, \gamma_{xs}, \kappa_{xx}, \kappa_{xs}]\). These are the wall terms whose section integrals give \(EA\), \(EI\), and the open-wall \(GJ = \sum G\,L\,t^3/3\) (the factor 4 comes from the twist curvature \(\kappa_{xs} = 2\phi'\) in the later section step).

Tests: `tests/test_thin_walled_laminate.py` (`test_isotropic_beam_wall_stiffness`, `test_condensed_equals_inverse_of_kept_compliance`, `test_condensation_never_stiffens`).
