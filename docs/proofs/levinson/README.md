# Reddy and Levinson beams: proofs

The derivation is in [element_library/higher_order_shear.md](../../element_library/higher_order_shear.md).
The statements below are checked in `tests/test_third_order_shear.py`.

## Reddy's cubic is the slaved Saint-Venant shear warping

For a rectangle of depth \(h\), \(z^3 = \tfrac{3h^2}{20}z + W_3\) with
\(W_3 \perp \{1, z\}\). Rewriting \(u_x = z\varphi - c_1z^3(\varphi + w')\) with
\(\bar\varphi = \varphi - 0.2(\varphi + w')\) gives \(u_x = z\bar\varphi + \Psi\bar\gamma\),
\(\bar\gamma = \bar\varphi + w'\), \(\Psi = -\tfrac{5}{3h^2}W_3\). The shear energy is
\(\tfrac{8}{15}GA(\varphi + w')^2 = \tfrac56 GA\bar\gamma^2\), and \(\Psi\) is the
Saint-Venant warping per unit shear strain of the rectangle. Hence
\(E\int\Psi^2 dA = EI/84\).

## Clamped cantilever, tip load

With \(H = EI/84\) and \(\chi = \gamma'\), the Euler-Lagrange equations of
\(\tfrac12\int(EI\theta'^2 + \kappa GA\gamma^2 + H\gamma'^2) - Pw(L)\) are
\(EI\theta'' = \kappa GA\gamma - H\gamma'' = P\). The boundary conditions are
\(\theta(0) = \gamma(0) = 0\), \(\theta'(L) = \gamma'(L) = 0\). So
\(\gamma = \tfrac{P}{\kappa GA}\big(1 - \cosh\lambda(L - x)/\cosh\lambda L\big)\) and

\[
w(L) = \frac{PL^3}{3EI} + \frac{P}{\kappa GA}\Big(L - \frac{\tanh\lambda L}{\lambda}\Big),
\qquad \lambda^2 = \frac{\kappa GA}{H}.
\]

## Free warping gives the Timoshenko beam

With the slopes free at the root the natural condition is \(H\gamma'(0) = 0\).
Under a tip load \(\gamma = P/\kappa GA\) is constant and solves everything, with
\(\chi = 0\), which is the Timoshenko solution. Under a uniform load the boundary
layers at the two ends integrate to zero, so the tip deflection is still
\(qL^4/8EI + qL^2/2\kappa GA\).

## Levinson equals Timoshenko for uncoupled sections

The Levinson weak form is \(\int B_6^{\mathsf T}(D_6 B_6 + D_{6\chi}B_\chi)\). It loses
its higher-order part when \(D_{6\chi} = 0\), which holds for any section whose slaved
warping carries no classical resultant (all isotropic sections, by the
deflation of \(\Psi\)).
