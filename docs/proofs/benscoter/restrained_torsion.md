# Restrained torsion of a Benscoter beam

For an uncoupled section (\(K_8\) diagonal in the torsion-warping block), the energy per
length of the torsion fields is

\[
\tfrac12\big[GJ\,\theta'^2 + S\,(\psi - \theta')^2 + E\Gamma\,\psi'^2\big],
\]

with \(\theta = \theta_x\), \(GJ = K_{55}\), \(S = K_{66}\), \(E\Gamma = K_{77}\).

**Equilibrium.** Varying \(\theta\) with a tip torque \(T\):
\(GJ\theta' - S(\psi - \theta') = T\), so

\[
\theta' = \frac{T + S\psi}{GJ + S}.
\]

Varying \(\psi\): \(S(\psi - \theta') - E\Gamma\,\psi'' = 0\). Substituting \(\theta'\):

\[
E\Gamma\,\psi'' - \frac{S\,GJ}{GJ + S}\,\psi = -\frac{S\,T}{GJ + S},
\qquad k^2 = \frac{S\,GJ}{(GJ + S)\,E\Gamma}.
\]

The particular solution is \(\psi = T/GJ\), the Saint-Venant twist rate.

**Cantilever.** Clamped and warping-restrained at \(x = 0\) (\(\psi = 0\)), free at
\(x = L\) (bimoment \(E\Gamma\psi' = 0\)):

\[
\psi = \frac{T}{GJ}\big(1 - \cosh kx + \tanh kL\,\sinh kx\big),
\]

\[
\theta(x) = \frac{1}{GJ + S}\Big[T x + \frac{S\,T}{GJ}\Big(x - \frac{\sinh kx}{k}
 + \tanh kL\,\frac{\cosh kx - 1}{k}\Big)\Big].
\]

**Limits.**

- \(S \to \infty\) (Vlasov): \(k^2 \to GJ/E\Gamma\) and \(\psi \to \theta'\).
- \(E\Gamma \to 0\) or \(\psi\) free: \(\theta' = T/GJ\) everywhere.
- Closed sections have a large \(GJ\), so \(1/k\) is of the order of the section size
  and the restraint changes the tip twist very little. Open sections
  have \(1/k\) of the order of the span, and the restraint stiffens them strongly.

**Discretisation.** With linear \(\theta\) and \(\psi\), \(d = \psi - \theta'\) is linear
while \(\theta'\) is constant per element. Full integration of \(S d^2\) forces \(\psi\)
to be constant per element as \(S\) grows (locking). One-point integration of the
\(d\) block removes the constraint at the element level and keeps the energy exact for
the constant-\(d\) states, as for the Timoshenko shear terms
([../timoshenko/coupled_selective_integration.md](../timoshenko/coupled_selective_integration.md)).
