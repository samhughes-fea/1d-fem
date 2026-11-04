# Single GBT mode on a cantilever

With the beam DOFs held, one mode with amplitude \(a(x)\) has the energy

\[
\Pi = \tfrac12\int_0^L \big(C a''^2 + D a'^2 + B a^2 + 2K_{20}\,a''a\big)\,\mathrm dx - P\,a(L),
\]

where \(C, D, B\) are the diagonal GBT terms and \(K_{20}\) the \(a''\)–\(a\) coupling
(the plate Poisson term \(D_{12}\int w\,w_{,ss}\,\mathrm ds\)). For isotropic walls the
\(a''\)–\(a'\) and \(a'\)–\(a\) terms vanish.

**Stationarity.** Integrating by parts:

\[
C a'''' - (D - 2K_{20})\,a'' + B a = 0 \quad\text{on } (0, L),
\]

with \(a(0) = a'(0) = 0\) at the clamped end and, at \(x = L\),

\[
C a''(L) + K_{20}\,a(L) = 0, \qquad -C a'''(L) + (D - K_{20})\,a'(L) = P .
\]

**Solution.** \(a = \sum_{i=1}^4 c_i e^{r_i x}\), where \(r_i\) are the roots of
\(C r^4 - (D - 2K_{20}) r^2 + B = 0\), and the four boundary conditions give \(c_i\). The
test solves this 4 × 4 system directly and compares it with the element at several mesh
sizes. The Hermite element converges as \(O(h^4)\) at the tip.
