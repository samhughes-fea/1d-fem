# Derivation of Tip Deflection and Rotation for Euler–Bernoulli Cantilever Beams

**Summary.** This note derives the analytical **tip** (free-end) transverse displacement $u_y(L)$ and cross-section rotation $\theta_z(L)$ for a prismatic Euler–Bernoulli cantilever of length $L$ and flexural rigidity $EI$, for the load cases used in verification: point loads at the end, mid-span, and quarter-span; and uniform, triangular, and parabolic distributed loads. The sign convention is downward load and downward displacement positive, consistent with the finite-element output and with `DEFORMATION_AND_GCI_NOTES.md`.

---

## 1. Governing relations

For a prismatic Euler–Bernoulli beam with coordinate $x$ from the fixed end, curvature $\kappa = M/(EI)$, slope $\theta(x) = \theta_z(x)$, and transverse displacement $u(x) = u_y(x)$ satisfy
$$
\frac{\mathrm{d}\theta}{\mathrm{d}x} = \kappa = \frac{M}{EI}, \qquad \frac{\mathrm{d}u}{\mathrm{d}x} = \theta,
$$
with $u(0) = \theta(0) = 0$. Hence
$$
\theta(x) = \int_0^x \frac{M(\xi)}{EI}\,\mathrm{d}\xi, \qquad u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi.
$$
Shear and moment are related to the applied load $q(x)$ (positive downward) by $\mathrm{d}V/\mathrm{d}x = -q$ and $\mathrm{d}M/\mathrm{d}x = V$, with $V(x) = \int_x^L q(\xi)\,\mathrm{d}\xi$ and $M(x) = \int_x^L V(\xi)\,\mathrm{d}\xi$. All tip formulae below follow from evaluating $\theta(L)$ and $u(L)$ from the moment distribution $M(x)$ for each load case.

---

## 2. Point load at $x = a$

A concentrated force $P$ (downward positive) at $x = a$ gives
$$
M(x) = P(a - x) \quad \text{for } x \leq a, \qquad M(x) = 0 \quad \text{for } x > a.
$$
For $x \leq a$,
$$
\theta(x) = \frac{P}{2EI}\,x(2a - x), \qquad u(x) = \frac{P}{6EI}\,x^2(3a - x).
$$
At the load point, $\theta(a) = \frac{Pa^2}{2EI}$ and $u(a) = \frac{Pa^3}{3EI}$. For $x > a$, the moment is zero so $\theta(x) = \theta(a)$ and $u(x) = u(a) + \theta(a)(x - a)$.

### 2.1 Tip values (always at $x = L$)

- **End load** ($a = L$): the tip is at the load. Using the $x \leq a$ formulae with $x = L$ and $a = L$:
  $$
  \theta_z(L) = \frac{PL^2}{2EI}, \qquad u_y(L) = \frac{PL^3}{3EI}.
  $$

- **Mid-span load** ($a = L/2$): the tip is beyond the load ($L > a$). So $\theta_z(L) = \theta(a)$ and $u_y(L) = u(a) + \theta(a)(L - a)$:
  $$
  \theta_z(L) = \frac{P(L/2)^2}{2EI} = \frac{PL^2}{8EI}, \qquad u_y(L) = \frac{Pa^3}{3EI} + \frac{Pa^2}{2EI}(L - a)\Big|_{a = L/2} = \frac{PL^3}{24EI} + \frac{PL^3}{16EI} = \frac{5PL^3}{48EI}.
  $$

- **Quarter-span load** ($a = L/4$): again $L > a$, so
  $$
  \theta_z(L) = \frac{P(L/4)^2}{2EI} = \frac{PL^2}{32EI}, \qquad u_y(L) = \frac{P(L/4)^3}{3EI} + \frac{P(L/4)^2}{2EI}\,\frac{3L}{4} = \frac{PL^3}{192EI} + \frac{3PL^3}{128EI} = \frac{11PL^3}{384EI}.
  $$

---

## 3. Distributed loads

For distributed load intensity $q(x)$ (positive downward), $V(x) = \int_x^L q(\xi)\,\mathrm{d}\xi$ and $M(x) = \int_x^L V(\xi)\,\mathrm{d}\xi$. Then $\theta(x) = \int_0^x (M/EI)\,\mathrm{d}\xi$ and $u(x) = \int_0^x \theta\,\mathrm{d}\xi$. The **tip** rotation and deflection are $\theta_z(L) = \theta(L)$ and $u_y(L) = u(L)$.

### 3.1 Uniform (UDL): $q(x) = w$

$V(x) = w(L - x)$, $M(x) = \tfrac{1}{2}w(L - x)^2$. Integrating,
$$
\theta(x) = \frac{w}{6EI}\bigl(3L^2 x - 3L x^2 + x^3\bigr), \qquad u(x) = \frac{w}{24EI}\bigl(6L^2 x^2 - 4L x^3 + x^4\bigr).
$$
At the tip ($x = L$):
$$
\theta_z(L) = \frac{wL^3}{6EI}, \qquad u_y(L) = \frac{wL^4}{8EI}.
$$

### 3.2 Triangular: $q(x) = w\,x/L$

Shear and moment are
$$
V(x) = \frac{w}{2L}(L^2 - x^2), \qquad M(x) = \int_x^L V(\xi)\,\mathrm{d}\xi = \frac{w}{6L}\bigl(2L^3 - 3L^2 x + x^3\bigr).
$$
Then
$$
\theta(x) = \frac{1}{EI}\int_0^x M(\xi)\,\mathrm{d}\xi = \frac{w}{6LEI}\biggl(2L^3 x - \frac{3L^2 x^2}{2} + \frac{x^4}{4}\biggr),
$$
$$
u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi = \frac{w}{6LEI}\biggl(L^3 x^2 - \frac{L^2 x^3}{2} + \frac{x^5}{20}\biggr).
$$
Evaluating at $x = L$:
$$
\theta_z(L) = \frac{w}{6LEI}\biggl(2L^4 - \frac{3L^4}{2} + \frac{L^4}{4}\biggr) = \frac{wL^3}{8EI}, \qquad u_y(L) = \frac{w}{6LEI}\biggl(L^5 - \frac{L^5}{2} + \frac{L^5}{20}\biggr) = \frac{11wL^4}{120EI}.
$$

### 3.3 Parabolic: $q(x) = w(x/L)^2$

$V(x) = \frac{w}{3L^2}(L^3 - x^3)$ and
$$
M(x) = \int_x^L V(\xi)\,\mathrm{d}\xi = \frac{w}{12L^2}\bigl(3L^4 - 4L^3 x + x^4\bigr).
$$
Then
$$
\theta(x) = \frac{1}{EI}\int_0^x M(\xi)\,\mathrm{d}\xi = \frac{w}{12L^2 EI}\biggl(3L^4 x - 2L^3 x^2 + \frac{x^5}{5}\biggr),
$$
$$
u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi = \frac{w}{12L^2 EI}\biggl(\frac{3L^4 x^2}{2} - \frac{2L^3 x^3}{3} + \frac{x^6}{30}\biggr).
$$
At the tip:
$$
\theta_z(L) = \frac{w}{12L^2 EI}\biggl(3L^5 - 2L^5 + \frac{L^5}{5}\biggr) = \frac{wL^3}{10EI}, \qquad u_y(L) = \frac{w}{12L^2 EI}\biggl(\frac{3L^6}{2} - \frac{2L^6}{3} + \frac{L^6}{30}\biggr) = \frac{13wL^4}{180EI}.
$$

---

## 4. Summary table

| Load case | $\theta_z(L)$ | $u_y(L)$ |
|-----------|----------------|----------|
| Point, end ($a = L$) | $\dfrac{PL^2}{2EI}$ | $\dfrac{PL^3}{3EI}$ |
| Point, mid-span ($a = L/2$) | $\dfrac{PL^2}{8EI}$ | $\dfrac{5PL^3}{48EI}$ |
| Point, quarter-span ($a = L/4$) | $\dfrac{PL^2}{32EI}$ | $\dfrac{11PL^3}{384EI}$ |
| UDL ($q = w$) | $\dfrac{wL^3}{6EI}$ | $\dfrac{wL^4}{8EI}$ |
| Triangular ($q = wx/L$) | $\dfrac{wL^3}{8EI}$ | $\dfrac{11wL^4}{120EI}$ |
| Parabolic ($q = w(x/L)^2$) | $\dfrac{wL^3}{10EI}$ | $\dfrac{13wL^4}{180EI}$ |

These expressions are the reference values used for tip verification in the grid-convergence table (`gci_richardson_roark_deflection_rotation.csv`) and are consistent with the full-field formulae in `DEFORMATION_AND_GCI_NOTES.md`.
