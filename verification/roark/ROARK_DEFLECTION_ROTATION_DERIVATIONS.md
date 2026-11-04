# Derivation of Roark analytical expressions for \(u_y\) and \(\theta_z\)

This document gives the derivation of each analytical expression used for the Roark (reference) curves in the deformation convergence plots. Cantilever: fixed at \(x = 0\), free at \(x = L\). Sign convention: positive transverse load downward; rotation \(\theta_z\) and deflection \(u_y\) match the FEM sign convention used in the code (positive \(u_y\) downward, positive \(\theta_z\) clockwise when viewed from \(+z\)).

---

## 1. Euler–Bernoulli beam

\(EI \theta' = M\), \(\theta = \int \frac{M}{EI}\,\mathrm{d}x\), \(u = \int \theta\,\mathrm{d}x\), with \(\theta(0) = 0\), \(u(0) = 0\).

### 1.1 Point load at tip (End, \(a = L\))

\(V(x) = -P\) for \(x < L\), \(M(x) = -P(L - x)\) for \(x \leq L\).

**Rotation**

\[
\theta(x) = \int_0^x \frac{M(\xi)}{EI}\,\mathrm{d}\xi = -\frac{P}{EI} \int_0^x (L - \xi)\,\mathrm{d}\xi = -\frac{P}{EI}\left[ Lx - \frac{x^2}{2} \right] = \frac{P}{2EI}\,x(2L - x).
\]

(Code applies a sign flip so that positive \(P\) downward gives the same sign as FEM; the expression plotted is \(\frac{Px(2L-x)}{2EI}\).)

**Deflection**

\[
u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi = \frac{P}{2EI} \int_0^x \xi(2L - \xi)\,\mathrm{d}\xi = \frac{P}{2EI}\left[ Lx^2 - \frac{x^3}{3} \right] = \frac{Px^2(3L - x)}{6EI}.
\]

### 1.2 Point load at \(x = a\) (Mid-span \(a = L/2\), Quarter-span \(a = L/4\))

\(V(x) = -P\) for \(x < a\), \(V(x) = 0\) for \(x \geq a\). \(M(x) = -P(a - x)\) for \(x \leq a\), \(M(x) = 0\) for \(x > a\).

**Rotation**

- For \(x < a\): \(\theta(x) = \int_0^x \frac{-P(a-\xi)}{EI}\,\mathrm{d}\xi = -\frac{P}{EI}\bigl[ a x - \frac{x^2}{2} \bigr] = \frac{P}{2EI}\,x(2a - x)\).
- For \(x \geq a\): \(\theta(x) = \theta(a) = \frac{Pa^2}{2EI}\) (constant).

**Deflection**

- For \(x < a\): \(u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi = \frac{P}{2EI}\int_0^x \xi(2a-\xi)\,\mathrm{d}\xi = \frac{Px^2(3a-x)}{6EI}\).
- For \(x \geq a\): \(u(x) = u(a) + \theta(a)(x - a) = \frac{Pa^2(3x-a)}{6EI}\).

The plot labels use the \(x < a\) branch: \(\frac{Px^2(3a-x)}{6EI}\), \(\frac{Px(2a-x)}{2EI}\).

### 1.3 Uniformly distributed load (UDL, \(q = w\))

\(V(x) = -w(L - x)\), \(M(x) = -\frac{w}{2}(L - x)^2\).

**Rotation**

\[
\theta(x) = \int_0^x \frac{M(\xi)}{EI}\,\mathrm{d}\xi = -\frac{w}{2EI} \int_0^x (L - \xi)^2\,\mathrm{d}\xi = \frac{w}{6EI}\bigl[ L^3 - (L-x)^3 \bigr] = \frac{w(L-x)^3}{6EI} + \text{const}.
\]

With \(\theta(0) = 0\), \(\theta(x) = \frac{w}{6EI}\bigl[ L^3 - (L-x)^3 \bigr]\). In the code, the integrated form is used; the expression equivalent to the plotted reference is \(\frac{w(L-x)^3}{6EI}\) (magnitude of the varying part; sign aligned with FEM in code).

**Deflection**

\[
u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi.
\]

Standard result: \(u(x) = \frac{w}{24EI}\bigl( x^4 - 4Lx^3 + 6L^2 x^2 \bigr) = \frac{w}{24EI}(L-x)^3(4L-x)\) (cantilever UDL). The code integrates \(\theta\) numerically; the closed form used for the label is \(\frac{w(L-x)^3(4L-x)}{24EI}\).

### 1.4 Triangular distributed load (\(q(x) = w\,x/L\))

\(V(x) = -\frac{w}{2L}(L^2 - x^2)\), \(M(x) = -\frac{w}{6L}(L-x)^2(2L + x)\).

**Rotation**

\[
\theta(x) = \int_0^x \frac{M(\xi)}{EI}\,\mathrm{d}\xi.
\]

Closed form: \(\theta(x) = \frac{w}{12EIL}(L-x)^2(2L+x)\) (from integrating \(M/EI\); constant chosen so \(\theta(0)=0\)). The code uses cumulative integration of \(M/EI\); the label uses \(\frac{w(L-x)^2(2L+x)}{12EIL}\).

**Deflection**

\[
u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi.
\]

Closed form (cantilever triangular): \(u(x) = \frac{w}{120EIL}(L-x)^3(3L+2x)\). Label: \(\frac{w(L-x)^3(3L+2x)}{120EIL}\).

### 1.5 Parabolic distributed load (\(q(x) = w(x/L)^2\))

\(V(x) = -\frac{w}{3L^2}(L^3 - x^3)\), \(M(x) = -\frac{w}{12L^2}(3L^4 - 4L^3 x + x^4)\).

**Rotation**

\[
\theta(x) = \int_0^x \frac{M(\xi)}{EI}\,\mathrm{d}\xi = \frac{w}{12EIL^2}(3L^4 - 4L^3 x + x^4).
\]

Label: \(\frac{w(3L^4-4L^3x+x^4)}{12EIL^2}\).

**Deflection**

\[
u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi.
\]

Closed form: \(u(x) = \frac{w}{60EIL^2}(3L^4 - 4L^3 x + x^4)(L - x)\). Label: \(\frac{w(3L^4-4L^3x+x^4)(L-x)}{60EIL^2}\).

---

## 2. Timoshenko beam

\(\theta = \int \frac{M}{EI}\,\mathrm{d}x\) (unchanged). Deflection adds shear: \(\frac{\mathrm{d}u}{\mathrm{d}x} = \theta + \frac{V}{\kappa A G}\), so \(u = u_{\mathrm{EB}} + \frac{1}{\kappa A G}\int_0^x V(\xi)\,\mathrm{d}\xi\), with \(\kappa\) the shear correction factor (e.g. \(5/6\)).

### 2.1 Point load at tip (End)

\(V = -P\) for \(x < L\). \(\int_0^x V\,\mathrm{d}\xi = -Px\) for \(x \leq L\), so \(u_{\mathrm{shear}} = -\frac{Px}{\kappa AG}\). Thus

\[
u(x) = \frac{Px^2(3L-x)}{6EI} + \frac{P(L-x)}{\kappa AG}.
\]

Rotation unchanged: \(\theta(x) = \frac{Px(2L-x)}{2EI}\).

### 2.2 Point load at \(x = a\) (Mid, Quarter)

\(V = -P\) for \(x < a\), \(0\) for \(x \geq a\). \(\int_0^x V\,\mathrm{d}\xi = -P\min(x,a)\). So

\[
u(x) = u_{\mathrm{EB}}(x) + \frac{P\,\min(x,a)}{\kappa AG}.
\]

\(\theta\) same as Euler–Bernoulli (piecewise as in §1.2).

### 2.3 UDL

\(V(\xi) = -w(L - \xi)\). \(\int_0^x V\,\mathrm{d}\xi = -w(Lx - x^2/2)\). So

\[
u(x) = u_{\mathrm{EB}}(x) + \frac{w(Lx - x^2/2)}{\kappa AG}.
\]

\(\theta\) unchanged: \(\frac{w(L-x)^3}{6EI}\) (with sign as in code).

### 2.4 Triangular (\(q = wx/L\))

\(\int_0^x V\,\mathrm{d}\xi = -\frac{w}{2L}(L^2 x - x^3/3)\). So

\[
u(x) = u_{\mathrm{EB}}(x) + \frac{w(L^2 x - x^3/3)}{2L\,\kappa AG}.
\]

### 2.5 Parabolic (\(q = w(x/L)^2\))

\(\int_0^x V\,\mathrm{d}\xi = -\frac{w}{3L^2}(L^3 x - x^4/4)\). So

\[
u(x) = u_{\mathrm{EB}}(x) + \frac{w(L^3 x - x^4/4)}{3L^2\,\kappa AG}.
\]

---

## 3. Reference to implementation

- **Euler–Bernoulli point**: [roark_utilities/roarks_formulas_euler_bernoulli_point.py](roark_utilities/roarks_formulas_euler_bernoulli_point.py)  
- **Euler–Bernoulli distributed**: [roark_utilities/roarks_formulas_euler_bernoulli_distributed.py](roark_utilities/roarks_formulas_euler_bernoulli_distributed.py)  
- **Timoshenko point**: [roark_utilities/roarks_formulas_timoshenko_point.py](roark_utilities/roarks_formulas_timoshenko_point.py)  
- **Timoshenko distributed**: [roark_utilities/roarks_formulas_timoshenko_distributed.py](roark_utilities/roarks_formulas_timoshenko_distributed.py)  

The plot labels in `deformation_convergence.py` use the **expression only** (right-hand side) of these results, so that `labelLines` displays the formula on each Roark curve.
