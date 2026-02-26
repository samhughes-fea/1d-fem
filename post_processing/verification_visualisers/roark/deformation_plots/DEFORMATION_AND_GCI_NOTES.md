# Verification of Euler–Bernoulli Beam Solutions: Reference Formulae, Spatial Comparison, and Grid Convergence

**Abstract.** This note supports the verification of a finite-element Euler–Bernoulli beam implementation by (i) deriving the analytical reference solutions used for comparison (Roark-type formulae for cantilever point and distributed loads), (ii) interpreting the full-span deformation plot that compares finite-element and reference transverse displacement and rotation, and (iii) explaining the Grid Convergence Index (GCI) and Richardson extrapolation table, including the role of the safety factor, the asymptotic relation between GCI values, and the circumstances under which reported GCI and observed order are not meaningful. The reference solutions and the GCI procedure are those implemented in the accompanying scripts and summarised in the file `gci_richardson_roark_deflection_rotation.csv` and the figure `deformation_convergence_uy_theta_all_loads.png`.

---

## 1. Introduction

Verification of a numerical beam model requires a reference solution and defined metrics. Here, the reference is provided by closed-form solutions for a linear, prismatic Euler–Bernoulli cantilever under point and distributed loads—formulae that appear in Roark’s *Formulas for Stress and Strain* and equivalent texts. Two verification outputs are discussed: a spatial comparison of transverse displacement $u_y$ and cross-section rotation $\theta_z$ along the beam (`deformation_convergence_uy_theta_all_loads.png`), and a grid-convergence and tip-verification table (`gci_richardson_roark_deflection_rotation.csv` and its LaTeX/PDF table). This document derives the reference relations, defines the quantities in the table, and clarifies the interpretation of the GCI and related entries (including when they are suppressed in the report).

---

## 2. Reference Solutions: Euler–Bernoulli Beam and Roark-Type Formulae

### 2.1 Governing relations

For a prismatic Euler–Bernoulli beam of length $L$, flexural rigidity $EI$, and coordinate $x$ from the fixed end, the bending moment $M(x)$, curvature $\kappa(x)$, slope $\theta(x) = \theta_z(x)$, and transverse displacement $u(x) = u_y(x)$ satisfy
$$
\kappa = \frac{M}{EI}, \qquad \frac{\mathrm{d}\theta}{\mathrm{d}x} = \kappa, \qquad \frac{\mathrm{d}u}{\mathrm{d}x} = \theta,
$$
with $u(0) = \theta(0) = 0$. Integrating from the fixed end yields
$$
\theta(x) = \int_0^x \frac{M(\xi)}{EI}\,\mathrm{d}\xi, \qquad u(x) = \int_0^x \theta(\xi)\,\mathrm{d}\xi.
$$
Shear and moment are related to the applied load $q(x)$ by $\mathrm{d}V/\mathrm{d}x = -q$ and $\mathrm{d}M/\mathrm{d}x = V$, with sign conventions such that a downward load and positive $u$ downward are consistent with the finite-element output.

### 2.2 Point load at $x = a$

A concentrated force $P$ (downward positive) at $x = a$ gives
$$
M(x) = P(a - x) \quad \text{for } x \leq a, \qquad M(x) = 0 \quad \text{for } x > a.
$$
For $x \leq a$,
$$
\theta(x) = \frac{P}{2EI}\,x(2a - x), \qquad u(x) = \frac{P}{6EI}\,x^2(3a - x).
$$
For $x > a$, the moment is zero so $\theta(x) = \theta(a)$ and $u(x) = u(a) + \theta(a)(x - a)$, with
$$
\theta(a) = \frac{Pa^2}{2EI}, \qquad u(a) = \frac{Pa^2}{6EI}(3a - a) = \frac{Pa^3}{3EI}.
$$
The implementations used for verification take $a = L$ (end), $a = L/2$ (mid-span), and $a = L/4$ (quarter-span); the same expressions are evaluated at the finite-element node positions for the deformation plot and on a fine 1000-point grid for tip values in the GCI table.

### 2.3 Distributed loads

For distributed load intensity $q(x)$ (positive downward), shear and moment are
$$
V(x) = \int_x^L q(\xi)\,\mathrm{d}\xi, \qquad M(x) = \int_x^L V(\xi)\,\mathrm{d}\xi.
$$
Then $\theta(x) = \int_0^x (M/EI)\,\mathrm{d}\xi$ and $u(x) = \int_0^x \theta\,\mathrm{d}\xi$. The three profiles used are:

- **Uniform (UDL):** $q(x) = w$. Then $V(x) = w(L - x)$, $M(x) = \tfrac{1}{2}w(L-x)^2$, $\theta(x) = \tfrac{w}{6EI}(3L^2 x - 3L x^2 + x^3)$, and $u(x) = \tfrac{w}{24EI}(6L^2 x^2 - 4L x^3 + x^4)$.
- **Triangular:** $q(x) = w\,x/L$. Then $V(x) = w(L^2 - x^2)/(2L)$, $M(x)$ and thence $\theta$, $u$ by integration.
- **Parabolic:** $q(x) = w(x/L)^2$. Then $V(x) = w(L^3 - x^3)/(3L^2)$, and $M$, $\theta$, $u$ follow accordingly.

In the code, $V$ and $M$ are implemented in closed form and $\theta$, $u$ are obtained by cumulative quadrature over a sorted $x$ grid so that the reference is consistent with the Euler–Bernoulli relations above. These reference fields are used for the deformation plot (evaluated or interpolated at the same $x$ as the finite-element solution) and for the tip values in the GCI table (evaluated on a 1000-point grid so the reference tip is effectively exact).

---

## 3. Spatial Verification: Deformation Plot

The figure `deformation_convergence_uy_theta_all_loads.png` presents a single row of two panels: transverse displacement $u_y$ (mm) and cross-section rotation $\theta_z$ (deg) along the beam axis $x$. The reference and the finite-element solution use the same units (mm for $u_y$, degrees for $\theta_z$), so the overlay is a direct comparison. For each of six load cases—end, mid-span, and quarter-span point loads, and UDL, triangular, and parabolic distributed loads—the finite-element solution from the **finest available mesh** (largest $n$ per case) is overlaid with the corresponding Roark-type reference derived in §2. The reference is evaluated (or interpolated) at the same abscissae as the finite-element solution.

Spatial agreement between the numerical and analytical curves over the full span confirms that the implemented shape functions, assembly, and boundary conditions reproduce the correct displacement and rotation fields. The quantities $u_y$ and $\theta_z$ are those for which the reference provides explicit or numerically integrated values; other degrees of freedom are not compared in this figure. This plot therefore provides **field-level** verification, complementary to the **point-wise** (tip) verification in the GCI table.

---

## 4. Grid Convergence Index and Richardson Extrapolation

### 4.1 Three-grid procedure and definitions

The table summarises a **three-grid** study with mesh sizes $h_1 < h_2 < h_3$ corresponding to $n_1 = 100$, $n_2 = 50$, and $n_3 = 25$ elements and refinement ratio $r = h_2/h_1 = h_3/h_2 = 2$. Denote the discrete solution (tip deflection or tip rotation) on each grid by $\phi_1$, $\phi_2$, $\phi_3$; thus $\phi_1$ is the solution on the finest grid ($n_1 = 100$) and $\phi_3$ on the coarsest ($n_3 = 25$). The **observed order** $p$ is obtained from the ratio of solution changes:
$$
\frac{\phi_3 - \phi_2}{\phi_2 - \phi_1} = r^p \quad \Rightarrow \quad p = \frac{\ln\bigl|\frac{\phi_3 - \phi_2}{\phi_2 - \phi_1}\bigr|}{\ln r}.
$$
The **Richardson extrapolant** is
$$
\phi_\mathrm{ext} = \phi_1 + \frac{\phi_1 - \phi_2}{r^p - 1}.
$$
The **Grid Convergence Index** (GCI) on the fine grid (GCI$_{12}$) and on the medium grid (GCI$_{23}$) is defined with a **safety factor** $F_s$ (Celik et al., 2008; ASME V&V 20–2009). For the three-grid case, the recommended value is $F_s = 1.25$. The formulae used in the report are
$$
\mathrm{GCI}_{12} = \frac{F_s}{r^p - 1}\,\left|\frac{\phi_2 - \phi_1}{\phi_1}\right| \times 100\,\%, \qquad \mathrm{GCI}_{23} = \frac{F_s}{r^p - 1}\,\left|\frac{\phi_3 - \phi_2}{\phi_2}\right| \times 100\,\%.
$$
Thus GCI$_{12}$ and GCI$_{23}$ are percentage estimates of the relative error on the fine and medium grids, inflated by $F_s$ to account for uncertainty in $p$ and departure from the asymptotic range. The factor 1.25 is the one recommended for three-grid studies so that the reported GCI is a conservative error bound when the observed order is used.

### 4.2 Asymptotic range and the ratio GCI$_{23}$/GCI$_{12}$

When the discretisation error is in the **asymptotic range**, the error on successive grids scales as $h^p$. With constant refinement ratio $r$, the ratio of the two GCI values satisfies
$$
\frac{\mathrm{GCI}_{23}}{\mathrm{GCI}_{12}} \approx r^p.
$$
For $r = 2$ and $p = 2$ (second-order scheme), $r^p = 4$. Therefore a ratio close to 4 supports that the solutions are in the asymptotic regime and that the observed order $p \approx 2$ is consistent with the scheme. The table reports GCI$_{23}$/GCI$_{12}$, $p$, and $r^p$ to allow this check. When the ratio is not close to $r^p$ (e.g. because the fine and medium grids already agree to machine precision), the ratio and $p$ can be meaningless or wildly large; in that case the report shows “---” for those columns (see §4.3).

### 4.3 When GCI and observed order are not reported

For **point-load** cases (end, mid-span, quarter-span), the solutions converge very quickly with mesh refinement. The differences $\phi_2 - \phi_1$ and $\phi_3 - \phi_2$ can then be at the level of numerical noise, so that:

- The ratio $(\phi_3 - \phi_2)/(\phi_2 - \phi_1)$ is ill-defined or enormous, giving an absurd “observed” $p$.
- GCI$_{12}$ and GCI$_{23}$ become negligible or numerically unreliable, and their ratio is not interpretable.

In those rows, the table leaves GCI, the ratio, $p$, and $r^p$ as “---” to avoid misleading values. The **primary verification metric** for all cases is the **percentage difference of the fine-grid solution from the Roark reference**:
$$
\Delta(\phi_1,\mathrm{Roark}) = 100\,\frac{\phi_1 - \phi_\mathrm{Roark}}{\phi_\mathrm{Roark}}\,\%.
$$
When this is near zero, the fine-grid tip value is in agreement with the analytical reference regardless of whether the three-grid GCI is meaningful. The column $\phi_\mathrm{rich}$ (Richardson extrapolant) is retained for completeness; when convergence is already very good, $\phi_\mathrm{rich}$ and $\phi_1$ are nearly equal.

### 4.4 Summary of the table

The table lists, per job and quantity of interest (tip $u_y$ or tip $\theta_z$): the fine-grid value $\phi_1$; GCI$_{12}$ and GCI$_{23}$ (or “---”); their ratio, $p$, and $r^p$ (or “---”); $\Delta(\phi_1,\mathrm{Roark})$; and $\phi_\mathrm{rich}$. For distributed loads (UDL, triangular, parabolic), the reported $p \approx 2$ and ratio $\approx 4$ indicate second-order convergence and asymptotic behaviour; for point loads, the “---” entries reflect that the error is already negligible and the three-grid procedure does not yield interpretable GCI or order. In all cases, small $\Delta(\phi_1,\mathrm{Roark})$ confirms agreement with the Roark reference at the tip.

---

## 5. Conclusions

The deformation plot and the GCI–Roark table together provide two levels of verification: (i) **spatial** agreement of $u_y(x)$ and $\theta_z(x)$ with the derived Roark-type solutions over the full beam, and (ii) **convergence and tip** agreement via the three-grid GCI procedure and the percentage difference from the reference. The use of the safety factor $F_s = 1.25$, the relation GCI$_{23}$/GCI$_{12}$ $\approx r^p$ in the asymptotic range, and the suppression of GCI and $p$ when the solution is already converged are consistent with standard practice and ensure that the table is interpretable and suitable for reporting in a verification context.

---

## References

- ASME V&V 20–2009. *Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer*. American Society of Mechanical Engineers, 2009.
- Celik, I. B., Ghia, U., Roache, P. J., Freitas, C. J., Coleman, H., & Raad, P. E. (2008). Procedure for estimation and reporting of uncertainty due to discretization in CFD applications. *Journal of Fluids Engineering*, 130(7), 078001.
- Young, W. C., Budynas, R. G., & Sadegh, A. M. (2012). *Roark’s Formulas for Stress and Strain* (8th ed.). McGraw-Hill.
