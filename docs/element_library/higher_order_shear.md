# Third-order shear: Reddy and Levinson elements

`LinearReddyBeamElement3D` and `LinearLevinsonBeamElement3D` share one kinematic
model and one element, `ThirdOrderShearBeamElement3D`
(`pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/third_order_beam.py`).
They differ only in the weak form. Both take the solid-rectangle inputs of
`section.txt` or a thin-walled section from `element_section.txt`.

## 1. Kinematics from the solid rectangle to walls

Reddy's beam (1984) puts a cubic through the depth of a rectangle:

\[
u_x = z\,\varphi - c_1 z^3(\varphi + w'), \qquad c_1 = \frac{4}{3h^2},
\]

so the shear strain \(\gamma_{xz} = (1 - 4z^2/h^2)(\varphi + w')\) is parabolic and
vanishes on the faces. Split the cubic into its part along \(z\) and the rest:
\(z^3 = \tfrac{3h^2}{20} z + W_3(z)\), where \(W_3 = z^3 - \tfrac{3h^2}{20}z\) is
orthogonal to \(1\) and \(z\). With \(\bar\varphi = \varphi - 0.2(\varphi + w')\),

\[
u_x = z\,\bar\varphi - c_1 W_3(z)\,(\varphi + w'), \qquad
\bar\gamma = \bar\varphi + w' = 0.8(\varphi + w').
\]

The shear energy is \(\int G(1 - 4z^2/h^2)^2 (\varphi + w')^2 dA = \tfrac56 GA\,\bar\gamma^2\),
so in these variables Reddy's beam is the **Timoshenko beam with \(\kappa = 5/6\)**
plus a warping field \(\Psi(z)\bar\gamma\), \(\Psi = -\tfrac{5}{3h^2}W_3\), whose
amplitude is tied to the shear strain. \(\Psi\) is exactly the Saint-Venant
warping of a rectangle per unit shear strain: under a constant shear force the
warping is constant along the beam and stretches nothing.

That gives the general definition used here, for any section:

- The section moves as the Saint-Venant (Timoshenko) section: plane rotation plus
  the Saint-Venant warping of its six strains, with the coupled 6 × 6 stiffness.
- The shear part of that warping, \(\Psi(s)\,\boldsymbol\gamma\), is slaved to the
  local shear strain. Where \(\boldsymbol\gamma\) changes along the beam faster
  than in the Saint-Venant solution it stretches the section:

\[
\varepsilon_{xx} \mathrel{+}= \Psi_y(s)\,\chi_y + \Psi_z(s)\,\chi_z, \qquad
\boldsymbol\chi = \boldsymbol\gamma' - R\,\boldsymbol\psi_6 .
\]

\(R\) (2 × 6) is the shear-strain rate of the Saint-Venant solution itself
(\(\boldsymbol\psi_6' = C\,T\,K\,\boldsymbol\psi_6\) with no distributed load);
it is zero for an uncoupled section and keeps Saint-Venant states free of
\(\chi\) for coupled ones.

On thin walls \(\Psi = W_0 K_{:,3:5}\) from the wall analysis (PR 2), with its
rigid part (\(1, y, z\)) removed in the walls' axial stiffness so the slaved
warping carries no \(N\), \(M_y\) or \(M_z\). On a thin web this reproduces
Reddy's cubic: the test `test_rectangle_web_gives_reddy_section_term` finds
\(EI/84\) to 0.3 %.

## 2. Section terms

`SectionAnalysis.third_order_stiffness` (8 × 8) acts on
\(\mathbf e_8 = [\varepsilon, \kappa_y, \kappa_z, \gamma_{xy}, \gamma_{xz}, \phi, \chi_y, \chi_z]\):

\[
\mathbf D_8 = \int \mathbf G_8(s)^{\mathsf T}\,\mathbf K_{wall}\,\mathbf G_8(s)\,ds,\qquad
\mathbf G_8 = [\,\Phi(s) K \;|\; \mathbf N_x(s)\Psi\,],
\]

where \(\Phi(s)K\) are the Saint-Venant wall strains per unit beam strain. The
leading 6 × 6 block is the Saint-Venant stiffness. For the solid rectangle of
`section.txt` the terms are Reddy's \(z^4\) and \(z^6\) integrals reduced to

\[
D_8[\chi_y,\chi_y] = \frac{EI_z}{84}, \qquad D_8[\chi_z,\chi_z] = \frac{EI_y}{84},
\]

from \(E\int\Psi^2 dA = E\,b\,(5/3h^2)^2\,h^7/2800 = E\,b h^3/1008\). Other solid
shapes are treated as rectangles (as Reddy's cubic assumes); a thin-walled
section gives its own terms. The old `[alpha]` column of `section.txt` is no
longer used.

## 3. Element

- **DOFs per node:** \([u, v, w, \theta_x, \theta_y, \theta_z, v', w']\). \(\chi\)
  needs \(\gamma' = v'' - \theta_z'\) and \(w'' + \theta_y'\), so \(v\) and \(w\) are
  cubic Hermite with their slopes as the two extra DOFs (`Q1`, `Q2`; the mesh
  stride becomes 8). \(u\) and the rotations are linear.
- **Strains** (element signs, \(\kappa_y\) negated as in every 1d-fem beam):
  \(\gamma_{xy} = v' - \theta_z\), \(\gamma_{xz} = w' + \theta_y\) (right-handed
  \(\theta_y\)), then \(\chi\).
- **Quadrature:** 4-point Gauss, exact. The Kirchhoff limit \(\gamma = 0\) leaves
  continuous piecewise-linear rotations, so nothing locks
  (`test_thin_beam_limit_is_euler_bernoulli`, \(L/h = 100\)).
- **Loads:** \(\int N_q^{\mathsf T} q\,dx\) with Hermite \(v\), \(w\), so a distributed
  load also loads the slopes; point loads by \(N_q(x_p)^{\mathsf T}P\).
- **Results:** the Gauss cache stores the six Timoshenko strains with the 6 × 6
  stiffness, so `section_forces` are the Timoshenko resultants. Wall and ply
  stresses (`wall_stresses.csv`) include the \(\Psi\chi\) axial strain.

### Reddy: variational

\(K_e = \int B_8^{\mathsf T} D_8 B_8\,dx\). The higher-order resultant
\(P = D_8[\chi,:]\,\mathbf e_8\) has its own equilibrium equation. A clamp fixes the
slopes as well (`Q1 = Q2 = 0`), which makes \(\gamma = 0\) there and restrains the
shear warping, Reddy's \(w = w' = \varphi = 0\). With the slopes free (warping
free) and no distributed load the element is the Timoshenko element, node for
node.

### Levinson: vectorial equilibrium

Levinson (1981) uses the same kinematics but writes only the classical
equilibrium equations, so the higher-order resultant has no equation of its own.
In weak form that is Petrov-Galerkin, classical virtual strains against the full
resultants:

\[
K_e = \int B_6^{\mathsf T}\,D_8[0{:}6,\,:]\,B_8\,dx .
\]

For any uncoupled section (all isotropic ones) \(D_8[0{:}6, \chi] = 0\) and this is the
Timoshenko beam with the section's shear factor, \(5/6\) for a rectangle. That is the
known result for Levinson's beam; its stresses still follow the higher-order
kinematics. For a section whose walls couple the slaved warping to the classical
strains (unbalanced laminates), \(K_e\) is not symmetric. Levinson's boundary
conditions are the classical ones, so leave `Q1`, `Q2` free at a clamp.

## 4. Checks

`tests/test_third_order_shear.py`:

| Check | Reference | Result |
|---|---|---|
| Thin web, \(\chi_z\) stiffness | \(EI/84\) (Reddy's cubic) | 0.3 % |
| Free slopes, tip load, isotropic and composite box | Timoshenko element | identical to round-off |
| Clamped Reddy, tip load, both planes | \(PL^3/3EI + (P/\kappa GA)(L - \tanh\lambda L/\lambda)\), \(\lambda^2 = 84\kappa GA/EI\) | < 2e-4 at 64 elements, converging |
| Uniform load, free warping | \(qL^4/8EI + qL^2/2\kappa GA\) | 3e-4 |
| Uniform load, clamped Reddy | closed form of \(\kappa GA\gamma - H\gamma'' = q(L - x)\) | 3e-4 |
| Thin limit \(L/h = 100\) | \(PL^3/3EI\) | 0.2 %, no locking |
| Unbalanced composite box | Reddy symmetric, Levinson not | yes |
| Wall stresses, free slopes | Saint-Venant recovery | 1e-6 |
| Thick cantilever \(L/h = 2\), fully clamped | 2D plane-stress Q8 model | Reddy 0.6 % above, Timoshenko 1.1 % above |

The 2D model clamps every point of the root section, restraining warping and
Poisson contraction; the clamped Reddy beam restrains the shear warping, so it
lands between Timoshenko and the 2D result.
