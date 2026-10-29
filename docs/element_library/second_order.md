# Second-order geometric stiffness, P-delta and linear buckling

Every linear beam element except GBT has a geometric stiffness
`linear_geometric_stiffness_matrix(U_e)` built from its current section forces
(`pre_processing/element_library/linear/beam/second_order.py`). It drives two analyses:

- **P-delta** (`[Static] second_order = p_delta`): `(K + K_σ(U)) U = F`, iterated from the
  linear solution (`processing/buckling/p_delta.py`, `LinearStaticSimulationRunner.solve_p_delta`).
- **Linear buckling** (`[Type] buckling`): `(K + λ K_σ(U_ref)) x = 0` about the linear
  state of the reference loads (`processing/buckling/linear_buckling.py`).

Rotations stay small; updating the geometry is the large-displacement analysis.

## 1. Strain energy of a section rigid in its plane

A point at `(y, z)` of the section moves by `u_y = v − z φ`, `u_z = w + y φ`, where `v`,
`w`, `φ = θ_x` are the right-handed DOFs of the reference axis. The quadratic part of the
Green strain is

```
ε_xx(2) = ½ [(v' − z φ')² + (w' + y φ')²]
γ_xs(2) = (v' − z φ')(−t_z φ) + (w' + y φ')(t_y φ)
```

along a wall with tangent `t` (the same with the normal `n` across it, for the Kirchhoff
transverse shear force `Q_x = dM_xx/dx`). Against the current axial stress `σ` and shear
flow `q` its work per length is `½ gᵀ S g` with `g = [v', w', φ', φ]`:

```
S = [[ N,    0,   −S_z, −V_z],
     [ 0,    N,    S_y,  V_y],
     [−S_z,  S_y,  W,    Q_r],
     [−V_z,  V_y,  Q_r,  0  ]]
```

| term | integral | meaning |
|---|---|---|
| `N` | `∫ σ` | axial force |
| `S_y`, `S_z` | `∫ σ y`, `∫ σ z` | first moments, `S_z = M_y`, `S_y = −M_z` |
| `W` | `∫ σ (y² + z²)` | Wagner term |
| `V_y`, `V_z` | `∫ q t_y`, `∫ q t_z` | shear forces |
| `Q_r` | `∫ q (y t_y + z t_z)` | shear-flow moment (zero for a closed cell under uniform flow) |

Everything is about the reference axis, so nothing is shifted by hand:

- **Shear-centre offset**: comes from the coupled Saint-Venant stiffness, whose shear-twist
  terms make the shear-rigid limit twist about the shear centre.
- **Tension-centre offset**: an axial force away from the tension centre is a moment about
  the axis, and the stresses carry it into `S_y`, `S_z` and `W`.
- **Wagner term**: `W` holds the polar term `N I_p/A`, its parallel-axis shift and the Wagner
  coefficients `β M` of the bending moments, all from the stresses.
- **Bimoment**: the restrained-warping stress `E ω ψ'` of a Benscoter element adds to `W`.

`ThinWalledSectionAnalysis.geometric_terms` (7 × 6) gives the terms per unit section force
from the Saint-Venant wall stresses, with `M_xx` through the thickness. Two further tables
give them per unit Benscoter strain `[d, ψ']` and per unit third-order strain `χ`.
`SectionAnalysis` holds all three element-signed. An isotropic `section.txt` section is
taken as doubly symmetric about its centroid (`W = N I_x / A`).

`K_σ = ∫ Gᵀ S G dx` (4-point Gauss) with `G` the element's interpolation of `g`:

| element | `v'`, `w'` | section forces |
|---|---|---|
| Euler-Bernoulli | Hermite slopes | `N`, `M`, `T` from the strains; `V = dM/dx` over the element |
| Timoshenko, Benscoter | chord slope plus the linear part of the rotation (the slope of the quadratic deflection the element implies) | shear strains (and Benscoter `d`) at the centre, as the reduced rule sees them |
| Reddy, Levinson | Hermite slopes (`Q1`, `Q2`) | Timoshenko forces plus the stress of `Ψ χ` |

The Timoshenko choices matter. Shear strains away from a linear element's centre carry the
parasitic part of the bending, and using them inflated the Michell cantilever load by 4 %
at 32 elements. The chord slope alone converged slowly in lateral buckling (−2.5 % at 32
elements).

Euler-Bernoulli elements have no shear strains, so they cannot place the shear centre off
the reference axis. Use Timoshenko or Benscoter elements for such sections. GBT elements
raise `NotImplementedError`, because their distortional modes need a wall-level `K_σ`
(local and distortional buckling).

## 2. P-delta

Fixed-point iteration `U_{k+1} = (K + K_σ(U_k))⁻¹ F`, which contracts at about `P / P_cr`.
`[Static]` keys are `second_order = none | p_delta`, `p_delta_tolerance` (default `1e-8`,
on `‖ΔU‖/‖U‖`) and `p_delta_max_iterations` (50). The job writes
`primary_results/p_delta_iterations.csv` and stops with an error if the iteration does not
converge. Past the first buckling load a fixed point can still exist, but it is unstable
(`K + K_σ` is indefinite).

## 3. Checks

`tests/test_second_order.py` uses the formulas evaluated from mid-line hand properties:

| case | reference | result |
|---|---|---|
| Pinned column, EB / Timoshenko / Reddy | Euler, Engesser | < 0.1 % |
| Channel in compression at the centroid, web, centroid or shear-centre axis | flexural-torsional cubic (Timoshenko & Gere 5.4), `Γ = 0` for Timoshenko, `Γ` for Benscoter | < 0.15 %, same for every axis |
| Doubly symmetric I, uniform moment, EB / Timoshenko / Benscoter | TG 6.3 | < 0.1 % |
| Monosymmetric I, both senses of moment | Trahair eq. 6.3 with the Wagner `β` | < 0.1 % |
| Narrow rectangle: cantilever tip load, simple beam central load | 4.013 and 16.94 `√(EI_z GJ)/L²` | < 0.1 % |
| Unbalanced composite channel | same loads for any reference axis | 1e-6 |
| Beam-column, central lateral load | `3(tan u − u)/u³` amplification | 0.2 % |
| Channel compressed at the web | secant formula with the tension-centre eccentricity | < 0.5 % |

The worked example [examples/thin_walled/channel_buckling](../../examples/thin_walled/channel_buckling/README.md)
runs the same channel through the buckling and P-delta job runners. Comparison with Abaqus
is part of the nonlinear validation.
