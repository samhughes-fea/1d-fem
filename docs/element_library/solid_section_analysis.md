# 2D section analysis (Giavotto / BECAS)

`SolidSectionAnalysis` gives the stiffness and properties of a beam section from a 2D
finite-element mesh of it (`solid_section_mesh.md`), with 3D anisotropic materials. It
works for any shape: solid, thick-walled, laminated, sandwich, or a mix. The method is
Giavotto et al. (1983) as implemented in BECAS (Blasques and Stolpe 2012). Code:
`pre_processing/section_integration/solid/section_analysis.py`.

It gives the same quantities as the thin-walled analysis, in the same order and right-handed
signs:

- `stiffness` (6 × 6, fully coupled Timoshenko), `flexibility`, `classical_stiffness`;
- `shear_centre`, `tension_centre`;
- `warping_constant` (`E Γ` from the torsional warping about the shear centre);
- `mass_properties`, `geometric_terms` (second-order terms), `axial_points`.

`SectionAnalysis.from_solid` maps them to element signs, so every element that takes a
thin-walled section also takes a 2D one, the Benscoter element included (restrained
warping, below). Third-order and GBT data are not built from 2D sections; cross-section
distortion stays in GBT.

## Kinematics

The section moves rigidly and warps in all three directions, `g = (g_x, g_y, g_z)`. With
`ψ = [ε, κ_y, κ_z, γ_xy, γ_xz, φ]` the 3D strains (Voigt `[xx, yy, zz, yz, xz, xy]`) are

```
ε_xx = ε + z κ_y − y κ_z            + g_x'
ε_yy = ∂g_y/∂y,   ε_zz = ∂g_z/∂z,   γ_yz = ∂g_y/∂z + ∂g_z/∂y
γ_xz = γ_xz + y φ + ∂g_x/∂z         + g_z'
γ_xy = γ_xy − z φ + ∂g_x/∂y         + g_y'
```

or `ε = Z ψ + B ĝ + N_x ĝ'`. The in-plane warping lets the section deform freely in its own
plane (Poisson contraction, anticlastic bending of walls across their width). The
thin-walled walls with `hoop="free"` do the same, and that is the option the 2D results
match.

## Solution

With `E = ∫BᵀCB`, `C_x = ∫BᵀCN_x`, `R = ∫BᵀCZ`, `L = ∫N_xᵀCZ` and `A = ∫ZᵀCZ`, the two
linear systems are those of the thin-walled analysis with its `Y` terms zero:

```
[E   R  D] [g1  ]   [0    ]        [E   R  D] [g0  ]   [−(C_x − C_xᵀ) g1 + L ψ1]
[Rᵀ  A  0] [ψ1  ] = [T F0 ]        [Rᵀ  A  0] [ψ0  ] = [F0 − Lᵀ g1             ]
[Dᵀ  0  0] [λ1  ]   [0    ]        [Dᵀ  0  0] [λ0  ]   [0                      ]
```

`D` holds the six warping fields equal to rigid motions: the three translations, the
in-plane rotation `(0, −z, y)`, and `g_x = z`, `g_x = −y`. `T` gives the force rates
(`M_y' = V_z`, `M_z' = −V_y`). With `ε = Φ F0` the flexibility is `∫ΦᵀCΦ dA`. Tied nodes
(where a thin-walled strip's end joins another wall) are eliminated before the solve. The
systems are sparse and solved by LU.

`axial_points` (pretwist coupling) use the axial modulus `E_x`: `C_xx` with the in-plane
stresses `σ_yy, σ_zz, τ_yz` free and the shear strains held, as the axial term of the
thin-walled wall law.

## Restrained warping (the 7th DOF)

`benscoter_stiffness` (8 × 8, order `[ψ6, d, ψ']` as the thin-walled one) adds the torsional
warping amplitude `ψ`. The warping is the Saint-Venant warping of `ψ6` plus
`ω (ψ − φ)`, with `ω` the torsional warping (all three components, about the shear centre,
its extension and bending parts removed):

- `d = ψ − φ` strains the section in shear, `ε = B ω d`;
- `ψ'` stretches it, `ε = N_x ω ψ' + B h ψ'`, where `h` lets the section relax in its own
  plane (Poisson contraction). Without `h` the axial stiffness of a ±45 ply would be the
  laterally constrained `C_xx`, 7 % too stiff on the CFRP I-section below.

`warping_constant` is `K[7, 7]`. `warping_inertia` is `∫ρ ω_x² dA`, and
`benscoter_geometric_terms` gives the bimoment's second-order terms.

The 2D warping also varies through the wall thickness and follows `ψ`, so `K[6, 6]` (warping
shear) is higher than the Kirchhoff walls' value, by 30 % on the I-sections below. Kirchhoff
walls cannot shear through their thickness. `K[7, 7]` agrees.

## 3D stresses

`point_stresses(forces=F)` (or `benscoter=e8`) gives a `PointStress` at every quadrature
point, and `stress_at(element, ξ, η, ...)` gives one anywhere in an element. Each holds:

- the strain and stress in section axes `[xx, yy, zz, yz, xz, xy]`;
- `stress_material`: ply or core axes `[11, 22, 33, 23, 13, 12]`;
- `stress_wall`: wall axes `[xx, ss, nn, sn, xn, xs]`; `interlaminar` is `(σ_nn, τ_sn, τ_xn)`.

For section forces, the forces vary along the beam as with no distributed load. The transverse
shear and interlaminar stresses of the bending gradient are therefore in the state. On a
`from_thin_walled` mesh `η` runs through the thickness, so a ply interface is `η = +1` of
one row and `η = −1` of the next.

### Ply interfaces

The raw stresses of a quadratic element on its edge are its least accurate values: with two
element rows per ply, `τ_xn` on a cross-ply interface is up to 15 % off, and six rows are
needed for 2 %. `interface_stresses(forces=F)` (or `benscoter=e8`) recovers them instead.
On a layered `Q8` mesh (`from_thin_walled`, `section_mesh_from_cut`) `ply_interfaces()`
lists the element pairs whose shared edge changes material, ply angle or wall direction.
On each side, the interlaminar stresses at the Gauss lines of the ply's two rows nearest
the interface are fitted by a quadratic in the distance from it. With one row per ply the
fit is linear. Each fit is taken to the interface, and the two sides are averaged.
`interface_stress(below, above, ξ)` gives one point. Two rows per ply are then exact for
the quadratic shear of a layered beam.

From a solved beam, `solid/element_recovery.py` takes the section state from Euler-Bernoulli,
Timoshenko or Benscoter elements, corotational ones included. Linear static jobs,
nonlinear static jobs whose elements are all corotational, and eigen jobs (one block of
rows per mode, column `mode`) write `secondary_results/section_stresses/section_stresses.csv`
and `interface_stresses.csv`. They are written at every element's mid-point, or at the
stations of `stress_stations.txt` (`post_processing/stations.py`): a row `[x]` is where the
axis crosses the plane `X = x`, and a row `[x] [y] [z]` is the nearest point of the axis.
The `x` column is the global `X` of the axis point, and `xi` is its element coordinate.
A linear Timoshenko element's shear strain is exact only at its mid-point, so put
stations for shear and interlaminar stresses at element mid-points.

## Checks (`tests/test_solid_section_analysis.py`)

| Case | Result |
|---|---|
| Isotropic 2:1 rectangle | `EA`, `EI` exact; `J = 0.22868 h b³` (series 0.228682) |
| `ν = 0` rectangle | shear coefficient `5/6` to 1e-4 |
| Homogeneous off-axis rectangle, offset | flexibility under `N`, `M_y`, `M_z` exact, including the extension-shear and bending-twist couplings; both centres at the centroid |
| CFRP I-section, `[0/45/−45/0]` flanges, `[45/−45]s` web | `EA`, `EI`, `GJ` within 0.2 % of thin-walled; `E Γ` within 2 %; flange-direction shear stiffer by 6 % (the 2D web shears through its thickness, Kirchhoff walls cannot) |
| CAS box `[30]₄` | the whole 6 × 6, couplings included, within 0.8 % of thin-walled |
| Timoshenko cantilever job on a coupled 2D section | tip deflection and twist from the section flexibility to 0.05 % |

Stresses and restrained warping (`tests/test_solid_section_stresses.py`):

| Case | Result |
|---|---|
| Rectangle, `M_y` and `V_z` (`ν = 0`) | `σ_xx = M z / I` exact; `τ_xz` the parabola to 1 % of its peak |
| Laminated I-section | the stresses integrate back to all six section forces |
| `[0/90/0/90]` strip under `V_z` | interlaminar `τ_xn` at each interface, from both plies, within 2 % of layered-beam theory |
| Sandwich, aluminium faces on a soft core | core shear uniform at `V / (b d)` within 3 %; face stress `M / (b t_f d)` |
| Isotropic I-section | `K[7, 7]` is Vlasov's `E t_f b³ h² / 24` to 0.1 % |
| CFRP I-section | `K[7, 7]` within 0.2 % of thin-walled; warping-shear coupling within 10 % |
| Benscoter cantilever under tip torque | tip twist within 2 % of the thin-walled section |

## Job input

`section_mesh_2d.txt` (format in `solid_section_mesh.md`) defines the 2D sections and
`element_section.txt` assigns them, alongside thin-walled sections if there are any (the
section ids must differ). Material ids not in its `[Solid_Material]` block come from the
`[Ply_Material]` block of `laminate.txt`. GBT elements and wall loads need thin-walled
sections.
