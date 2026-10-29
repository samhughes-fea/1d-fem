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
thin-walled section also takes a 2D one. Benscoter, third-order and GBT data are not
built from 2D sections yet.

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

`E Γ` uses the axial modulus `E_x`: `C_xx` with the in-plane stresses
`σ_yy, σ_zz, τ_yz` free and the shear strains held. This matches the axial term of the
thin-walled wall law. The same modulus weights `axial_points` (pretwist coupling).

## Checks (`tests/test_solid_section_analysis.py`)

| Case | Result |
|---|---|
| Isotropic 2:1 rectangle | `EA`, `EI` exact; `J = 0.22868 h b³` (series 0.228682) |
| `ν = 0` rectangle | shear coefficient `5/6` to 1e-4 |
| Homogeneous off-axis rectangle, offset | flexibility under `N`, `M_y`, `M_z` exact, including the extension-shear and bending-twist couplings; both centres at the centroid |
| CFRP I-section, `[0/45/−45/0]` flanges, `[45/−45]s` web | `EA`, `EI`, `GJ` within 0.2 % of thin-walled; `E Γ` within 2 %; flange-direction shear stiffer by 6 % (the 2D web shears through its thickness, Kirchhoff walls cannot) |
| CAS box `[30]₄` | the whole 6 × 6, couplings included, within 0.8 % of thin-walled |
| Timoshenko cantilever job on a coupled 2D section | tip deflection and twist from the section flexibility to 0.05 % |

## Job input

`section_mesh_2d.txt` (format in `solid_section_mesh.md`) defines the 2D sections and
`element_section.txt` assigns them, alongside thin-walled sections if there are any (the
section ids must differ). Material ids not in its `[Solid_Material]` block come from the
`[Ply_Material]` block of `laminate.txt`. GBT elements and wall loads need thin-walled
sections.
