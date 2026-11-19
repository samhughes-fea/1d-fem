# Ply failure criteria

`anisotropic_section_model/ply_failure.py` checks ply stresses against four
failure criteria. `src/onedfem/post_processing/failure/ply_failure.py` applies them at the Gauss points
of a solved job.

## Exposure, reserve factor and margin

Every criterion gives an **exposure** \(f\): how far the stress state is along a
proportional load path towards the failure surface (\(f = 1\) on the surface).

- The reserve factor is \(RF = 1/f\) and the margin of safety is \(MS = 1/f - 1\).
- For a criterion that is homogeneous of degree one in the stresses (maximum stress,
  Hashin fibre compression, Puck), \(f\) is the criterion itself.
- For a quadratic one, \(f\) is found from \(a\,RF^2 + b\,RF = 1\), where \(a\) holds the
  quadratic terms and \(b\) the linear ones (Tsai-Wu, Hashin matrix compression). For the
  purely quadratic Hashin modes this is the square root of the failure index.

Either way \(f\) scales linearly with the load. In a linear analysis \(1/f\) at the worst
point is therefore the first-ply-failure load factor.

Stresses are in ply axes (1 fibre, 2 transverse, 3 normal), Voigt order
\([\sigma_{11}, \sigma_{22}, \sigma_{33}, \tau_{23}, \tau_{13}, \tau_{12}]\).

## Criteria

`PlyStrength(Xt, Xc, Yt, Yc, S12, S23=None, Zt=None, Zc=None, F_star=-0.5, p_TL_t=0.30,
p_TL_c=0.25, p_TT=0.25)`. Compressive strengths are magnitudes.

- `S23` defaults to Puck's \(R^A_{\perp\perp} = Y_c / (2(1 + p_{\perp\perp}))\).
- `Zt` and `Zc` default to `Yt` and `Yc` (a transversely isotropic ply).

| Criterion | Form | Modes |
|---|---|---|
| Maximum stress | each stress over its strength | `s11t`, `s11c`, `s22t`, … `t12` |
| Tsai-Wu | \(F_i\sigma_i + F_{ij}\sigma_i\sigma_j = 1\), with \(F_{ij} = F^*\sqrt{F_{ii}F_{jj}}\) for the normal-stress pairs | none |
| Hashin (1980, 3D) | fibre tension or compression, matrix tension or compression | `FT`, `FC`, `MT`, `MC` |
| Puck | fibre failure \(\sigma_{11}/X\); inter-fibre failure on the most stressed action plane | `FFT`, `FFC`, `IFF-A`, `IFF-B`, `IFF-C` |

### Puck's inter-fibre failure

Puck's inter-fibre failure is evaluated on action planes parallel to the fibres at angle
\(\theta\), with \(\sigma_n\), \(\tau_{nt}\) and \(\tau_{n1}\) the stresses on that plane
(Puck and Schürmann 2002; VDI 2014):

- For \(\sigma_n \ge 0\): \(f = \sqrt{[(1/R^+_\perp - p^+_{\perp\psi}/R^A_{\perp\psi})\sigma_n]^2 + (\tau_{nt}/R^A_{\perp\perp})^2 + (\tau_{n1}/R_{\perp\parallel})^2} + (p^+_{\perp\psi}/R^A_{\perp\psi})\,\sigma_n\).
- For \(\sigma_n < 0\): \(f = \sqrt{(\tau_{nt}/R^A_{\perp\perp})^2 + (\tau_{n1}/R_{\perp\parallel})^2 + (p^-_{\perp\psi}\sigma_n/R^A_{\perp\psi})^2} + (p^-_{\perp\psi}/R^A_{\perp\psi})\,\sigma_n\).
- In both, \(p_{\perp\psi}/R^A_{\perp\psi} = (p_{\perp\perp}/R^A_{\perp\perp})\cos^2\psi + (p_{\perp\parallel}/R_{\perp\parallel})\sin^2\psi\).

The fracture plane is found by searching \(\theta\) in 1° steps and refining the maximum
with a parabola. In plane stress the search reproduces Puck's closed-form modes A, B and C
(tested to 2e-4). Mode C is consistent only when \(p^-_{\perp\parallel}/R_{\perp\parallel} =
p_{\perp\perp}/R^A_{\perp\perp}\), Puck's own coupling. Under pure transverse compression
the fracture angle is \(\arccos\sqrt{1/(2(1 + p_{\perp\perp}))}\), about 51° for
\(p_{\perp\perp} = 0.25\).

## In a job

Add `ply_strength.txt`:

```
[Ply_Strength]
[material_id] [Xt] [Xc] [Yt] [Yc] [S12] [p_TL_t] [p_TL_c] [p_TT]
1 1500e6 1200e6 50e6 200e6 70e6 0.35 0.30 0.25
```

- The first six columns are required. Any of `[S23] [Zt] [Zc] [F_star] [p_TL_t] [p_TL_c]
  [p_TT]` may follow, in any order, and `-` keeps a default.
- Material ids are those of `laminate.txt` (`[Ply_Material]`) and `section_mesh_2d.txt`
  (`[Solid_Material]`).
- Materials without a row (cores, for example) are not checked.

A static job (linear, or nonlinear with corotational elements) then writes, in
`secondary_results/ply_failure/`:

- `ply_failure.csv`: one row per ply point per Gauss point. It gives the element, Gauss
  point, `xi` and global `x`; the wall, `s`, ply and `n` (thin-walled) or the section
  element (2D); `y` and `z`; the material; the six ply stresses; each criterion's exposure
  and mode; Puck's fracture angle; and the largest exposure, its margin and the governing
  criterion.
- `ply_failure_summary.csv`: each element's worst point per criterion.

### Where the stresses come from

The Gauss points are those of each element's stiffness integration, read from the element
objects cached at the formulation step (`gauss_data`).

- **Thin-walled sections:** the bottom and top of every ply, at 5 points along each wall (2
  on sections of more than 64 walls), as `wall_stresses.csv`. The transverse shears are
  turned from wall to ply axes (\(\tau_{13}, \tau_{23}\) from \(\tau_{xn}, \tau_{sn}\)).
- **2D sections:** every quadrature point of the section mesh, in material axes.

See [examples/thin_walled/ply_failure](../../examples/thin_walled/ply_failure/README.md).
On a composite box spar, the thin-walled and 2D sections give the same first-ply failure
within 2 % for maximum stress, Hashin and Puck. Tsai-Wu differs by 10 %, because its worst
point is the cap–web corner.

## Under section forces from elsewhere, with an envelope

When the section forces come from another tool (an aeroelastic run, a set of design load
cases), add `section_forces.txt` beside `ply_strength.txt`:

```
[Section_Forces]
[load_case] [element_id] [x] [frame] [f1] [f2] [f3] [f4] [f5] [f6]
dlc1.1_12ms  0  0.25  engine   2.0e3 -1.07e5 1.88e4 3.0e3 2.0e4 5.0e2
dlc6.1_50yr  0  0.25  beamdyn  ...
```

- `element_id` picks the section (the element's, from `element_section.txt`); `x` labels
  the station in the output.
- `frame` is that of the six forces: `engine` is 1d-fem's right-handed
  `[N, M_y, M_z, V_y, V_z, T]`; `public` and `beamdyn` are anisotropic-section-model's and
  BeamDyn's `[F_x, F_y, N, M_x, M_y, M_z]` ([frames](../conventions/SECTION_STIFFNESS_ORDER.md)).
  The forces are about the beam axis through the element nodes.
- 1d-fem's own `tertiary_results/.../section_forces_elem_*.csv` writes `M_y` with the
  element sign: negate it to use those rows here.

Every job with the file (any analysis type: the beam solve plays no part) writes, in
`secondary_results/section_force_checks/`:

- `section_force_checks.csv`: for every row and criterion (and `exposure`, the largest of
  them), the most exposed point: exposure, margin, mode, wall, ply, `n`, `(y, z)`, 2D
  element and material;
- `section_force_envelope.csv`: for every element, station and criterion, the largest
  exposure over the load cases, the load case that sets it and where.

The plies are checked by the section provider's `SectionStrength` at the same points as
`ply_failure.csv` (5 points along each wall, the bottom and top of every ply; every
quadrature point of a 2D section), so a solved job's own section forces given back here
reproduce its `ply_failure.csv` exactly. Each section is solved once for six unit forces,
and a load case then costs a matrix product and the criteria. Table sections
(`section_matrix.txt`) carry no plies and are refused.
