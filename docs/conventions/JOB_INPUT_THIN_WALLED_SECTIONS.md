# Job input: thin-walled anisotropic sections

Three optional job files describe anisotropic thin-walled sections from their walls and plies and assign them to elements. They must be given together. A job with none of them runs exactly as before. `section.txt` is still required for every element and still supplies the properties of elements that are not assigned a thin-walled section.

## `laminate.txt`

```text
[Ply_Material]
[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]
1  140e9 10e9 5e9 0.3 5e9 3.5e9 1600

[Laminate]
[laminate_id] [ply] [material_id] [angle_deg] [thickness]
1  1  1   0.0  0.125e-3
1  2  1  45.0  0.125e-3
```

- `[rho]` is optional (omit the column on every row, or give it on every row).
- `E1`, `E2`, `G12`, `G13`, `G23` must be positive and `nu12**2 * E2 / E1 < 1`.
- Plies are numbered `1..n` within each laminate. Ply 1 is the bottom ply (at \(n=-h/2\)); rows may come in any order.
- `angle_deg` is measured from the beam axis \(x\) toward the wall direction \(s\).
- Each laminate gets transverse shear factor 5/6.

Parser: `pre_processing.parsing.laminate_parser.parse_laminate_file` returns `{laminate_id: WallLaminate}`.

## `section_geometry.txt`

```text
[Section_Point]
[section_id] [point_id] [y] [z]
1  1  0.0  0.0
1  2  0.4  0.0
1  3  0.4  0.2
1  4  0.0  0.2

[Section_Wall]
[section_id] [wall_id] [point_start] [point_end] [laminate_id]
1  1  1  2  1
1  2  2  3  1
1  3  3  4  1
1  4  4  1  1
```

- Points are mid-line points in the section plane \((y, z)\). Walls are straight segments between them.
- A wall runs from `point_start` to `point_end`; that sets its \(s\) direction and so the side its plies stack toward (\(\mathbf{e}_n = \mathbf{e}_x\times\mathbf{e}_s\)). On a counter-clockwise cell the bottom ply is on the outside.
- Rejected: unknown point or laminate ids, zero-length walls, two walls on the same pair of points, points on no wall, duplicate ids.
- Closed cells are counted from the wall graph (`walls - points + components`).

Parser: `pre_processing.parsing.section_geometry_parser.parse_section_geometry_file` returns `{section_id: ThinWalledSection}`.

## `element_section.txt`

```text
[Element_Section]
[element_id] [section_id]
0  1
1  1
```

- Each listed element uses the thin-walled section `section_id`. Each section is analysed once, however many elements use it.
- Section coordinates \((y, z)\) are measured from the beam reference line through the element nodes, in the global \(y\) and \(z\) directions. The linear elements have no element frames yet, so the beam axis is global \(x\).
- Unknown element or section ids are errors.
- Element types that do not yet build their stiffness from a `SectionAnalysis` raise `NotImplementedError` when assigned a section, rather than silently using `section.txt`.
- `LinearGBTBeamElement3D` rows carry the section's warping modes. Optional columns set how many distortional modes (`[modes]`) and extra shear-lag warping modes (`[shear_modes]`) they carry. Each node then has `n_sv + shear_modes + 2 × modes` extra DOFs `Q1 … Qn`. Here `n_sv` is the number of Saint-Venant warping modes of the section (3 for a rectangular box). The warping amplitudes come first (`Q1 … Q(n_sv + shear_modes)`), then each mode's amplitude and slope. Warping restraint is all warping amplitudes zero. A diaphragm (rib) is the mode amplitudes zero; also fixing the slopes stops the mode warping. See [../element_library/gbt_cross_section_modes.md](../element_library/gbt_cross_section_modes.md).

```text
[Element_Section]
[element_id] [section_id] [modes] [shear_modes]
0  1  2  1
1  1  2  1
```

## `[Section_Mesh]` (optional, in `section_geometry.txt`)

```text
[Section_Mesh]
[section_id] [n_sub]
1  4
```

Sets how many cubic elements each wall of a section is split into (default 1). More elements give more warping and distortional modes.

## `wall_load.txt` (optional)

```text
[Wall_Load]
[element_id] [wall_id] [p_x] [p_s] [p_n]
4  3  0.0  0.0  -2000.0
```

A uniform surface load (force per area) on a wall over an element: `p_x` along the beam axis, `p_s` along the wall direction and `p_n` along the wall normal. Rows for the same element and wall add up. GBT elements get the load on all their DOFs, so an eccentric load also drives distortion. Other elements get the resultant forces and moments (and the Benscoter warping load).

The parsed result reaches every element as `element.section_analysis` (a `SectionAnalysis`, see [../element_library/section_analysis.md](../element_library/section_analysis.md)).

Theory: [../element_library/thin_walled_section_walls.md](../element_library/thin_walled_section_walls.md).
