# 2D section meshes and 3D materials

The 2D section solver (Giavotto, as in BECAS; see `solid_section_analysis.md`) works on a
finite-element mesh of the cross-section. This page covers its inputs: the mesh and the
materials. Code:
`pre_processing/section_integration/solid/`.

## Elements

`T6` (6-node triangle) and `Q8` (8-node quadrilateral) quadratic elements, with Gmsh's
node order: corners counter-clockwise first, then the mid-side nodes. They integrate with
a degree-5 rule: 7 points on `T6`, 3 × 3 Gauss on `Q8`. Inverted or degenerate elements
are refused.

## Materials and orientation

Each element carries a `SolidMaterial`: nine orthotropic engineering constants in its own
axes. `SolidMaterial.from_ply` builds one from a thin-walled `PlyMaterial`, with
`E3 = E2`, `nu13 = nu12` and `nu23 = E2 / (2 G23) − 1`. Each element also carries two
angles:

- `alpha`: the wall direction `s` in the section, from `y` toward `z`;
- `theta`: the ply angle, from the beam axis `x` toward `s`, as in `laminate.txt`.

The material axes are then `e1 = cos θ x + sin θ s`, `e2 = −sin θ x + cos θ s` and
`e3 = n = x × s`. The stiffness in section axes is `C = T C₁₂₃ Tᵀ`, with `T` the Bond
matrix, in the Voigt order `[xx, yy, zz, yz, xz, xy]`.

The test `test_plane_stress_of_the_3d_material_is_the_ply_stiffness` condenses the
through-thickness stresses out of `C` in wall axes. It recovers the thin-walled ply's
`Q̄` for any `θ` and `α`, so the two section models put the fibres in the same place.

## Meshes

- `SectionMesh.rectangle(b, h, n_y, n_z, material)`: a solid rectangle.
- `SectionMesh.from_thin_walled(section, laminates)`: each wall of a thin-walled section
  becomes a strip of `Q8` elements, with each ply a row of elements through the thickness.
  - Each strip is its wall's mid-line length long, so the junctions count twice, as in the
    mid-line model.
  - Strips are joined at the section points: at each point the first wall listed is the
    master. The end-edge nodes of the later walls are tied into the master strip, as the
    interpolation of its nearest element (`SectionMesh.ties`).
- `SectionMesh.from_gmsh(path, groups, materials)`: a Gmsh 4.1 ASCII mesh. Each physical
  surface maps to `(material, theta, alpha)`, and its `T6` and `Q8` elements are kept. The
  mesh's `(x, y)` coordinates are the section's `(y, z)`.
- `section_mesh_2d.txt` in a job (`parse_section_mesh_2d`):

```
[Section_Node]
[section_id] [node_id] [y] [z]
[Section_Element]
[section_id] [element_id] [type] [nodes...] [material_id] [theta] [alpha]
[Solid_Material]          (optional; other material ids come from laminate.txt)
[material_id] [E1] [E2] [E3] [G12] [G13] [G23] [nu12] [nu13] [nu23] [rho]
```
