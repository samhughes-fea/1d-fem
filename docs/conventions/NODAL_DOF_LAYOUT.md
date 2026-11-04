# Nodal DOF layout

Source: `pre_processing/element_library/dof_layout.py`.

## Mesh stride

Every node carries the six standard DOFs and then the same number of extra DOFs:

| local index | name | meaning |
|---|---|---|
| 0–5 | `UX UY UZ RX RY RZ` | translations and rotations |
| 6 | `Q1` (alias `W`, `CHI`) | warping amplitude (Vlasov χ or Benscoter ψ) |
| 7 … | `Q2 … Qn` | further extra DOFs (cross-section mode amplitudes) |

The number of extra DOFs is the largest any element row asks for: one for
`[warping] = 1`, or `element_dictionary["extra_dofs"][row]` when that key is set.
`mesh_dof_per_node(element_dictionary)` returns the stride, and every runner uses it
for `total_dof = n_nodes × stride` and `global_dof = node_id × stride + local_index`.

## Element DOF order

An element's `K_e`, `F_e` and `U_e` are ordered as

```text
[node 1: 0..5, node 2: 0..5, node 1: extra, node 2: extra]
```

`Element1DBase.assemble_global_dof_indices` maps that order onto the mesh stride. An
element with fewer extra DOFs than the stride (a truss or plain beam on a warping mesh)
uses only the first ones, and DOFs that no element touches must be restrained.

## Input names

`prescribed_displacement.txt` accepts the standard names, `Q1 … Qn`, and `W` / `CHI` for `Q1`.
A name that does not exist on the mesh stride is skipped with a warning.
