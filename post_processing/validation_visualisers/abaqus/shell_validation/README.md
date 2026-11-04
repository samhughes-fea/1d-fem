# Shell validation of the second-order and corotational beams

Abaqus shell models of five thin-walled beam problems, and the comparison of each with
the 1d-fem beam that should reproduce it. The target is 3 % on the compared quantity.

| case | shell analysis | 1d-fem beam | compared |
|---|---|---|---|
| `channel_ft` | `*BUCKLE`, channel column, fork ends, compression at the tension centre | `LinearBenscoterBeamElement3D` buckling | lowest flexural and flexural-torsional loads |
| `ibeam_ltb_buckle` | `*BUCKLE`, I-beam, fork ends, uniform major-axis moment | `LinearBenscoterBeamElement3D` buckling | lateral-torsional buckling moment |
| `ibeam_ltb_path` | NLGEOM static to 0.95 `M_cr0`, with a lateral force `0.01 * 4M/L` on the compression flange at midspan | `CorotationalBenscoterBeamElement3D` | midspan lateral deflection and twist along the path |
| `ibeam_large_twist` | NLGEOM static, cantilever, root clamped (warping restrained), tip torque to 6 kN m | `CorotationalBenscoterBeamElement3D` | tip twist at 1, 2, 4 and 6 kN m |
| `composite_box` | NLGEOM static, bend-twist coupled CFRP box cantilever, tip flapwise force (linear tip deflection 20 % of the span) | `CorotationalTimoshenkoBeamElement3D` | tip shortening, deflection and twist |

The sections, materials and supports are in `cases.py`. The channel and the I-beams are the
sections of `examples/thin_walled/channel_buckling` and `examples/thin_walled/large_twist`.

## Running it

From the repository root, with Abaqus on the PATH (or `ABAQUS_CAE_ROOT` set, as for the
other Abaqus validation scripts):

```
python -m post_processing.validation_visualisers.abaqus.shell_validation.run --solver abaqus --cpus 4
```

This writes each input deck, runs it (`abaqus job=<case> interactive`), extracts the
station displacements from the ODB (`abaqus python extract_odb.py`) and writes
`summary_abaqus.json`, `report_abaqus.md` and one plot per nonlinear case
(`<case>_abaqus.png`) in this folder. The decks and Abaqus files are under
`post_processing/validation_visualisers/abaqus_results/shell_validation/abaqus/<case>/`.

To run the decks somewhere else, write them with `--generate-only`, run them, put each
`<case>.csv` from `extract_odb.py` next to its deck and compare with `--compare-only`.
`--cases` picks cases and `--fine` scales the shell mesh (1 is the validation mesh; 2
halves the element size, to check convergence).

The 1d-fem results are cached in `beam_results.json`; `--rerun-beam` recomputes them.

## How the shell model matches the beam

- **Walls.** Every wall mid-line becomes a strip of S4R shells, and walls that meet share
  their nodes, as the section analysis assumes. Element nodes run along `x` then along the
  wall direction `s`, so the shell normal is `x cross s`, the section analysis's wall
  normal. Each ply has its own orientation with the fibre at its angle from `x` toward `s`,
  and the plies are listed from the bottom of the wall, so a laminate means the same thing
  in both models.
- **Fork ends.** The end section is tied to a reference node through DOFs 2 and 3 only
  (`*Coupling`, `*Kinematic`, `2, 3`): rigid in its plane, free to warp. The reference node
  holds `v`, `w` and the twist.
- **Clamped root.** Every node of the root section is fixed, which also restrains warping.
- **Tip torque.** The same in-plane coupling, with the torque on the reference node.
- **Rigid tip (box).** A full kinematic coupling. The box is closed, so the restrained tip
  warping is a local effect.
- **End forces and moments.** Axial forces and bending moments are nodal forces with the
  linear stress distribution of a plane section, weighted by each node's share of the axial
  membrane stiffness, so they do not restrain warping either.
- **Beam displacements of the shell.** At each output station the section's rigid in-plane
  motion is fitted to its nodes (a finite rotation for the nonlinear runs, an infinitesimal
  one for the buckling modes). The fit's residual tells a global buckling mode from a local
  plate mode, and only global modes are compared.

## CalculiX cross-check

The same decks run in CalculiX with `--solver calculix` (quadratic S8R shells; CalculiX
takes composite sections on those only, and a rigid end becomes a `*RIGID BODY`).
CalculiX expands shells into solids and ties them at the wall junctions, and with NLGEOM
those ties stop its Newton iterations from converging. So for CalculiX the nonlinear cases
run linearly (`--linear`, compared with the linear beam at load level 1). That checks the
decks' geometry, laminates, supports and loads, but not the nonlinear response, which
needs Abaqus. `report_calculix_linear.md` has the results.
