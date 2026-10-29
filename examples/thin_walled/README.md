# Thin-walled composite section examples

Worked examples that go from a laminated thin-walled section to beam
deflections and then back down to wall, ply and interlaminar stresses.

| Example | What it shows |
|---|---|
| [box_spar](box_spar/README.md) | CFRP box spar: four beam theories against the exact Saint-Venant beam, Bredt torsion, shear lag at a restrained root, ply and interlaminar stresses |
| [channel_buckling](channel_buckling/README.md) | Channel column with offset shear and tension centres: flexural-torsional buckling and P-delta through the job runners |
| [three_cell_blade](three_cell_blade/README.md) | Three-cell blade section with spar caps and two webs: offset shear centre, multi-cell Bredt-Batho shear flows, cap and web stresses |

Run each from the repository root, for example:

```
python examples/thin_walled/box_spar/run.py
```

`common.py` holds the shared pieces: it writes the 1d-fem job files of a
cantilever (`grid.txt`, `element.txt`, `laminate.txt`,
`section_geometry.txt`, `element_section.txt`, loads and supports), runs it
with `process_job`, and reads back `U_global.csv` and `wall_stresses.csv`. It
also has `bredt_multicell`, the classical multi-cell Bredt-Batho torsion
solution used as an independent check.

Method: [docs/element_library/thin_walled_stress_recovery.md](../../docs/element_library/thin_walled_stress_recovery.md).
