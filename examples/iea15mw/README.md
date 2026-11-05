# IEA 15 MW blade stations

Shared input for the blade examples ([blade_beam](../blade_beam/README.md),
[blade_section_2d](../blade_section_2d/README.md)) and the section validation in
`validation/iea15mw_sections`. Nothing here runs on its own.

- The stations themselves are blade input data, so they are not in this public
  repository. `stations.py` reads them from the private reference-blades catalogue
  (`iea15/derived/sections/precomp/precomp_sections.json`: the IEA Wind 15 MW blade as PreComp sees it,
  taken from the turbine's ontology file via WISDEM) through `reference_blades.py`. Clone
  the catalogue to `reference_blades/` at the repository root (git-ignored) or set
  `REFERENCE_BLADES_DIR`. Without it the blade examples, the section validation and their
  tests are skipped; their committed figures and results stay. The examples use the root
  (r = 4.0 m, chord 5.20 m), maximum chord (r = 28.2 m, 5.77 m) and 70 % span
  (r = 84.7 m, 3.27 m) stations.
- `stations.py`: each station as an implicit-geometry-kernel section (`geometry`), its
  wall laminates (`laminates`) and ply materials (`ply_materials`). The module docstring
  states the idealisation every model shares: one skin layup all round, carbon caps between
  the webs, sandwich webs, no trailing-edge reinforcement.
- `segment.py`: the outboard 20.2 m (r = 64.5 m to 84.7 m) as one tapered, twisted loft
  of the 70 % span section, or a prismatic one (`geometry(prismatic=True)`).

The geometry functions need implicit-geometry-kernel
(`pip install -r requirements-geometry.txt`, Python 3.11 or later).
