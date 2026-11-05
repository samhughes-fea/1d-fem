# 1d-fem

Beam finite elements in 3D for slender composite structures: blades, wings and spars.
A job is a directory of text files; 1d-fem solves it and writes the results, from nodal
displacements down to ply stresses and failure margins at every Gauss point.

- **Elements:** Euler-Bernoulli, Timoshenko, Levinson and Reddy (third-order shear),
  Benscoter (restrained warping) and GBT (cross-section distortion) beams, plus bars and
  trusses; total-Lagrangian, corotational and geometrically exact beams for large
  rotations. Tapered and pretwisted elements. [Element types](docs/api/elements.md).
- **Sections:** anisotropic thin-walled sections from walls and ply layups (multi-cell,
  open or closed), and 2D layered finite-element sections for thick walls. Each gives the
  coupled 6 x 6 stiffness, mass, warping and GBT modes, and recovers wall, ply and
  interlaminar stresses. [Thin-walled](docs/element_library/thin_walled_section_analysis.md),
  [2D](docs/element_library/solid_section_analysis.md).
- **Analyses:** linear and nonlinear static (Newton, arc length), eigen, linear buckling
  (global, local and distortional with GBT), P-delta, harmonic and transient. Rotating
  structures (centrifugal stiffening, spin softening, gyroscopic terms) and follower loads.
- **Composite design:** ply failure criteria (maximum stress, Tsai-Wu, Hashin, Puck) and
  exact sensitivities of displacements and frequencies to ply angles and thicknesses.
- **Blades:** a beam model of a blade straight from an implicit-geometry-kernel loft.

## Install

Python 3.10 or later. 1d-fem is the `onedfem` package (`src/onedfem`). There is no release
on PyPI yet; install it from GitHub:

```bash
pip install "onedfem @ git+https://github.com/samhughes-fea/1d-fem"
pip install "onedfem[geometry] @ git+https://github.com/samhughes-fea/1d-fem"   # Python 3.11+: blade geometry
```

or, to run the examples, validation and tests, from a clone:

```bash
pip install -e ".[test]"
pip install -e ".[geometry]"   # optional, Python 3.11+
```

The repository keeps the code in `src/onedfem` and, outside the package, the examples, the
validation and verification cases with their reference data (`validation/`,
`verification/`), the tracked result trees (`results/`) and the job directories (`jobs/`).

## Running a job

```python
from onedfem.workflow_orchestrator.run_job import process_job, setup_job_results_directory

results = setup_job_results_directory("my_job")    # under results/
process_job("path/to/my_job", results, {}, {})
```

`onedfem-run-jobs` (or `python -m onedfem.workflow_orchestrator.run_job`) runs every
`jobs/job_*` directory of the project root: the repository in a clone, otherwise the working
directory or `ONEDFEM_PROJECT_ROOT`. A job
directory holds:

| File | Content |
|---|---|
| `grid.txt` | Node coordinates |
| `element.txt` | Element connectivity and type |
| `material.txt`, `section.txt` | Isotropic material and section properties per element |
| `simulation_settings.txt` | Analysis type and solver options ([taxonomy](docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md)) |
| `point_load.txt`, `distributed_load.txt` | Loads |
| `prescribed_displacement.txt` | Supports |
| `laminate.txt`, `section_geometry.txt`, `element_section.txt` | Optional: composite thin-walled or 2D sections ([format](docs/conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md)) |
| `section_taper.txt`, `precurvature.txt` | Optional: taper and pretwist |
| `follower_load.txt` | Optional: follower loads ([format](docs/conventions/JOB_INPUT_FOLLOWER_LOADS.md)) |
| `ply_strength.txt` | Optional: ply failure margins ([ply_failure](docs/element_library/ply_failure.md)) |
| `layup_sensitivity.txt` | Optional: layup sensitivities ([layup_sensitivity](docs/element_library/layup_sensitivity.md)) |

The examples write their jobs with `examples/thin_walled/common.py`, which is the quickest
way to see each file's format.

## Examples

| Example | What it shows |
|---|---|
| [thin_walled](examples/thin_walled/README.md) | Composite box spars, channels, I-beams and a three-cell blade section: beam theories against each other, buckling, large twist, ply failure, GBT local buckling and layup sensitivities |
| [blade_beam](examples/blade_beam/README.md) | A beam of the outboard IEA 15 MW blade built from a loft, against a shell model |
| [blade_section_2d](examples/blade_section_2d/README.md) | 2D layered section of an IEA 15 MW blade section: ply and interlaminar stresses |
| [rotating_beam](examples/rotating_beam/README.md) | Rotating beams: centrifugal stiffening, spin softening and the Campbell diagram |
| [follower_loads](examples/follower_loads/README.md) | Follower forces and moments at large rotation |
| [iea15mw](examples/iea15mw/README.md) | The IEA 15 MW blade stations shared by the blade examples |

Every figure is drawn from the FEA data itself: element shape-function interpolation along
the beam, Gauss-point values and their shape-function projections on sections.

## Related repositories

| Repository | Owns | Relation to 1d-fem |
|---|---|---|
| implicit-geometry-kernel | blade and section geometry (airfoils, lofts, `section_at(z)`) | optional; `src/onedfem` reads its geometry duck typed and never imports it |
| laminate-constitutive-model | materials, plies, laminates, ABD, failure criteria | not used yet (1d-fem has its own laminate code) |
| shell-fem | shell and solid FE, blade shell models (linear static) | shell reference in examples and validation only |
| anisotropic-section-model | blade 6×6 / 7×7 section properties, small exact beam references | cross-checks in examples and validation only |
| 1d-fem | beam solver for all analysis types, with its own thin-walled, 2D and GBT section analysis | depends on none of the others |

`tests/test_repo_boundaries.py` fails if library code imports any of them. Tests that use
them skip when they are not installed.

## Documentation

- [docs/api](docs/api/README.md): API reference, generated from the docstrings.
- [docs/element_library](docs/element_library/): theory of each element, section model and
  feature, with its checks.
- [docs/conventions](docs/conventions/README.md): job file formats, DOF layout and coding
  standards.
- [docs/CHANGELOG.md](docs/CHANGELOG.md).

## Testing

```bash
pip install -e ".[test]"
python -m pytest tests -q -n auto
```

CI ([.github/workflows/pytest.yml](.github/workflows/pytest.yml)) runs on Python 3.10 and
3.11 with `FEM_FORMULATION_CACHE_STRICT_SHAPE=1`, so a missing Gauss `shape_functions` cache
fails fast. Set it locally to match:

```bash
export FEM_FORMULATION_CACHE_STRICT_SHAPE=1
```

Smoke jobs with warping are under [`jobs/`](jobs/README_JOBS.md)
(`job_smoke_nl_eb_warp`, `job_smoke_nl_ts_warp`).


## License

MIT, see [LICENSE](LICENSE).
