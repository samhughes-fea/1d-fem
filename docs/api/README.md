# API reference

Generated from the modules' docstrings by `python scripts/build_api_reference.py`; do not edit by hand. A module is listed when it declares `__all__` or a `Public API` section. The theory behind each part is in [element_library](../element_library/) and the job file formats in [conventions](../conventions/).

| Page | Covers |
|---|---|
| [Section analysis](sections.md) | Thin-walled and 2D (solid) cross-section analysis: stiffness, mass, warping, GBT modes, stress recovery, ply failure and layup sensitivities. |
| [Elements](elements.md) | Beam, bar and truss elements and their shared utilities. |
| [Job input](job_input.md) | Parsers for the job directory's text files, meshes and loads. |
| [Solvers](solvers.md) | Static, nonlinear, eigen, buckling, harmonic, transient and spectral solvers and their runners. |
| [Post-processing](post_processing.md) | Result files written after a solve: stresses, failure margins, sensitivities and field sampling for plots. |
| [Interfaces](interfaces.md) | Adapters from other repositories' data (a blade's geometry, duck typed) to sections and jobs. The core does not depend on them. |
| [Running jobs](workflow.md) | Running a job directory end to end. |
