# API reference

Generated from the modules' docstrings by `python scripts/build_api_reference.py`; do not edit by hand. A module is listed when it declares `__all__` or a `Public API` section. The theory behind each part is in [element_library](../element_library/) and the job file formats in [conventions](../conventions/).

| Page | Covers |
|---|---|
| [Section analysis](sections.md) | The section provider seam and the section input and recovery that 1d-fem keeps. The cross-section analysis itself (thin-walled and 2D sections, GBT modes, stress recovery, ply failure, layup sensitivities) comes from the section provider, anisotropic-section-model by default. |
| [Elements](elements.md) | Beam, bar and truss elements and their shared utilities. |
| [Job input](job_input.md) | Parsers for the job directory's text files, meshes and loads. |
| [Solvers](solvers.md) | Static, nonlinear, eigen, buckling, harmonic, transient and spectral solvers and their runners. |
| [Post-processing](post_processing.md) | Result files written after a solve: stresses, failure margins, sensitivities and field sampling for plots. |
| [Running jobs](workflow.md) | Running a job directory end to end. |
