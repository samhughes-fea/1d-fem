# Abaqus validation directory standard

This document standardizes how Abaqus-based validation is organized across simulation families.

## Shared vs family-local ownership

### Shared Abaqus layer

Keep shared under [`post_processing/validation_visualisers/abaqus/`](../../post_processing/validation_visualisers/abaqus/):

- Abaqus launcher/config
- job parsing
- CAE script generation
- ODB extraction helpers
- common raw Abaqus result storage

### Family-local validation layer

Each simulation family owns its own comparison and reporting logic under:

- [`post_processing/validation_visualisers/static/`](../../post_processing/validation_visualisers/static/)
- [`post_processing/validation_visualisers/eigen/`](../../post_processing/validation_visualisers/eigen/)
- [`post_processing/validation_visualisers/buckling/`](../../post_processing/validation_visualisers/buckling/)
- [`post_processing/validation_visualisers/transient/`](../../post_processing/validation_visualisers/transient/)
- [`post_processing/validation_visualisers/harmonic/`](../../post_processing/validation_visualisers/harmonic/)

## Required family-local subfolders

Each family directory should contain:

- `reference_cases/`
- `plots/`
- `tables/`
- `output/`

and one family README:

- `README_VALIDATION_<FAMILY>.md`

## File ownership rules

### Shared layer owns
- script generation
- Abaqus invocation
- ODB export contracts
- raw `abaqus_results/<job_name>/...`

### Family layer owns
- FEM vs Abaqus comparisons
- family-specific metrics
- plots/tables
- benchmark summaries
- tolerance policies and case documentation

## Naming guidance

- raw Abaqus outputs remain in [`abaqus_results/`](../../post_processing/validation_visualisers/abaqus/config.py)
- family summaries go in each family `output/`
- generated figures go in each family `plots/`
- tabulated benchmark results go in each family `tables/`

## Promotion path

Each family should grow through the benchmark ladder defined in [`BENCHMARK_LADDER_STANDARD.md`](BENCHMARK_LADDER_STANDARD.md):

- smoke
- pinned acceptance
- calibrated repository reference
- external-reference validated
