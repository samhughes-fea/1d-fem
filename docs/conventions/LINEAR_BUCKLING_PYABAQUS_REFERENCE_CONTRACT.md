# Linear buckling pyAbaqus reference contract

This document defines the initial pyAbaqus external-reference contract for linear buckling validation using the shared workflow under [`post_processing/validation_visualisers/abaqus/`](../../post_processing/validation_visualisers/abaqus/).

## Objective

Promote linear buckling validation from repository-pinned benchmark checks to canonical external-reference comparison against Abaqus buckling factors and selected buckling mode-shape data.

This contract complements the pinned benchmark posture in [`EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md`](EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md).

## Canonical benchmark job

Initial pinned job:

- [`jobs/job_benchmark_linear_buckling_column/`](../../jobs/job_benchmark_linear_buckling_column/)

## Script-generation contract

The shared Abaqus generator in [`job_to_abaqus_script.py`](../../post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py) must accept canonical linear buckling jobs and generate a CAE script that:

- creates the buckling extraction step required for the benchmark job
- applies the prescribed BCs and reference loading
- exports the Abaqus reference artifacts defined below

## Abaqus reference artifact contract

Raw Abaqus reference artifacts should be written under:

- `post_processing/validation_visualisers/abaqus_results/<job_name>/`

Required files:

- `buckling_load_factors.csv`
- `buckling_mode_shapes.csv`

Minimum `buckling_load_factors.csv` columns:

- `mode_index`
- `load_factor`

Minimum `buckling_mode_shapes.csv` columns should include:

- `mode_index`
- node or DOF identifiers
- exported displacement components for selected coordinates

## FEM artifact contract

The FEM side should provide buckling results sufficient to compare against the Abaqus load-factor and selected mode-shape exports for the same benchmark job.

At minimum, the family comparator should be able to read:

- buckling load factors
- selected mode-shape coordinates or DOF values

## Required comparison quantities

At minimum, the linear buckling family should compare:

- first selected buckling load factors
- selected normalized mode-shape coordinates for those modes

Later phases may add:

- higher-mode comparisons
- broader mode-shape correlation metrics

## Family-local ownership

The buckling family layer under [`post_processing/validation_visualisers/buckling/`](../../post_processing/validation_visualisers/buckling/) should own:

- FEM-vs-Abaqus comparison logic
- buckling summary tables and plots
- tolerance policy for the pinned benchmark set
- family summary outputs

## Scope note

This phase defines the external-reference contract only. It does not yet implement the full linear-buckling FEM-vs-Abaqus comparator.
