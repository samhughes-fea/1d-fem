# Eigen pyAbaqus reference contract

This document defines the initial pyAbaqus external-reference contract for eigen validation using the shared workflow under [`post_processing/validation_visualisers/abaqus/`](../../post_processing/validation_visualisers/abaqus/).

## Objective

Promote eigen validation from repository-pinned benchmark checks to canonical external-reference comparison against Abaqus frequencies and selected mode-shape data.

This contract complements the pinned benchmark posture in [`EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md`](EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md).

## Canonical benchmark job

Initial pinned job:

- [`jobs/job_benchmark_eigen_cantilever/`](../../jobs/job_benchmark_eigen_cantilever/)

## Script-generation contract

The shared Abaqus generator in [`job_to_abaqus_script.py`](../../post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py) must accept canonical eigen jobs and generate a CAE script that:

- creates the modal or frequency extraction step required for eigen comparison
- applies the prescribed BCs for the benchmark job
- exports the Abaqus reference artifacts defined below

## Abaqus reference artifact contract

Raw Abaqus reference artifacts should be written under:

- `post_processing/validation_visualisers/abaqus_results/<job_name>/`

Required files:

- `eigen_frequencies.csv`
- `mode_shapes.csv`

Minimum `eigen_frequencies.csv` columns:

- `mode_index`
- `frequency_hz`

Minimum `mode_shapes.csv` columns should include:

- `mode_index`
- node or DOF identifiers
- exported displacement components for selected coordinates

## FEM artifact contract

The FEM side should provide modal results sufficient to compare against the Abaqus frequency and selected mode-shape exports for the same benchmark job.

At minimum, the family comparator should be able to read:

- mode frequencies
- selected mode-shape coordinates or DOF values

## Required comparison quantities

At minimum, the eigen family should compare:

- frequencies for the first selected modes
- selected normalized mode-shape coordinates for the same modes

Later phases may add:

- modal assurance style metrics
- broader mode-shape field comparisons

## Family-local ownership

The eigen family layer under [`post_processing/validation_visualisers/eigen/`](../../post_processing/validation_visualisers/eigen/) should own:

- FEM-vs-Abaqus comparison logic
- mode tables and plots
- family summary outputs
- tolerance policy for the pinned benchmark set

## Scope note

This phase defines the external-reference contract only. It does not yet implement the full eigen FEM-vs-Abaqus comparator.
