# Harmonic pyAbaqus reference contract

This document defines the initial pyAbaqus external-reference contract for harmonic validation using the shared workflow under [`post_processing/validation_visualisers/abaqus/`](../../post_processing/validation_visualisers/abaqus/).

## Objective

Extend harmonic validation beyond analytical-only checks so canonical harmonic benchmark jobs can be compared against Abaqus-exported frequency-response data.

This contract complements the analytical benchmark note in [`HARMONIC_SDOF_REFERENCE_VALIDATION.md`](HARMONIC_SDOF_REFERENCE_VALIDATION.md).

## Canonical benchmark job

The harmonic family should begin with one pinned benchmark job under [`jobs/`](../../jobs/) suitable for frequency-domain comparison.

The initial benchmark should define:

- forcing frequency range
- number of frequency points
- damping settings
- selected response DOFs for comparison

## Script-generation contract

The shared Abaqus generator in [`job_to_abaqus_script.py`](../../post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py) must accept canonical harmonic jobs and generate a CAE script that:

- creates the harmonic analysis step for the requested frequency range
- applies the prescribed BCs and harmonic loading
- writes the raw Abaqus reference artifacts described below

## Abaqus reference artifact contract

Raw Abaqus reference artifacts should be written under:

- `post_processing/validation_visualisers/abaqus_results/<job_name>/`

Required files:

- `frequency_response.csv`

Required columns should include at minimum:

- `frequency_hz`
- selected displacement response values for the chosen DOFs

If complex response components are exported explicitly, later phases may also require:

- `component_name`
- real component values
- imaginary component values

## FEM artifact contract

The FEM side should provide harmonic response artifacts sufficient to compare against the Abaqus frequency-response export for the same benchmark job.

At minimum, the family comparator should be able to read:

- frequency samples
- selected DOF harmonic responses

## Required comparison quantities

At minimum, the harmonic family should compare:

- response amplitude or complex component response versus frequency at selected DOFs

Later phases may add:

- resonance peak location comparison
- phase comparison
- real and imaginary component parity

## Family-local ownership

The harmonic family layer under [`post_processing/validation_visualisers/harmonic/`](../../post_processing/validation_visualisers/harmonic/) should own:

- FEM-vs-Abaqus comparison logic
- plots and tables
- tolerance policy for canonical benchmark cases
- summary outputs

## Scope note

This phase defines the harmonic pyAbaqus file contract and family ownership only. It does not yet implement the full harmonic FEM-vs-Abaqus comparator.
