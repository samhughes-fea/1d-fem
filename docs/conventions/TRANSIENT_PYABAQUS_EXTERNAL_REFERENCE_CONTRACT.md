# Transient pyAbaqus external-reference contract

This document defines the file and script-generation contract for transient external-reference validation using the existing pyAbaqus workflow under [`post_processing/validation_visualisers/abaqus/`](../../post_processing/validation_visualisers/abaqus/).

## Target benchmark

Pinned structural case:

- [`jobs/job_benchmark_transient_cantilever_multidof/`](../../jobs/job_benchmark_transient_cantilever_multidof/)

## Script-generation contract

The transient external-reference workflow must be able to generate an Abaqus CAE script from the pinned job directory using the same mechanism already used for linear static validation.

Minimum contract:

1. `_parse_job(...)` accepts the transient benchmark job.
2. `_generate_script_content(...)` embeds:
   - job name
   - output directory
   - prescribed BCs
   - point-load data
3. the parsed payload exposes the simulation settings path for downstream transient-specific export extensions.

## Reference-file contract

The external-reference wave will compare FEM transient output against Abaqus-exported time-history data at selected DOFs.

Expected reference location pattern:

- `post_processing/validation_visualisers/abaqus_results/<job_name>/`

Expected future transient reference files:

- selected-DOF transient response CSV(s)
- optional metadata file describing time grid and DOF mapping

## Scope note

This phase only establishes the pyAbaqus export and reference-file contract. It does not yet implement the final FEM-vs-Abaqus time-history comparison.
