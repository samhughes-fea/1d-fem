# Linear static pyAbaqus reference contract

This document defines the canonical pyAbaqus external-reference contract for linear static validation using the shared workflow under [`post_processing/validation_visualisers/abaqus/`](../../post_processing/validation_visualisers/abaqus/).

## Objective

Normalize linear static validation into the same contract-driven structure now used by newer simulation families.

The linear static family should remain the baseline reference family for cantilever load-case comparisons while using the shared Abaqus directory and artifact conventions from [`ABAQUS_VALIDATION_DIRECTORY_STANDARD.md`](ABAQUS_VALIDATION_DIRECTORY_STANDARD.md).

## Canonical load families

The family should cover the six canonical cantilever load cases:

1. tip point load
2. midspan point load
3. quarter-point point load
4. UDL
5. triangular distributed load
6. parabolic distributed load

## Canonical benchmark jobs

The pinned linear-static suite should be expressed through canonical jobs under [`jobs/`](../../jobs/).

At minimum, each load family should have:

- one pinned benchmark root
- a mesh ladder suitable for FEM-vs-Abaqus comparison

## Script-generation contract

The shared Abaqus generator in [`job_to_abaqus_script.py`](../../post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py) must accept canonical linear static jobs and generate a CAE script that:

- creates the beam geometry and section
- applies the prescribed displacements
- applies point or distributed loads as defined by the job
- runs a linear static Abaqus step
- exports the reference artifacts defined below

## Abaqus reference artifact contract

Raw Abaqus reference artifacts should be written under:

- `post_processing/validation_visualisers/abaqus_results/<job_name>/`

Required files:

- `U_global.csv`
- `rotation_source.txt`

Optional future files:

- selected summary displacement CSV
- deformed-shape sample CSV

## FEM artifact contract

The FEM side should provide artifacts sufficient to compare against the Abaqus exports for the same pinned job.

At minimum, the family comparator should be able to read:

- the FEM displacement output corresponding to `U_global.csv`
- any family-local summary metric used for pinned acceptance

## Required comparison quantities

At minimum, the linear static family should compare:

- nodal displacement parity at exported DOFs
- tip displacement parity for the canonical load case

Later phases may add:

- deformed-shape path parity
- rotation parity where available
- reaction-force parity where practical

## Family-local ownership

The linear static family layer under [`post_processing/validation_visualisers/static/`](../../post_processing/validation_visualisers/static/) should own:

- FEM-vs-Abaqus comparison logic
- family plots and tables
- pass/fail evaluation
- family summary output

## Scope note

This document defines the contract normalization target for linear static pyAbaqus validation. It does not itself implement the full family-local comparator rollout for all six load cases.
