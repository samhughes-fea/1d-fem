# pyAbaqus validation suite implementation plan

## Objective

Build a canonical pyAbaqus external-reference validation suite across the major simulation families, using a shared Abaqus generation and extraction layer under [`post_processing/validation_visualisers/abaqus/`](post_processing/validation_visualisers/abaqus/) and family-local comparison/reporting layers under [`post_processing/validation_visualisers/`](post_processing/validation_visualisers/).

This plan is intended as a future reference for implementation sequencing, scope control, and benchmark-family alignment.

## Target simulation families

The suite should cover the following families:

1. linear static
2. nonlinear static
3. transient
4. harmonic
5. eigen
6. linear buckling
7. nonlinear buckling

These families do not all start at the same maturity level, so implementation should proceed in waves.

## Current maturity snapshot

### Most mature

- nonlinear static has the strongest external-reference framing in [`docs/conventions/NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md`](docs/conventions/NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md)
- transient already has an explicit script-generation and file-contract document in [`docs/conventions/TRANSIENT_PYABAQUS_EXTERNAL_REFERENCE_CONTRACT.md`](docs/conventions/TRANSIENT_PYABAQUS_EXTERNAL_REFERENCE_CONTRACT.md)

### Partially mature

- eigen and linear buckling have pinned repository benchmark definitions in [`docs/conventions/EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md`](docs/conventions/EIGEN_AND_LINEAR_BUCKLING_BENCHMARKS.md), but not yet a full pyAbaqus contract
- harmonic has analytical validation in [`docs/conventions/HARMONIC_SDOF_REFERENCE_VALIDATION.md`](docs/conventions/HARMONIC_SDOF_REFERENCE_VALIDATION.md), but not a pyAbaqus external-reference contract

### Least mature for external-reference validation

- linear static likely has legacy workflow support but lacks a newly normalized family contract in the current conventions set
- nonlinear buckling has an MVP continuation contract in [`docs/conventions/NONLINEAR_BUCKLING_CONTINUATION.md`](docs/conventions/NONLINEAR_BUCKLING_CONTINUATION.md) and benchmark posture in [`docs/conventions/NONLINEAR_BUCKLING_BENCHMARKS.md`](docs/conventions/NONLINEAR_BUCKLING_BENCHMARKS.md), but external-reference comparison remains a later-stage problem

## Design principles

### 1. Shared Abaqus layer, family-owned comparisons

The shared layer under [`post_processing/validation_visualisers/abaqus/`](post_processing/validation_visualisers/abaqus/) should own:

- job parsing
- simulation-type dispatch
- CAE script generation
- Abaqus execution hooks
- ODB extraction helpers
- canonical raw result layout in `abaqus_results/<job_name>/`

Each family layer under [`post_processing/validation_visualisers/`](post_processing/validation_visualisers/) should own:

- comparison metrics
- tables and plots
- family summaries
- tolerance policies
- pass/fail evaluation

This division is consistent with [`docs/conventions/ABAQUS_VALIDATION_DIRECTORY_STANDARD.md`](docs/conventions/ABAQUS_VALIDATION_DIRECTORY_STANDARD.md).

### 2. Canonical benchmark job per family before broad rollout

Each family should first lock one canonical job and one benchmark contract before expanding to a full mesh ladder or full load family matrix.

### 3. Contract first, implementation second

Before adding extractor or comparator code, define the required Abaqus-side and FEM-side artifacts in conventions docs.

### 4. Promote through the benchmark ladder

Each family should progress through the benchmark ladder defined in [`docs/conventions/BENCHMARK_LADDER_STANDARD.md`](docs/conventions/BENCHMARK_LADDER_STANDARD.md):

- smoke
- pinned acceptance
- calibrated repository reference
- external-reference validated

## Required suite-wide contract matrix

Each simulation family should eventually define the following fields:

| Field | Meaning |
|---|---|
| canonical benchmark job | pinned job directory for the family |
| Abaqus generation contract | what the generated pyAbaqus script must create or export |
| Abaqus extracted artifacts | CSV or metadata files written from ODB |
| FEM artifacts | files written by the in-repo runner for comparison |
| comparison quantities | scalar, path, field, or history metrics |
| tolerance policy | exact comparison tolerance ownership |
| comparator entrypoint | family-local script or runner |
| family summary artifact | final pass/fail and summary report location |

## Family-by-family implementation plan

### Wave 1 — normalize the highest-value external-reference families

This wave should establish a coherent cross-family suite shape using the families that are either already partially scaffolded or easiest to compare.

#### 1. Linear static

Linear static should be formalized as a first-class pyAbaqus validation family.

Planned deliverables:

- add [`docs/conventions/LINEAR_STATIC_PYABAQUS_REFERENCE_CONTRACT.md`](docs/conventions/LINEAR_STATIC_PYABAQUS_REFERENCE_CONTRACT.md)
- optionally add [`docs/conventions/LINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md`](docs/conventions/LINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md)
- define canonical load-family set:
  - tip point load
  - midspan point load
  - quarter-point point load
  - UDL
  - triangular distributed load
  - parabolic distributed load
- define required Abaqus exports:
  - `U_global.csv`
  - `rotation_source.txt`
  - optional family-specific summary CSV if already used by the existing comparator flow
- define FEM outputs to compare against the Abaqus exports
- identify comparator entrypoints under [`post_processing/validation_visualisers/static/`](post_processing/validation_visualisers/static/)

Success criteria:

- one canonical linear static job parses and generates a valid Abaqus script
- one family-local comparator is pinned and CI-visible
- contract docs describe the full load-family rollout path

#### 2. Harmonic

Harmonic is currently validated analytically, but not yet via pyAbaqus external reference.

Planned deliverables:

- add [`docs/conventions/HARMONIC_PYABAQUS_REFERENCE_CONTRACT.md`](docs/conventions/HARMONIC_PYABAQUS_REFERENCE_CONTRACT.md)
- define a canonical harmonic benchmark job
- define Abaqus exports for:
  - frequency samples
  - selected DOF displacement response per frequency
  - optional real and imaginary components where applicable
- define FEM-side outputs using the current harmonic runner artifacts
- define comparator and family summary structure under [`post_processing/validation_visualisers/harmonic/`](post_processing/validation_visualisers/harmonic/)

Success criteria:

- harmonic job can be dispatched through the shared Abaqus generator
- extracted Abaqus frequency-response data matches the contract shape
- at least one scalar or curve comparison is defined

#### 3. Eigen

Eigen needs a pyAbaqus contract for frequencies and selected mode-shape comparison.

Planned deliverables:

- add [`docs/conventions/EIGEN_PYABAQUS_REFERENCE_CONTRACT.md`](docs/conventions/EIGEN_PYABAQUS_REFERENCE_CONTRACT.md)
- use [`jobs/job_benchmark_eigen_cantilever/`](jobs/job_benchmark_eigen_cantilever/) as the initial pinned job
- define Abaqus exports for:
  - natural frequencies
  - selected normalized mode-shape coordinates or nodal DOFs
- define FEM artifacts already produced by the eigen path
- define a family comparator under [`post_processing/validation_visualisers/eigen/`](post_processing/validation_visualisers/eigen/)

Success criteria:

- first eigen benchmark produces a pinned modal comparison table
- frequencies can be compared within documented tolerance

#### 4. Linear buckling

Linear buckling should be formalized next to eigen because both are spectral families.

Planned deliverables:

- add [`docs/conventions/LINEAR_BUCKLING_PYABAQUS_REFERENCE_CONTRACT.md`](docs/conventions/LINEAR_BUCKLING_PYABAQUS_REFERENCE_CONTRACT.md)
- use [`jobs/job_benchmark_linear_buckling_column/`](jobs/job_benchmark_linear_buckling_column/) as the initial pinned job
- define Abaqus exports for:
  - load factors
  - selected buckling mode-shape data
- define FEM artifacts already produced by the buckling path
- define family comparator under [`post_processing/validation_visualisers/buckling/`](post_processing/validation_visualisers/buckling/)

Success criteria:

- pinned load-factor comparison exists
- first mode-shape comparison contract is defined

### Wave 2 — deepen already-started families

#### 5. Nonlinear static

Nonlinear static already has the best contract base and should be extended from partial implementation to full suite completeness.

Planned deliverables:

- finish rollout of all six canonical load families from [`docs/conventions/NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md`](docs/conventions/NONLINEAR_STATIC_ABAQUS_VALIDATION_SUITE.md:9)
- standardize mesh ladders for all families
- complete comparison tooling for:
  - tip displacement vs load factor
  - optional interior displacement history
  - optional final deformed-shape parity
- consolidate family summaries under [`post_processing/validation_visualisers/static/`](post_processing/validation_visualisers/static/)

Success criteria:

- all six nonlinear-static families exist with pinned jobs and family-local outputs
- multiple calibrated Abaqus reference pairs are documented and executable

#### 6. Transient

Transient has contract scaffolding and should next gain comparison implementation.

Planned deliverables:

- extend Abaqus extraction for selected DOF time histories
- define DOF mapping metadata format
- add FEM-vs-Abaqus comparator under [`post_processing/validation_visualisers/transient/`](post_processing/validation_visualisers/transient/)
- define canonical transient summary tables and plots

Success criteria:

- transient benchmark job produces both FEM and Abaqus time-history files
- comparison runner can evaluate selected structural response histories

### Wave 3 — advanced path-dependent external-reference validation

#### 7. Nonlinear buckling

Nonlinear buckling should be planned now but implemented after simpler external-reference families are stable.

Why later:

- continuation response is path-dependent
- imperfection amplitude and source matter
- path comparison is more complex than scalar or frequency-table comparison

Planned deliverables:

- add [`docs/conventions/NONLINEAR_BUCKLING_PYABAQUS_REFERENCE_CONTRACT.md`](docs/conventions/NONLINEAR_BUCKLING_PYABAQUS_REFERENCE_CONTRACT.md)
- define the canonical imperfect-column reference case from [`jobs/job_benchmark_nl_buckling_imperfect_column/`](jobs/job_benchmark_nl_buckling_imperfect_column/)
- define Abaqus-side continuation history contract
- define FEM-side continuation history comparison contract
- specify path-level tolerance policy for:
  - load factor vs tip displacement
  - critical-region sampled agreement
  - convergence and branch interpretation limits

Success criteria:

- one imperfect-column external-reference path comparison is reproducible
- continuation summary and comparison plots are family-local and documented

## Shared code implementation work

The following shared files are expected to carry most of the cross-family pyAbaqus implementation burden:

- [`post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py`](post_processing/validation_visualisers/abaqus/job_to_abaqus_script.py)
- [`post_processing/validation_visualisers/abaqus/extract_odb_results.py`](post_processing/validation_visualisers/abaqus/extract_odb_results.py)
- [`post_processing/validation_visualisers/abaqus/simulation_type_dispatch.py`](post_processing/validation_visualisers/abaqus/simulation_type_dispatch.py)
- [`post_processing/validation_visualisers/abaqus/config.py`](post_processing/validation_visualisers/abaqus/config.py)

Shared implementation responsibilities:

- dispatch by canonical simulation type
- embed simulation-family-specific CAE step setup
- export family-specific result artifacts from ODB
- write raw results to the standard Abaqus result directory
- keep extraction logic small and contract-driven

## Family-local code implementation work

Each family layer should own its own comparator entrypoints and summary generation.

Planned family-local homes:

- [`post_processing/validation_visualisers/static/`](post_processing/validation_visualisers/static/)
- [`post_processing/validation_visualisers/transient/`](post_processing/validation_visualisers/transient/)
- [`post_processing/validation_visualisers/harmonic/`](post_processing/validation_visualisers/harmonic/)
- [`post_processing/validation_visualisers/eigen/`](post_processing/validation_visualisers/eigen/)
- [`post_processing/validation_visualisers/buckling/`](post_processing/validation_visualisers/buckling/)

Each family should eventually provide:

- one main suite runner
- one or more comparison helpers
- summary output tables
- plots
- a pass/fail evaluation function

## Testing plan

### Contract-level tests

Add or extend tests so CI can lock the suite shape before all external references are complete.

Test responsibilities:

- canonical benchmark jobs exist
- directory standards exist
- the shared Abaqus generator accepts canonical jobs
- extracted artifact names match the documented contract
- family-local comparator entrypoints exist

Likely primary test files:

- [`tests/test_validation_visualisers.py`](tests/test_validation_visualisers.py)
- [`tests/test_validation_directory_standard.py`](tests/test_validation_directory_standard.py)

### Incremental validation tests

As family implementations land, add tests for:

- file generation contract
- result-shape sanity
- benchmark scalar parity where pinned
- summary pass/fail logic

## Suggested implementation order

1. linear static contract normalization
2. harmonic pyAbaqus contract
3. eigen pyAbaqus contract
4. linear buckling pyAbaqus contract
5. shared Abaqus generator and extractor extensions for those families
6. nonlinear static suite completion
7. transient comparator implementation
8. nonlinear buckling external-reference design and implementation

## Implementation checkpoints

### Checkpoint A — contract completeness

All target families have a conventions document describing:

- pinned job
- raw Abaqus artifact set
- FEM artifact set
- comparison quantity
- summary ownership

### Checkpoint B — shared generator coverage

The shared Abaqus generation and extraction layer can dispatch:

- linear static
- nonlinear static
- transient
- harmonic
- eigen
- linear buckling

Nonlinear buckling may remain planned or partial at this checkpoint.

### Checkpoint C — family comparator coverage

Each primary family has at least one family-local comparator entrypoint and one pinned validation case.

### Checkpoint D — advanced external-reference closure

Nonlinear buckling obtains a documented path-comparison policy and at least one external-reference benchmark.

## Out-of-scope concerns for the first implementation wave

The following should not block initial suite rollout:

- full nonlinear buckling branch switching
- advanced bifurcation tracking
- rich mode-shape correlation metrics beyond an initial selected-DOF contract
- large-scale reference corpus population for every family before the suite shape is stable

## Recommended immediate next action

The next implementation wave should begin by adding the missing family contract documents for:

- linear static
- harmonic
- eigen
- linear buckling

Then extend the shared pyAbaqus generator and ODB extraction flow to satisfy those contracts before building the family-local comparator runners.
