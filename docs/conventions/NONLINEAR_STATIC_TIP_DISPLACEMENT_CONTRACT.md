# Nonlinear static validation: tip displacement vs load contract

Benchmark job:

- [`jobs/job_benchmark_nl_static_cantilever_tip/`](../../jobs/job_benchmark_nl_static_cantilever_tip/)

## Comparison quantity

The first nonlinear-static Abaqus comparison contract is the tip vertical displacement versus load level for the cantilever benchmark.

## Phase 2 scope

This phase defines the contract only. It does not yet implement the full FEM-vs-Abaqus comparison.

Required comparison quantity:

- tip displacement at the loaded node DOF (`UY`)
- evaluated at pinned load levels / nonlinear increments

## Acceptance policy

For the current phase, success means:

1. the benchmark job and validation README exist
2. the comparison contract doc exists
3. regression coverage asserts the contract assets are present

## Next phase

The next phase should add:

- Abaqus export for selected nonlinear load steps
- FEM-vs-Abaqus tip-displacement table/plot generation
- tolerance-backed comparison at pinned load levels
