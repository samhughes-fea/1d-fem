# Transient multi-DOF reference validation

This document defines the first structural multi-DOF transient benchmark beyond SDOF analytical checks.

## Benchmark case

Pinned validation job: [`jobs/job_benchmark_transient_cantilever_multidof/`](../../jobs/job_benchmark_transient_cantilever_multidof/)

The case is a cantilever beam-line model using the transient runner with a short-duration harmonic-type forcing history.

## Comparison quantities

Use the following repository-calibrated quantities:

- tip displacement history at the final node translational DOF
- maximum absolute tip displacement over the time window
- one interior-node displacement history sample or peak

## Phase 1 acceptance policy

For the current phase, this is a pinned structural benchmark with repository-calibrated quantities, not yet an external-reference benchmark.

## Tolerance policy

Use benchmark-specific tolerances recorded next to the pinned quantities when the first calibrated run is frozen.
