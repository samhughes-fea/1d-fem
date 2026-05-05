# Eigen and linear buckling benchmark suite

Pinned benchmark jobs:

- [`jobs/job_benchmark_eigen_cantilever/`](../../jobs/job_benchmark_eigen_cantilever/)
- [`jobs/job_benchmark_linear_buckling_column/`](../../jobs/job_benchmark_linear_buckling_column/)

## Purpose

This suite promotes eigen and linear buckling from smoke-only execution toward acceptance-level benchmark coverage.

## Acceptance criteria

### Eigen

- frequencies file exists
- mode-shapes file exists
- primary artifact manifest exists
- frequencies are positive
- first benchmark frequencies match the pinned repository reference values within the documented tolerance

### Linear buckling

- buckling load-factors file exists
- buckling mode-shapes file exists
- primary artifact manifest exists
- load factors are finite and positive
- the first benchmark load factor matches the pinned reference value within the documented tolerance

## Scope note

These are repository acceptance benchmarks, not yet literature-tolerance closures. They lock artifact contract and baseline positivity/shape expectations while full external-reference benchmarks are still being curated.
