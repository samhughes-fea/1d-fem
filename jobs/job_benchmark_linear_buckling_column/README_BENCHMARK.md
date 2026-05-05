# Pinned linear buckling benchmark: Euler-type column

This job pins the Euler-type cantilever/column buckling case as the repository benchmark asset for linear buckling validation.

## Reference model

Use the classical Euler cantilever critical load

\[
P_{cr} = \frac{\pi^2 E I}{4 L^2}
\]

with the same fixed-free boundary condition convention used by the benchmark job.

## Pinned repository reference scalar

For the current pinned job configuration, the first buckling load factor target is:

- `λ1 = 23.526405`

This repository-pinned value is used for regression calibration of the current benchmark job configuration.

## Acceptance criteria

1. buckling load-factor and mode-shape primary artifacts are written
2. the artifact manifest exists
3. the first buckling load factor matches the pinned reference value within tolerance

## Tolerance policy

Use a relative tolerance of **5%** on the first buckling load factor for the current benchmark job configuration.
