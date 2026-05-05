# Transient reference validation: Rayleigh-damped cantilever case

This note defines the next transient reference-grade validation step beyond the SDOF analytical benchmark.

## Reference case

Use the existing Rayleigh-damped cantilever transient regression built from the shared cantilever modal case helper used in [`tests/test_transient_forcing_rayleigh.py`](../../tests/test_transient_forcing_rayleigh.py).

Pinned configuration:

- `time_step = 0.01`
- `end_time = 0.03`
- `rayleigh_alpha = 0.02`
- `rayleigh_beta = 1e-6`

## Pinned repository reference quantity

For the current benchmark configuration, use the maximum absolute displacement from the transient history as the pinned scalar.

Pinned repository target:

- `max_abs_displacement = 9.517172091466624e-10`

## Acceptance criteria

1. transient displacement history exists
2. transient primary summary exists
3. summary reports `damping_source = rayleigh`
4. `max_abs_displacement` in the summary matches the pinned reference quantity within tolerance

## Tolerance policy

Use absolute tolerance `1e-12` for the pinned summary scalar in the current repository benchmark configuration.
