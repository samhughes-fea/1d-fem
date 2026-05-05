# Harmonic reference validation: SDOF analytical benchmark

This note defines the first reference-grade harmonic validation target for [`simulation_runner/harmonic/harmonic_simulation.py`](../../simulation_runner/harmonic/harmonic_simulation.py) and the underlying frequency-response kernels in [`processing/harmonic/frequency_response.py`](../../processing/harmonic/frequency_response.py).

## Reference case

Use an undamped SDOF oscillator with:

- mass `m = 1`
- stiffness `k = 100`
- damping `c = 0`
- harmonic forcing amplitude `F = 1`
- forcing frequency `ω = 10 rad/s`

The analytical steady-state response is:

\[
u(\omega) = \frac{F}{k - \omega^2 m}
\]

For the pinned benchmark case above:

- `u = 1 / (100 - 100)`, which is singular at resonance, so use a nearby reference case instead.

Pinned repository validation case:

- `ω = 8 rad/s`
- `u = 1 / (100 - 64) = 1 / 36`

## Acceptance quantity

Compare the computed harmonic displacement against the analytical SDOF displacement magnitude/value for the pinned scalar case.

## Tolerance policy

For the current repository benchmark:

- absolute tolerance `1e-10`
- relative tolerance `1e-9`
