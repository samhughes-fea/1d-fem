# Transient reference validation: SDOF Newmark benchmark

This note defines the first reference-grade transient validation target for [`TransientSimulationRunner`](../../simulation_runner/transient/dynamic_simulation.py) and the underlying [`newmark_integrate()`](../../processing/dynamic/time_integration.py).

## Reference case

Use an undamped SDOF oscillator with:

- mass `m = 1`
- stiffness `k = 100`
- damping `c = 0`
- initial displacement `u(0) = 1`
- initial velocity `v(0) = 0`
- zero external forcing

The analytical solution is:

\[
u(t) = \cos(\omega_n t), \quad \omega_n = \sqrt{k/m}
\]

## Acceptance quantity

Compare the numerical displacement history against the analytical displacement history on the same time grid using maximum absolute error.

## Tolerance policy

For the current repository benchmark:

- `max_abs_error <= 2e-2`

This is the first reference-grade transient stability/accuracy gate for the Newmark implementation used in the project.
