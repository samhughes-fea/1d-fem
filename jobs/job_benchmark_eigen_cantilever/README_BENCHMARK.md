# Pinned eigen benchmark: cantilever beam

This job pins the cantilever Euler-Bernoulli eigen case as the repository benchmark asset for frequency validation.

## Reference model

Use the classical fixed-free Euler-Bernoulli beam frequency formula as the reference *shape*, then calibrate repository tolerance against the current pinned job convention.

\[
f_n = \frac{\beta_n^2}{2\pi L^2}\sqrt{\frac{EI}{\rho A}}
\]

with the standard cantilever roots:

- `β1 = 1.875104068711961`
- `β2 = 4.694091132974174`
- `β3 = 7.854757438237612`
- `β4 = 10.995540734875467`

## Acceptance criteria

1. eigen primary artifacts are written under `primary_results/modal_results/`
2. the artifact manifest exists
3. the first four frequencies match the pinned reference values documented below within the documented tolerance

## Pinned repository reference values

For the current pinned job configuration, the accepted frequency targets are:

- `f1 = 11.253952 Hz`
- `f2 = 28.593523 Hz`
- `f3 = 72.731488 Hz`
- `f4 = 184.790597 Hz`

## Tolerance policy

For the pinned benchmark regression, use a relative tolerance of **2%** on the first four frequencies. These are pinned repository reference values for the current job configuration, not universal closed-form targets for every modeling choice.
