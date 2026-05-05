# Non-static production-readiness checklist

This checklist defines the minimum standard required for non-static analysis families to be treated as production-ready at the same operational level as the static stack rooted in [`processing/static/results/RESULTS_DESIGN.md`](../../processing/static/results/RESULTS_DESIGN.md).

It applies to:

- eigen
- buckling
- transient
- harmonic
- shared spectral / transient processing backends that support those runners

## 1. Canonical taxonomy and package ownership

- Canonical runner path is documented and stable.
- Canonical processing package path is documented and stable.
- Deprecated aliases are explicitly listed with removal policy.
- README, docstrings, and import examples use canonical names.

## 2. Primary artifact contract

Each family must define, document, and regression-protect:

- primary output directory
- filenames and shapes
- units / normalization convention
- machine-readable artifact manifest behavior
- diagnostics log locations
- summary artifact expectations

## 3. Staged processing contract

Each non-static family must expose a stable staged-processing contract with:

- clearly named stage classes
- stable `run(...)` entry points
- deterministic telemetry stage names
- documented BC handling
- explicit solver-path selection semantics

## 4. Post-processing contract

Each family must explicitly declare one of the following:

- full secondary/tertiary parity with static
- snapshot-based formulation-cache post-processing with documented limits
- no secondary/tertiary support yet, with justification

If snapshot-driven post-processing is used, the family must define:

- how the snapshot `U_global` is chosen
- whether results are authoritative or convenience-only
- how multiple snapshots/frequencies/modes are laid out on disk

## 5. Validation ladder

Each family must maintain all of the following where applicable:

- unit tests for kernels/utilities
- smoke/integration tests for runner dispatch and artifact creation
- pinned benchmark jobs
- tolerance-backed analytical or external-reference checks
- regression coverage for deprecated compatibility paths where still supported

## 6. Diagnostics and summaries

Each family must provide:

- telemetry under `diagnostics/`
- per-stage logs under `logs/`
- stable summary artifact(s) for operators
- enough metadata to reconstruct the solver path used

## 7. Declared limitations

Each family must document:

- supported physics scope
- unsupported regimes
- approximations or reduced-fidelity post-processing paths
- benchmark classes still pending

## Family-specific artifact expectations

### Eigen

- primary frequencies and mode shapes
- normalization policy documented
- modal lightweight secondary metrics documented where written
- pinned acceptance benchmark job and artifact checks

### Buckling

- linear buckling load factors and mode shapes
- prestress source documented (`linear_static`, `nonlinear_static`, etc.)
- nonlinear continuation history and summary when nonlinear buckling is used
- pinned linear-buckling acceptance benchmark job and artifact checks

### Transient

- time grid
- displacement / velocity / acceleration histories
- damping path metadata (assembled `C` vs Rayleigh fallback)
- snapshot-index policy for optional secondary/tertiary post

### Harmonic

- frequency grid
- complex displacement outputs (real / imag / abs / phase)
- solve-path metadata (direct vs modal superposition)
- damping metadata and post-processing mode declaration
- explicit contract for harmonic complex snapshot export (`real`, `imag`, `both`, or `complex_components`)

## Shared backend expectations

### Spectral backend

- stable stage interfaces for eigen and buckling
- cross-family regression ownership
- diagnostics independent of runner-specific post-processing

### Transient backend

- stable assembly / BC / integration interfaces
- explicit damping precedence and diagnostics
- artifact contract aligned with runner output docs

## Release gate

A non-static family should not be called “production-ready” until every section above is either:

- implemented and tested, or
- explicitly documented as deferred with a scoped limitation note.
