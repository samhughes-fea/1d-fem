# Deprecated and removed element type strings

Authoritative list of **public** `element_dictionary["types"]` / `element.txt` class names by lifecycle bucket. See [`JOB_INPUT_BEAM_WARPING.md`](JOB_INPUT_BEAM_WARPING.md) for the preferred baseline + `[warping]` workflow.

## Bucket A — Removed thin linear `LinearWarping*` aliases

These were thin subclasses of the baseline beam element; behaviour is unchanged using the baseline type with `[warping] = 1` (or `element_dictionary["warping"]`).

| Removed type string | Use instead |
|---------------------|-------------|
| `LinearWarpingEulerBernoulliBeamElement3D` | `LinearEulerBernoulliBeamElement3D` + `[warping]` |
| `LinearWarpingTimoshenkoBeamElement3D` | `LinearTimoshenkoBeamElement3D` + `[warping]` |
| `LinearWarpingLevinsonBeamElement3D` | `LinearLevinsonBeamElement3D` + `[warping]` |
| `LinearWarpingReddyBeamElement3D` | `LinearReddyBeamElement3D` + `[warping]` |

`ElementFactory` raises `ValueError` with a migration message if a removed string appears.

## Bucket B — Removed: nonlinear warping Euler–Bernoulli alias

TL EB with Vlasov warping is implemented in a single class: `NonlinearEulerBernoulliBeamElement3D` (14 local DOFs when `[warping]` enables the mesh).

| Removed type string | Use instead |
|---------------------|-------------|
| `NonlinearWarpingEulerBernoulliBeamElement3D` | `NonlinearEulerBernoulliBeamElement3D` + `[warping]` |

## Bucket C — Deregistered: nonlinear warping Timoshenko (stub)

| Type string | Status |
|-------------|--------|
| `NonlinearWarpingTimoshenkoBeamElement3D` | **Not registered.** (Legacy stub class file removed; no warping NL Timoshenko implementation.) Use `NonlinearTimoshenkoBeamElement3D` for 12-DOF TL Timoshenko until a warping-capable NL Timoshenko ships. |

The package directory `pre_processing/element_library/nonlinear/timoshenko_with_warp/` is reserved (no implementation module); the factory does not load a warping NL Timoshenko class.

## Legacy inference in `beam_warping.py`

If `element_dictionary["warping"]` is absent, `mesh_uses_warping_dof` / `element_warping_stiffness_on` still treat any element *type* string containing `"Warping"` as enabling warping (legacy jobs). Prefer explicit `[warping]` so behaviour does not depend on substring checks. Removed public aliases above are rejected by the factory with a migration error.
