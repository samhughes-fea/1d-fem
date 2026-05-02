# Large-rotation beam kinematics vs Total Lagrangian (TL)

This project keeps **resultant beam strains** in **Voigt form on the line** (see `docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md`). Full rank-2 Green–Lagrange storage at every Gauss point as in solid elements is **not** the default beam architecture.

## Total Lagrangian (`NonlinearTimoshenkoBeamElement3D`, `NonlinearEulerBernoulliBeamElement3D`)

- Geometry and strain are referred to the **fixed mesh** (reference chord).
- Large rotations appear through **Green–Lagrange** strain and stress resultants built from that fixed reference.

## Co-rotational (`CorotationalBeamElement3D`)

- A **moving** orthonormal basis follows the **current chord** between nodes.
- A **small-strain** linear Timoshenko or Euler–Bernoulli stiffness is formed in that basis and rotated to global DOFs.
- Distinct from TL: no TL Green–Lagrange assembly in global axes; objective split between rigid rotation of the chord and small strain in the corotated frame.

## Geometrically exact shear-deformable registration (`GeometricallyExactShearDeformableBeam3D`)

- **Currently** uses the **same TL Timoshenko weak form and tangent** as `NonlinearTimoshenkoBeamElement3D` for a working nonlinear path.
- **Future:** director-based classical GEBT weak form with consistent tangent; until then, compare against TL for regression in the small-strain limit.

## Practical comparison

| Aspect | TL nonlinear beams | Co-rotational |
|--------|---------------------|---------------|
| Reference | Fixed undeformed geometry | Current chord each iteration |
| Strain measure | Green–Lagrange on reference line | Infinitesimal in corotated frame |
| Typical use | Same TL pipeline as documented | Large rigid rotation of element with small incremental strain in local frame |

Verification tests live under `tests/test_large_rotation_beam_kinematics.py`.

## Co-rotational tangent: `finite_difference` vs `elastic_material`

[`CorotationalBeamElement3D`](../../pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_3D.py) exposes `tangent_stiffness_mode` (constructor / optional job setting `simulation_settings["nonlinear"]["corotational_tangent_mode"]`).

| Mode | Mechanism | When to use |
|------|-----------|-------------|
| **`finite_difference`** (default) | Central differences on `internal_force_vector` — **consistent** with the force path; includes implicit dependence of the corotated frame on `U_e` via those probes. | **Default for production.** Use when chord rotations are large, Newton stalls or diverges, or when validating against reference solutions / other codes. Full spin stiffness is captured indirectly through the force sampling. |
| **`elastic_material`** | Analytic **`Tᵀ K_local(L) K`** only (symmetrized); **no** explicit corotational **spin** stiffness. Cheaper (no 12 extra force calls per tangent). | Use when speed matters **and** incremental chord motion / strain in the corotated frame stay moderate; validate convergence on representative jobs before relying on it globally. |

If Newton iterations explode or residual stagnates after switching to `elastic_material`, revert to **`finite_difference`** for that case.

Job-level smoke (minimal `static_nonlinear`, both tangent modes, `newton_converged` in `primary_summary.csv`): [`tests/test_corotational_nonlinear_job_tangent_modes.py`](../../tests/test_corotational_nonlinear_job_tangent_modes.py).

See also [`docs/conventions/API_STANDARDS.md`](../conventions/API_STANDARDS.md) (simulation runners).
