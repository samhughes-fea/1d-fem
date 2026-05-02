# Classical geometrically exact shear-deformable beam (planned)

This document scopes a **future** native implementation of a shear-deformable geometrically exact beam (GEBT) with **finite rotations** and director-based kinematics (Simo & Vu-Quoc–type), distinct from the chord-frame **Total Lagrangian Timoshenko** element (`NonlinearTimoshenkoBeamElement3D`) which uses fixed linear `B`, Green–Lagrange strain, and `assemble_timoshenko_K0` for the material stiffness.

## Relation to existing elements

| Aspect | TL Timoshenko (`NonlinearTimoshenkoBeamElement3D`) | Classical shear-deformable GEBT (this track) |
|--------|-----------------------------------------------------|-----------------------------------------------|
| Frame | Reference chord; small-strain operator `B` | Moving/director frame or equivalent parametrisation |
| Rotations | Incremental / Green–Lagrange pipeline | Finite rotation update (e.g. quaternion or Rodrigues) |
| Weak form | Documented in `total_lagrangian_beam_formulation.md` | Weak form from beam GE literature |
| Quadrature | `TimoshenkoQuadratureOrders` + `loop_order` for TL loops | To be fixed once kernels are implemented |

## Planned ingredients

- **DOFs:** Standard 6 DOF per node placeholder (12 per element) may evolve with rotation parametrisation choice.
- **Strains and resultants:** Work-conjugate pairs per chosen beam GE theory (stretch, shear strains, curvatures; axial force, shear forces, moments, torque).
- **Weak form:** Virtual work statement leading to `F_int` and tangent `K_T` with consistent linearisation.
- **References:** J. C. Simo & L. Vu-Quoc, *Computer Methods in Applied Mechanics and Engineering* (1986-onwards); G. Jelenić & M. A. Crisfield, *Computers & Structures* (1999); optional K. J. Bathe beam formulations.

## Implementation status

`GeometricallyExactShearDeformableBeam3D` is registered under `pre_processing/element_library/nonlinear/large_rotations/geometrically_exact_shear_deformable_beam/`. The **locked weak form**, Gauss assembly for `F_int` and `K_T`, and Voigt strain hook are documented in [`gesdb_weak_form.md`](gesdb_weak_form.md). The strain map routes through `gesdb_kinematics.gesdb_director_voigt_strain_timoshenko_12`, which for the present 2-node Timoshenko interpolation matches chord-frame Green–Lagrange Voigt strains used by `NonlinearTimoshenkoBeamElement3D`, so regressions against the TL parent remain tight.

Optional **`gesdb_tl_fallback`** bypasses the GESDB strain hook; **`gesdb_kernel`** selects **`tl_locked`** (chord-frame TL equivalence) versus **`native`** (engineering axial stretch + linear Timoshenko rows; see [`gesdb_weak_form.md`](gesdb_weak_form.md)).

**Tests:** `tests/test_gesdb_milestone_kinematics.py`.

See also `docs/element_library/large_rotation_vs_total_lagrangian.md`.

## Follow-on milestones

Future work can replace the strain map with a distinct Simo–Vu-Quoc `B` operator while keeping the same Gauss structure documented in [`gesdb_weak_form.md`](gesdb_weak_form.md).

Regression coverage: `tests/test_large_rotation_beam_kinematics.py`, `tests/test_gesdb_milestone_kinematics.py`.
