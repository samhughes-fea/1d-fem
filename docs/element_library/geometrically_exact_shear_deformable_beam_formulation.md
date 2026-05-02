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

Stub class: `GeometricallyExactShearDeformableBeam3D` in `pre_processing/element_library/nonlinear/geometrically_exact_shear_deformable_beam/`. Assembly methods raise `NotImplementedError` until the formulation above is coded.
