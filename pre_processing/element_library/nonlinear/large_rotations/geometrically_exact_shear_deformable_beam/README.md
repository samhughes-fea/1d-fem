# Geometrically exact shear-deformable beam (stub)

This package reserves the factory type `GeometricallyExactShearDeformableBeam3D` for a **classical** shear-deformable geometrically exact beam (finite rotations, director-based kinematics), not the chord Total Lagrangian Timoshenko implementation.

## Milestones

1. Weak form and strain measures per Simo & Vu-Quoc (and related) literature — see `docs/element_library/geometrically_exact_shear_deformable_beam_formulation.md`.
2. Consistent tangent and internal force at Gauss points.
3. Quadrature policy derived from that formulation.

## Migration

For nonlinear analysis that previously used `GEBTShearBeamElement3D` (TL Timoshenko stack), use `NonlinearTimoshenkoBeamElement3D` with mesh `integration_orders`; material stiffness uses `assemble_timoshenko_K0` with the same selective quadrature as linear Timoshenko.
