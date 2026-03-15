# GEBT shear beam formulation (Phase 3a)

Shear-deformable Geometrically Exact Beam Theory element: **GEBTShearBeamElement3D**.

## Interface

- **tangent_stiffness_matrix(U_e):** Tangent stiffness K_T = K_0 + K_σ(U_e), shape (12, 12).
- **internal_force_vector(U_e):** Internal force F_int = ∫ B^T S dx, shape (12,).
- **element_stiffness_matrix():** Returns ElementObject with K_e = tangent at U_e=0 (same as linear Timoshenko K_e for same integration rule).
- **strain_at_gauss_points(U_e):** Strain E_lin + E_nl at each Gauss point (for converged cache).

## Limit at U=0

At zero displacement, K_T equals the linear Timoshenko element stiffness K_e when both use the same section, material, length, and **selective integration** (1-point shear, bending order from element integration_orders). This is enforced in `tests/test_gebt_shear_initial_stiffness_vs_linear.py` and `scripts/verify_gebt_shear_initial_stiffness.py`.

## Relation to Total Lagrangian Timoshenko

The current implementation uses the same kinematics as the Total Lagrangian nonlinear Timoshenko beam (Green–Lagrange strain, geometric stiffness K_σ). So at U≠0 the response matches NonlinearTimoshenkoBeamElement3D. A full current-configuration GEBT formulation (exact director update, strain in deformed config) can be added later; this element satisfies the Phase 3a requirement that at U=0 the tangent matches linear Timoshenko and that the nonlinear runner can use it for Newton–Raphson.

## References

- Plan: Phase 3 GEBT and systematic tests
- Total Lagrangian: [docs/element_library/total_lagrangian_beam_formulation.md](total_lagrangian_beam_formulation.md)
- Nonlinear runner: simulation_runner/static/nonlinear_static_simulation.py
