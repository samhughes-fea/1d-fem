# Higher-order shear (Phase 2b)

This document describes the **higher-order shear** beam formulations implemented in the 1D element library: **Reddy** and **Levinson**. Both use third-order shear kinematics (no shear correction factor κ; shear stiffness GA).

## Implemented elements

- **Reddy:** `LinearReddyBeamElement3D` — 2-node, 12 DOF; same kinematics as Levinson (quintic transverse displacement, cubic rotation); γ = ∂u/∂x − θ + α ∂²θ/∂x²; D-matrix uses GA (no κ).
- **Levinson:** Same strain and material formulation as Reddy; see `pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/levinson/` and `.../reddy/` (Reddy re-exports Levinson B and D).

## Tests

- **Shape functions:** `tests/test_higher_order_shear_shape_functions.py` — interpolation and continuity for Reddy shape functions.
- **B/D dimensions and kinematics:** `tests/test_higher_order_shear_B_D_dimensions.py` — B (6×12), D (6×6), shear rows use higher-order kinematics (displacement gradient, θ, and optionally α d²θ/dx²).
- **Thick beam:** `tests/test_higher_order_shear_thick_beam.py` — cantilever with L/h small; compares Reddy and Timoshenko tip deflection (both valid; may differ in thick regime).

See also existing Reddy/Levinson unit tests under `tests/test_reddy_*.py` and `tests/test_levinson_*.py`.
