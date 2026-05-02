# Removed: `GEBTShearBeamElement3D`

The type string **`GEBTShearBeamElement3D`** and package `nonlinear/gebt_shear/` have been **removed**. That class was a Total Lagrangian Timoshenko stack with selective material stiffness — not a classical finite-rotation geometrically exact beam.

## Migration

| Old | New |
|-----|-----|
| `GEBTShearBeamElement3D` | `NonlinearTimoshenkoBeamElement3D` |

Keep the same mesh **`integration_orders`** (axial, bending_y, bending_z, shear_y, shear_z, torsion). Material stiffness **`K_0`** / linear **`K_e`** both use **`assemble_timoshenko_K0`** with **`TimoshenkoQuadratureOrders`** resolved from the element row (default **shear block** = 1 Gauss point for the shear stiffness rows, matching prior selective integration).

TL internal force and **`K_sigma`** use a single Gauss rule of order **`loop_order`** (default `max` of those mesh orders, at least 2).

## Classical GEBT (future)

For a **literature** shear-deformable geometrically exact beam (finite rotations), see [geometrically_exact_shear_deformable_beam_formulation.md](geometrically_exact_shear_deformable_beam_formulation.md) and factory type **`GeometricallyExactShearDeformableBeam3D`** (stub).

Total Lagrangian reference: [total_lagrangian_beam_formulation.md](total_lagrangian_beam_formulation.md).
