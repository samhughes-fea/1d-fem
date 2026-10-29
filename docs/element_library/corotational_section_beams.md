# Corotational beams on anisotropic sections (large displacement)

`CorotationalEulerBernoulliBeamElement3D`, `CorotationalTimoshenkoBeamElement3D` and
`CorotationalBenscoterBeamElement3D`
(`pre_processing/element_library/nonlinear/large_rotations/corotational/corotational_section.py`)
take any linear element with its `SectionAnalysis` as a local core. They give large
displacements and rotations with the coupled section stiffness, shear-centre and
tension-centre offsets, restrained warping and large twist. Run them with
`[Type] static_nonlinear`.

## 1. Kinematics

- **Nodal rotations** are total rotation vectors, `R_i = exp(θ_i) E_0`, where `E_0` is
  the undeformed element frame. The additive Newton update of the nonlinear runner is
  consistent with this.
- **Element frame** `R_r`: its first axis follows the current chord. The second is the mean
  of the nodal second axes, made orthogonal to the chord (Battini and Pacoste 2002).
- **Local deformation** `p`, in the core element's DOF order:
  - the chord stretch `l − l0`;
  - the nodal rotations relative to the frame, `θ̄_i = log(R_rᵀ R_i)`;
  - the extra DOFs (the Benscoter warping amplitude `ψ`), which are scalars of the section
    and pass through unchanged.

  A rigid motion gives `p = 0` exactly.

## 2. Local energy

```
Π(p) = ½ pᵀ K p + ½ pᵀ K_σ(p) p + Π_4(p)
```

- `K` is the core's linear stiffness: coupled section, warping for Benscoter.
- `K_σ(p)` is the second-order geometric stiffness of [second_order.md](second_order.md),
  linear in `p`. The cubic term is the work of the stresses on the quadratic strains of the
  rigid in-plane section motion. It carries the offsets, the Wagner term and the bimoment.
- `Π_4 = ½ ∫ Σ_k k ε₂²` with `ε₂ = ½[(v' − zφ')² + (w' + yφ')²]` completes the axial strain of
  large twist. The sum runs over the section's axial points (`SectionAnalysis.axial_points`:
  the wall Gauss points with `K_xx ds`, or a 3 × 3 rectangle for `section.txt` sections).

The local force `f = ∂Π/∂p` and Hessian `H` are exact. The global force is `F = Bᵀ f`
with `B = ∂p/∂U` from the complex step, so the residual is exact to round-off. The tangent
is `Bᵀ H B + Σ f_i ∂²p_i/∂U²`, with the second derivatives from central differences.

Not covered yet: precurvature (prebend, presweep and pretwist come in the next PR), GBT,
and the quartic terms of the shear strains.

## 3. Checks (`tests/test_large_displacement.py`)

| case | reference | result |
|---|---|---|
| Rigid motion of 2 rad | zero strain and force | exact |
| Tangent | central differences of the force | 1e-6 |
| Tip load `PL²/EI = 1` | Bisshopp-Drucker elastica 0.30172 L | 0.01 % (EB) |
| Tip moment `πEI/L` | semicircle | 0.1 % at 32 elements |
| I-beam, uniform twist, free warping | `T = GJk + ½E(I_rr − I_p²/A)k³` | 1e-6 |
| I-beam, warping restrained, 6 kN m | nonlinear Vlasov equation (collocation) | +0.25 % (Benscoter warping shear) |
| I-beam LTB, uniform moment | Trahair 6.24 with prebuckling deflection | < 1 % (EB and Benscoter) |

The LTB check follows the nonlinear path with a small lateral imperfection and finds the
load where the tangent stops being positive definite. That load is above the classical
`M_cr`, by the prebuckling-deflection factor (3.6 % here), as it should be. The classical
value is the linear buckling of [second_order.md](second_order.md).

The example [examples/thin_walled/large_twist](../../examples/thin_walled/large_twist/README.md)
runs the restrained large-twist case through the nonlinear job runner.
