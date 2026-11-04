# Job input: follower loads (`follower_load.txt`)

A follower load is a nodal force or moment that turns with the node it acts on, such as
a thrust that stays along a deforming wing or a pressure resultant on a turning section.
It is given in the node's undeformed orientation, in global components:

```
[Follower load]
# x  y  z   F_x F_y F_z   M_x M_y M_z
1.0  0  0   0   250 0     0   0   0
```

- Each row has nine numbers. The point must be a node (to 1e-6 of the model size).
- The load scales with the load factor of the nonlinear static runner, like `point_load.txt`,
  and the two files can be used together (dead loads in `point_load.txt`).
- Only nonlinear static jobs (`static_nonlinear`) take it. Other job types stop with an error,
  because a follower load in a linear job is the same as a dead load.

## What the runner does

With `theta` the node's total rotation vector and `R = exp(theta)`, the generalised forces are

- `R F_0` on the node's translations;
- `J(theta)^T R M_0` on its rotation vector, with `J = I + (1 - cos t)/t^2 K + (t - sin t)/t^3 K^2`
  (`K = [theta x]`, `t = |theta|`), because a spatial moment does work on the spin `J delta_theta`.

The load depends on the displacement, so every Newton iteration also adds the load
stiffness `K_L = -d g / dU` to the tangent. It is unsymmetric, because a follower load is
not conservative. The force part is closed form (`[R F_0 x] J`), and the moment part is a
central difference of the closed-form moment. Code: `processing/static/follower_loads.py`.

The rotation vector is the corotational elements' nodal rotation, so use
`CorotationalEulerBernoulliBeamElement3D`, `CorotationalTimoshenkoBeamElement3D` or
`CorotationalBenscoterBeamElement3D` for large rotations.

## Validation

`tests/test_follower_loads.py` and [examples/follower_loads](../../examples/follower_loads/README.md):
a cantilever under a tip force that stays normal to the tip, against the elastica with a
follower force, and the same beam under a dead load against the dead-load elastica.
