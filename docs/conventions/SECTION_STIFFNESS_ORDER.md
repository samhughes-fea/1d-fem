# Section stiffness order and BeamDyn

1d-fem works in its own beam axes throughout, and so does its section provider
(`onedfem.sections`, anisotropic-section-model by default):

- beam axis **x**, section axes **(y, z)**;
- generalised strains `[eps, kappa_y, kappa_z, gamma_xy, gamma_xz, phi]`;
- section forces `[N, M_y, M_z, V_y, V_z, T]`.

BeamDyn (and the windIO and anisotropic-section-model's top-level `ThinWalledSection`) uses
beam axis **z**, section axes **(x, y)**, and the order
`[Fx, Fy, N, Mx, My, Mz]` against `[gamma_x, gamma_y, eps, kappa_x, kappa_y, kappa_z]`.

The two axis systems are a cyclic relabelling, 1d-fem (x, y, z) = BeamDyn (z, x, y), so
both are right handed and no signs change. A 6 x 6 stiffness converts by reordering:

```python
import numpy as np

p = [3, 4, 0, 1, 2, 5]          # BeamDyn position -> 1d-fem index
K_beamdyn = K_onedfem[np.ix_(p, p)]

q = np.argsort(p)               # [2, 3, 4, 0, 1, 5], the way back
K_onedfem = K_beamdyn[np.ix_(q, q)]
```

anisotropic-section-model exports `p` as `thin_wall.BEAMDYN_FROM_ONEDFEM`. A 6 x 6 mass
matrix converts the same way.

BeamDyn takes the stiffness about its reference axis; move the 1d-fem stiffness there first
(`anisotropic_section_model.thin_walled.stiffness_about(K, y0, z0)`) and then reorder.
