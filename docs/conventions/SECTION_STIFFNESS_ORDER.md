# Section stiffness order and BeamDyn

1d-fem works in its own beam axes throughout, and so does its section provider
(`onedfem.sections`, anisotropic-section-model by default):

- beam axis **x**, section axes **(y, z)**;
- generalised strains `[eps, kappa_y, kappa_z, gamma_xy, gamma_xz, phi]`;
- section forces `[N, M_y, M_z, V_y, V_z, T]`.

The provider hands 1d-fem right-handed matrices (section contract 2.0,
`SectionAnalysis.frame == "engine"`). 1d-fem's beam `B` matrices have `kappa_y = -theta_y'`,
so an element re-signs its section once, when it takes it
(`onedfem.sections.to_element_section`, `ELEMENT_SIGN = diag(1, -1, 1, 1, 1, 1)`):
`K_element = ELEMENT_SIGN @ K @ ELEMENT_SIGN`, and the copy has `frame == "element"`. The
flip is its own inverse (`to_engine_section`). Tapers and GBT alignment work on the
provider's sections before any element sees them; pretwist goes through
`element_pretwist`, which re-signs back and forth around the provider's `with_pretwist`.

Two other frames appear where 1d-fem meets other tools. anisotropic-section-model's `axes`
module holds every conversion; the summary:

**anisotropic-section-model's public frame** (its top-level `ThinWalledSection`,
`SolidSection`, `TimoshenkoBeam`): beam axis **z**, section axes **(x, y)** with x to the
trailing edge and y to the suction side, and BeamDyn's order `[Fx, Fy, N, Mx, My, Mz]`
against `[gamma_x, gamma_y, eps, kappa_x, kappa_y, kappa_z]`. It is a cyclic relabelling,
1d-fem (x, y, z) = public (z, x, y): both are right-handed, no signs change, and a 6 x 6
stiffness converts by reordering:

```python
import numpy as np

p = [3, 4, 0, 1, 2, 5]          # public position -> 1d-fem index (axes.FROM_ENGINE)
K_public = K_onedfem[np.ix_(p, p)]

q = np.argsort(p)               # [2, 3, 4, 0, 1, 5], the way back (axes.TO_ENGINE)
K_onedfem = K_public[np.ix_(q, q)]
```

A 6 x 6 mass matrix is on the displacements and rotations, whose relabelling is
`[2, 0, 1, 5, 3, 4]` (1d-fem's `[u_x, u_y, u_z, theta_x, theta_y, theta_z]` from the public
`[u_z, u_x, u_y, theta_z, theta_x, theta_y]`).

**BeamDyn's blade frame** (a BeamDyn blade file, windIO): X towards the suction side and Y
towards the trailing edge, so it is the mirror image of the public frame, not a rotation of
it. Going there swaps x and y and flips the sign of the moments and rotations
(`axes.MIRROR`, `axes.to_beamdyn_frame`). Twist needs no conversion: BeamDyn's
`initial_twist` is the public twist.

BeamDyn takes the stiffness about its reference axis: move the 1d-fem stiffness there first
(`anisotropic_section_model.axes.stiffness_about(K, reference)`) and then convert.

Sections given as tables in any of the three frames can go straight into a job:
`section_matrix.txt` names the frame of each section and the provider converts it
(`SectionAnalysis.from_matrices`; [format](JOB_INPUT_THIN_WALLED_SECTIONS.md#section_matrixtxt-optional)).
