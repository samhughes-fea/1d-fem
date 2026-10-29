# pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/levinson/linear_levinson_3D.py
"""
2-node 3D Levinson beam: third-order shear kinematics with vectorial equilibrium.

The kinematics, DOFs and operators are those of ``ThirdOrderShearBeamElement3D``
(``../third_order_beam.py``), the same as the Reddy element. Levinson (1981) takes
the stress resultants of the higher-order kinematics but writes only the classical
equilibrium equations (``N' = 0``, ``M' = V``, ``V' = -q``, ``T' = 0``); the
higher-order resultant has no equation of its own.

**Weak form (Petrov-Galerkin):** ``K_e = int B6^T D8[:6, :] B8 dx``: classical virtual
strains against the full resultants. ``K_e`` is not symmetric when the section
couples the slaved warping to the classical strains. For an uncoupled section (any
isotropic one) it reduces to Timoshenko with the section's shear factor (``5/6``
for a rectangle), the known result for Levinson's beam; the higher-order kinematics
still set the stresses. Levinson's boundary conditions are the classical ones, so
the slopes ``Q1``, ``Q2`` are usually left free at a clamp.

See ``docs/element_library/higher_order_shear.md``.
"""

from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.third_order_beam import (
    ThirdOrderShearBeamElement3D,
)


class LinearLevinsonBeamElement3D(ThirdOrderShearBeamElement3D):
    """Levinson third-order shear beam (isotropic or thin-walled sections)."""

    element_type_name = "Levinson-3D"
    variational = False
