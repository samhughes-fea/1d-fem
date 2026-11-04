# pre_processing/element_library/linear/beam/third_order_shear_deformation_theory/reddy/linear_reddy_3D.py
"""
2-node 3D Reddy beam: third-order shear, variationally consistent.

The kinematics, DOFs and operators are those of ``ThirdOrderShearBeamElement3D``
(``../third_order_beam.py``): Timoshenko section motion plus the shear warping
slaved to the shear strain, ``eps_xx += Psi(s) . chi`` with
``chi = gamma' - R psi6``. Per node ``[u, v, w, theta_x, theta_y, theta_z, v', w']``.

**Weak form:** ``K_e = int B8^T D8 B8 dx`` (4-point Gauss), so the higher-order
resultant ``P = D8[6:, :] e8`` has its own equilibrium equation (Reddy 1984). At a
clamp fix ``Q1`` and ``Q2`` (``v' = w' = 0``) as well as the six standard DOFs to
restrain the shear warping (``w = w' = phi = 0``).

See ``docs/element_library/higher_order_shear.md``.
"""

from pre_processing.element_library.linear.beam.third_order_shear_deformation_theory.third_order_beam import (
    ThirdOrderShearBeamElement3D,
)


class LinearReddyBeamElement3D(ThirdOrderShearBeamElement3D):
    """Reddy third-order shear beam (isotropic or thin-walled sections)."""

    element_type_name = "Reddy-3D"
    variational = True
