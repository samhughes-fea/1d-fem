# Nonlinear (Total Lagrangian) 3D Levinson third-order beam (16 DOF).

from onedfem.pre_processing.element_library.nonlinear.tosdt.nonlinear_third_order_beam import (
    NonlinearThirdOrderShearBeamElement3D,
)


class NonlinearLevinsonBeamElement3D(NonlinearThirdOrderShearBeamElement3D):
    """TL Levinson: Petrov–Galerkin test on six resultants plus GL increments."""

    element_type_name = "Levinson-3D-Nonlinear"
    variational = False
