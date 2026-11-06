# Nonlinear (Total Lagrangian) 3D Reddy third-order beam (16 DOF).

"""Linear tosdt Hermite ``B``/``D`` plus Timoshenko Green–Lagrange on the 12 standard DOFs.

See ``nonlinear_third_order_beam.NonlinearThirdOrderShearBeamElement3D``.
"""

from onedfem.pre_processing.element_library.nonlinear.tosdt.nonlinear_third_order_beam import (
    NonlinearThirdOrderShearBeamElement3D,
)


class NonlinearReddyBeamElement3D(NonlinearThirdOrderShearBeamElement3D):
    """TL Reddy: variational ``K = int B8.T D B8`` plus GL ``K_delta`` and Hessian ``K_sigma``."""

    element_type_name = "Reddy-3D-Nonlinear"
    variational = True
