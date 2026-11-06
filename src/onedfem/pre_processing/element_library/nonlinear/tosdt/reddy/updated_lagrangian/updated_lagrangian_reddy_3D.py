# Updated-Lagrangian 3D Reddy third-order beam (16 DOF).

from onedfem.pre_processing.element_library.nonlinear.tosdt.updated_lagrangian_third_order import (
    UpdatedLagrangianThirdOrderShearBeamElement3D,
)


class UpdatedLagrangianReddyBeamElement3D(UpdatedLagrangianThirdOrderShearBeamElement3D):
    """UL Reddy: variational ``K_0`` on the current chord plus beam-column ``K_sigma``."""

    element_type_name = "Reddy-3D-UpdatedLagrangian"
    variational = True
