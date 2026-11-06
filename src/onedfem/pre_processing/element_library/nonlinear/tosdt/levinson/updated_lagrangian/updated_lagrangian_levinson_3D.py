# Updated-Lagrangian 3D Levinson third-order beam (16 DOF).

from onedfem.pre_processing.element_library.nonlinear.tosdt.updated_lagrangian_third_order import (
    UpdatedLagrangianThirdOrderShearBeamElement3D,
)


class UpdatedLagrangianLevinsonBeamElement3D(UpdatedLagrangianThirdOrderShearBeamElement3D):
    """UL Levinson: Petrov–Galerkin ``K_0`` on the current chord plus beam-column ``K_sigma``."""

    element_type_name = "Levinson-3D-UpdatedLagrangian"
    variational = False
