# pre_processing/element_library/nonlinear/__init__.py
# Nonlinear (geometric) element formulations. Classes added as implemented.

from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.total_lagrangian.nonlinear_euler_bernoulli_3D import (
    NonlinearEulerBernoulliBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.total_lagrangian.nonlinear_timoshenko_3D import (
    NonlinearTimoshenkoBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.tosdt.levinson.nonlinear_levinson_3D import (
    NonlinearLevinsonBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.tosdt.reddy.nonlinear_reddy_3D import (
    NonlinearReddyBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.euler_bernoulli.updated_lagrangian.updated_lagrangian_euler_bernoulli_3D import (
    UpdatedLagrangianEulerBernoulliBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.timoshenko.updated_lagrangian.updated_lagrangian_timoshenko_3D import (
    UpdatedLagrangianTimoshenkoBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.tosdt.levinson.updated_lagrangian.updated_lagrangian_levinson_3D import (
    UpdatedLagrangianLevinsonBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.tosdt.reddy.updated_lagrangian.updated_lagrangian_reddy_3D import (
    UpdatedLagrangianReddyBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_3D import (
    CorotationalBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section import (
    CorotationalBenscoterBeamElement3D,
    CorotationalEulerBernoulliBeamElement3D,
    CorotationalGBTBeamElement3D,
    CorotationalLevinsonBeamElement3D,
    CorotationalReddyBeamElement3D,
    CorotationalTimoshenkoBeamElement3D,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.geometrically_exact_shear_deformable_beam.geometrically_exact_shear_deformable_beam_3D import (
    GeometricallyExactShearDeformableBeam3D,
)
from onedfem.pre_processing.element_library.nonlinear.large_rotations.gebt_unshearable.gebt_unshearable_3D import (
    GEBTUnshearableBeamElement3D,
)

__all__ = [
    "NonlinearEulerBernoulliBeamElement3D",
    "NonlinearTimoshenkoBeamElement3D",
    "NonlinearLevinsonBeamElement3D",
    "NonlinearReddyBeamElement3D",
    "UpdatedLagrangianEulerBernoulliBeamElement3D",
    "UpdatedLagrangianTimoshenkoBeamElement3D",
    "UpdatedLagrangianLevinsonBeamElement3D",
    "UpdatedLagrangianReddyBeamElement3D",
    "CorotationalBeamElement3D",
    "CorotationalEulerBernoulliBeamElement3D",
    "CorotationalTimoshenkoBeamElement3D",
    "CorotationalBenscoterBeamElement3D",
    "CorotationalGBTBeamElement3D",
    "CorotationalReddyBeamElement3D",
    "CorotationalLevinsonBeamElement3D",
    "GeometricallyExactShearDeformableBeam3D",
    "GEBTUnshearableBeamElement3D",
]
