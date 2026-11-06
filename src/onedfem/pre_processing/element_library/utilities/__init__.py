# pre_processing\element_library\utilities\__init__.py

# Import shape function sets
from onedfem.pre_processing.element_library.utilities.dof_mapping import expand_stiffness_matrix, expand_force_vector
from onedfem.pre_processing.element_library.utilities.gauss_quadrature import get_gauss_points, integrate_vector, integrate_matrix
from onedfem.pre_processing.element_library.utilities.jacobian import compute_jacobian_matrix, compute_jacobian_determinant, general_jacobian_and_determinant
# Legacy TL names: prefer nonlinear.{euler_bernoulli,timoshenko}.total_lagrangian.utilities.
from onedfem.pre_processing.element_library.utilities.total_lagrangian_beam import (
    GreenLagrangeStrainOperator,
    StressResultantOperator,
    GeometricStiffnessOperator,
)

# Define explicitly exported members
__all__ = [
    "expand_stiffness_matrix",
    "expand_force_vector",
    "get_gauss_points",
    "integrate_vector",
    "integrate_matrix",
    "compute_jacobian_matrix",
    "compute_jacobian_determinant",
    "general_jacobian_and_determinant",
    "GreenLagrangeStrainOperator",
    "StressResultantOperator",
    "GeometricStiffnessOperator",
]