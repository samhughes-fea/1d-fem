# pre_processing/element_library/nonlinear/timoshenko/utilities/geometric_stiffness.py
"""
\\(\\mathbf{K}_\\sigma\\) for nonlinear Timoshenko — same weak-form implementation as EB (**re-export**).

**Tangent context:** Parent uses \\(\\mathbf{K}_T = \\mathbf{K}_0 + \\mathbf{K}_\\delta + \\mathbf{K}_\\sigma\\); ``GeometricStiffnessOperator`` supplies **only** \\(\\mathbf{K}_\\sigma\\).
Pass ``dN_dx`` from **Timoshenko** shape functions (linear Lagrange on displacements and rotations). Definitions and Gauss terms:
``nonlinear/euler_bernoulli/utilities/geometric_stiffness.py``.
"""
from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)

__all__ = ["GeometricStiffnessOperator"]
