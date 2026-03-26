# pre_processing/element_library/nonlinear/timoshenko/utilities/geometric_stiffness.py
"""
Re-export ``GeometricStiffnessOperator`` for nonlinear Timoshenko.

``K_sigma`` is built as a Gauss sum ``+= ... * w_g * detJ``; pass ``dN_dx`` from **Timoshenko** shape functions (linear Lagrange on ``u`` and ``theta``).
Definitions: ``nonlinear/euler_bernoulli/utilities/geometric_stiffness.py`` module docstring.
"""
from pre_processing.element_library.nonlinear.euler_bernoulli.utilities.geometric_stiffness import (
    GeometricStiffnessOperator,
)

__all__ = ["GeometricStiffnessOperator"]
