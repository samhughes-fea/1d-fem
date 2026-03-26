# pre_processing/element_library/linear/reddy/utilities/shape_functions.py
"""
Reddy shapes: re-export Levinson ``ShapeFunctionOperator`` — **``N``**, **``dN_dξ``**, **``d2N_dξ2``** shapes unchanged (**(n_gp, 12, 6)** etc.).
"""

from pre_processing.element_library.linear.levinson.utilities.shape_functions import (
    ShapeFunctionOperator,
)

__all__ = ["ShapeFunctionOperator"]
