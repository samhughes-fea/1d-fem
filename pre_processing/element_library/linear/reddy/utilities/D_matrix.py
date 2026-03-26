# pre_processing/element_library/linear/reddy/utilities/D_matrix.py
"""
Reddy ``D``: re-export of Levinson ``MaterialStiffnessOperator`` — ``D`` (6, 6), ``S = D @ eps``, ``G*A`` shear (no ``kappa``).
"""

from pre_processing.element_library.linear.levinson.utilities.D_matrix import (
    MaterialStiffnessOperator,
)

__all__ = ["MaterialStiffnessOperator"]
