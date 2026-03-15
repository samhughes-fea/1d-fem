# pre_processing/element_library/linear/reddy/utilities/D_matrix.py
"""
Reddy beam material matrix: same as Levinson (6×6, GA for shear, no κ).
Re-export Levinson MaterialStiffnessOperator.
"""

from pre_processing.element_library.linear.levinson.utilities.D_matrix import (
    MaterialStiffnessOperator,
)

__all__ = ["MaterialStiffnessOperator"]
