# pre_processing/element_library/linear/reddy/utilities/B_matrix.py
"""
Reddy beam strain-displacement matrix: same strain definitions as Levinson
(γ_xy = du_y/dx - θ_z + α d²θ_z/dx², κ_z = dθ_z/dx, κ_y = dθ_y/dx).
Re-export Levinson StrainDisplacementOperator.
"""

from pre_processing.element_library.linear.levinson.utilities.B_matrix import (
    StrainDisplacementOperator,
)

__all__ = ["StrainDisplacementOperator"]
