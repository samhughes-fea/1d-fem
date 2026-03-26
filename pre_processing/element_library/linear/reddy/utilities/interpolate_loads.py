# pre_processing/element_library/linear/reddy/utilities/interpolate_loads.py
"""Re-export Levinson ``LoadInterpolationOperator`` — ``F_dist += w_g * N.T @ q * detJ`` (``detJ = L/2``)."""

from pre_processing.element_library.linear.levinson.utilities.interpolate_loads import (
    LoadInterpolationOperator,
)

__all__ = ["LoadInterpolationOperator"]
