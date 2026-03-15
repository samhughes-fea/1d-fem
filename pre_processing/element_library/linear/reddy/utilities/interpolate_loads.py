# pre_processing/element_library/linear/reddy/utilities/interpolate_loads.py
"""Reddy beam load interpolation: same as Levinson. Re-export LoadInterpolationOperator."""

from pre_processing.element_library.linear.levinson.utilities.interpolate_loads import (
    LoadInterpolationOperator,
)

__all__ = ["LoadInterpolationOperator"]
