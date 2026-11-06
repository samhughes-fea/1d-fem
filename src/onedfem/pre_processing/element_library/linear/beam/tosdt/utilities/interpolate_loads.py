# pre_processing/element_library/linear/beam/tosdt/utilities/interpolate_loads.py
"""Distributed-load interpolation for third-order beams — delegates to Timoshenko."""

from onedfem.pre_processing.element_library.linear.beam.fosdt.timoshenko.utilities.interpolate_loads import (
    LoadInterpolationOperator,
)

__all__ = ["LoadInterpolationOperator"]
