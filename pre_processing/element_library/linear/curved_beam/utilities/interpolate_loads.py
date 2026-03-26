# Distributed loads: same operator as straight Timoshenko (x-based table).
"""
Re-export Timoshenko ``LoadInterpolationOperator`` for curved element: ``F_dist += w_g * N.T @ q * detJ``.
"""
from pre_processing.element_library.linear.timoshenko.utilities.interpolate_loads import LoadInterpolationOperator
__all__ = ["LoadInterpolationOperator"]
