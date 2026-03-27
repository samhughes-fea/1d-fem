# Distributed loads: same operator as straight Timoshenko (x-based table).
"""
Re-export Timoshenko ``LoadInterpolationOperator`` for curved element: ``F_dist += w_g * N.T @ q * detJ``.
"""
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.interpolate_loads import LoadInterpolationOperator
__all__ = ["LoadInterpolationOperator"]
