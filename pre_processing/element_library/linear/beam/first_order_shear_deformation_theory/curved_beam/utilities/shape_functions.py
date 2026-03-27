# Curved element uses registry LinearTimoshenkoBeamElement3D; local module for four-file layout.
"""
Re-export Timoshenko ``ShapeFunctionOperator``: same ``N`` (12, 6) per GP as straight Timoshenko.

Used by ``linear_curved_timoshenko_3D`` on the chord map; curvature ``kappa0`` enters ``CurvedStrainDisplacementOperator``, not ``N``.
Weak form: ``F_dist += w_g * N.T @ q * detJ``, ``detJ = L/2``.
"""
from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.utilities.shape_functions import ShapeFunctionOperator
__all__ = ["ShapeFunctionOperator"]
