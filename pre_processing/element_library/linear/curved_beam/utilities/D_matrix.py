# pre_processing/element_library/linear/curved_beam/utilities/D_matrix.py
"""
Re-exports ``MaterialStiffnessOperator`` from ``linear/timoshenko/utilities/D_matrix.py``.

``D`` (6, 6) and Voigt layout are straight Timoshenko; curvature ``kappa0`` enters only ``CurvedStrainDisplacementOperator`` (``B``).
Weak form remains ``K_e += B.T @ D @ B * w_g * detJ`` in ``linear_curved_timoshenko_3D.py``.
"""
from pre_processing.element_library.linear.timoshenko.utilities.D_matrix import MaterialStiffnessOperator

__all__ = ["MaterialStiffnessOperator"]
