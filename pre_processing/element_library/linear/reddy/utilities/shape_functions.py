# pre_processing/element_library/linear/reddy/utilities/shape_functions.py
"""
Reddy shapes: re-export Levinson ``ShapeFunctionOperator`` (same ``N``, ``dN_dxi``, ``d2N_dxi2``, shape ``(n_gp, 12, 6)``).

Distributed loads: ``F_dist += w_g * N.T @ q * detJ``. Strain operator differs via ``alpha`` in Reddy ``B_matrix``.
See ``linear_reddy_3D.py``.
"""

from pre_processing.element_library.linear.levinson.utilities.shape_functions import (
    ShapeFunctionOperator,
)

__all__ = ["ShapeFunctionOperator"]
