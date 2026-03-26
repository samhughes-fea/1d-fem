# pre_processing/element_library/linear/reddy/utilities/B_matrix.py
"""
Reddy beam ``B``: re-export of Levinson ``StrainDisplacementOperator``.

``B`` (6, 12) per Gauss point; ``eps`` (6,) Voigt order matches Levinson (see that module). Reddy differs in ``alpha`` terms inside ``B``.
Parent ``linear_reddy_3D`` uses ``K_e += B.T @ D @ B * w_g * detJ`` with selective integration like Levinson.
"""

from pre_processing.element_library.linear.levinson.utilities.B_matrix import (
    StrainDisplacementOperator,
)

__all__ = ["StrainDisplacementOperator"]
