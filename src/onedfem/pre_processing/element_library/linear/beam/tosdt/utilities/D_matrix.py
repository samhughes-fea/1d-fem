# pre_processing/element_library/linear/beam/tosdt/utilities/D_matrix.py
"""Material stiffness **D** (8, 8) — third-order shear (Reddy/Levinson).

**Constitutive:** ``S = D @ e8`` with
``e8 = [eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x, chi_y, chi_z]``.
``D`` is ``SectionAnalysis.third_order_stiffness``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MaterialStiffnessOperator:
    D: np.ndarray

    def assembly_form(self) -> np.ndarray:
        return np.asarray(self.D, dtype=np.float64)

    def timoshenko_block(self) -> np.ndarray:
        return self.assembly_form()[:6, :6]
