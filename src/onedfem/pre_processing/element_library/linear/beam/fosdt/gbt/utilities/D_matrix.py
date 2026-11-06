# pre_processing/element_library/linear/beam/fosdt/gbt/utilities/D_matrix.py
"""Material stiffness **D** — GBT section matrices (infinitesimal strain).

**Constitutive:** ``S = D @ eps`` with the GBT strain packing
``[eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x, b'(p), b(p), a''(m), a'(m), a(m)]``.
``D`` is element-signed ``SectionAnalysis.gbt_matrices``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """Wraps the GBT stiffness (and optional mass) for assembly."""

    D: np.ndarray
    M: np.ndarray | None = None

    def assembly_form(self) -> np.ndarray:
        return np.asarray(self.D, dtype=np.float64)

    def mass_form(self) -> np.ndarray | None:
        return None if self.M is None else np.asarray(self.M, dtype=np.float64)
