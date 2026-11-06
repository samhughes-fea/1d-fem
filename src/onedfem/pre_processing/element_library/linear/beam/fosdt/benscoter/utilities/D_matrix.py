# pre_processing/element_library/linear/beam/fosdt/benscoter/utilities/D_matrix.py
"""Material stiffness **D** (8, 8) — Benscoter section law (infinitesimal strain).

**Constitutive:** ``S = D @ eps`` with
``eps = [eps_x, kappa_y, kappa_z, gamma_xy, gamma_xz, phi_x, d, psi']``.
``D`` is the element-signed ``SectionAnalysis.benscoter_stiffness``. Leading 6×6 is
Saint-Venant; ``D[6,6]`` warping shear; ``D[7,7] = E Γ``.

**Stiffness:** ``K_e += B.T @ D @ B * w_g * detJ`` (selective reduced rows in the orchestrator).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MaterialStiffnessOperator:
    """Wraps the 8×8 Benscoter section stiffness for assembly."""

    D: np.ndarray

    def assembly_form(self) -> np.ndarray:
        return np.asarray(self.D, dtype=np.float64)

    def postprocess_form(self) -> np.ndarray:
        return self.assembly_form()
