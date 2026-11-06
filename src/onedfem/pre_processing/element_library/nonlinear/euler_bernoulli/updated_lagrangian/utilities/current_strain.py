# Small (Almansi/Cauchy-reduced) strain on the **current** UL chord: discrete ``E = B U_e``.

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CurrentConfigurationStrainOperator:
    """``eps = B @ U_e`` with linear ``B`` evaluated on the current length."""

    def strain(self, B: np.ndarray, U_e: np.ndarray) -> np.ndarray:
        return (np.asarray(B, dtype=np.float64) @ np.asarray(U_e, dtype=np.float64).ravel()).ravel()
