# UL stress resultants: ``S = D @ eps`` on the current configuration.

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class StressResultantOperator:
    def from_strain(self, D: np.ndarray, eps: np.ndarray) -> np.ndarray:
        return np.asarray(D, dtype=np.float64) @ np.asarray(eps, dtype=np.float64).ravel()
