from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class HarmonicComplexSnapshot:
    component_name: str
    U_global: np.ndarray
    results_subdir: str


def build_harmonic_complex_snapshots(
    *,
    U_column: np.ndarray,
    frequency_index: int,
    displacement_component: str,
) -> list[HarmonicComplexSnapshot]:
    comp = str(displacement_component).strip().lower()
    if comp not in {"real", "imag", "both", "complex_components"}:
        raise ValueError(
            "harmonic displacement component must be 'real', 'imag', 'both', or 'complex_components'"
        )
    Ucol = np.asarray(U_column, dtype=np.complex128).ravel()
    out: list[HarmonicComplexSnapshot] = []
    if comp in {"real", "both", "complex_components"}:
        out.append(
            HarmonicComplexSnapshot(
                component_name="real",
                U_global=np.real(Ucol),
                results_subdir=f"harmonic_post/freq_{frequency_index:04d}",
            )
        )
    if comp in {"imag", "both", "complex_components"}:
        out.append(
            HarmonicComplexSnapshot(
                component_name="imag",
                U_global=np.imag(Ucol),
                results_subdir=f"harmonic_post/freq_{frequency_index:04d}_imag",
            )
        )
    return out
