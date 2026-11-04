"""
Hooks for harmonic tertiary extensions (complex displacement field).

Today, formulation-cache post uses real snapshots from ``real(U)`` / ``imag(U)``
(see ``RESULTS_DESIGN.md``). Native complex-valued stress recovery would extend
orchestrators to consume :math:`\\mathbf{u}_r + i \\mathbf{u}_i` coherently; this
module holds small, dependency-free helpers for that migration.
"""

from __future__ import annotations

import numpy as np


def static_recovery_pair_from_complex_column(U_col: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Return ``(real(u), imag(u))`` as separate float vectors for two-pass static recovery.

    Parameters
    ----------
    U_col
        One column of the harmonic displacement matrix (complex, length n_dof).
    """
    z = np.asarray(U_col, dtype=np.complex128).ravel()
    return np.real(z), np.imag(z)
