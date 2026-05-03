"""Regression: ill-scaled pencil still solves; dense_threshold is honored."""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix, eye

from processing.spectral.operations.solve_generalized_eigenproblem import SolveGeneralizedEigenproblem


def test_solve_generalized_eigenproblem_ill_scaled_pencil():
    """Large stiffness scale vs identity mass — should return finite frequencies."""
    n = 24
    scale = 1e18
    K = csr_matrix(scale * np.eye(n))
    M = csr_matrix(np.eye(n))
    ev, vec, fh = SolveGeneralizedEigenproblem(
        num_modes=3,
        dense_threshold=n + 10,
        job_results_dir=None,
    ).run(K, M)
    assert ev.shape == (3,)
    assert vec.shape == (n, 3)
    assert fh.shape == (3,)
    assert np.all(np.isfinite(fh))
    assert np.all(fh >= 0.0)
