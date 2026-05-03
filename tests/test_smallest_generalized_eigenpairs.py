"""``processing.eigen.smallest_generalized_eigenpairs`` — parity with dense pencil."""

from __future__ import annotations

import numpy as np
from scipy.sparse import eye

from processing.eigen.smallest_generalized_eigenpairs import smallest_generalized_eigenpairs


def test_smallest_generalized_sparse_branch_matches_dense() -> None:
    n = 24
    K = eye(n, format="csr") * 55.0
    M = eye(n, format="csr") * 0.4
    lam_d, Phi_d = smallest_generalized_eigenpairs(K, M, 5, dense_threshold=10_000)
    lam_s, Phi_s = smallest_generalized_eigenpairs(K, M, 5, dense_threshold=0)
    np.testing.assert_allclose(lam_d, lam_s, rtol=1e-6, atol=1e-8)
