"""``onedfem.processing.eigen.smallest_generalized_eigenpairs`` — parity with dense pencil."""

from __future__ import annotations

import numpy as np
from scipy.sparse import eye

from onedfem.processing.eigen.smallest_generalized_eigenpairs import smallest_generalized_eigenpairs


def test_smallest_generalized_sparse_branch_matches_dense() -> None:
    n = 24
    K = eye(n, format="csr") * 55.0
    M = eye(n, format="csr") * 0.4
    lam_d, Phi_d = smallest_generalized_eigenpairs(K, M, 5, dense_threshold=10_000)
    lam_s, Phi_s = smallest_generalized_eigenpairs(K, M, 5, dense_threshold=0)
    np.testing.assert_allclose(lam_d, lam_s, rtol=1e-6, atol=1e-8)


def test_massless_dofs_do_not_spoil_the_lowest_modes() -> None:
    """DOFs without mass (GBT warping amplitudes) give infinite eigenvalues, not low ones."""
    import scipy.sparse as sp
    from scipy.linalg import eigh

    rng = np.random.default_rng(0)
    n, n0 = 700, 60
    A = rng.standard_normal((n, n))
    K = A @ A.T + n * np.eye(n)
    m = rng.uniform(0.5, 2.0, n)
    m[:n0] = 0.0                                    # massless DOFs
    lam, _ = smallest_generalized_eigenpairs(sp.csr_matrix(K), sp.diags(m).tocsr(), 6, dense_threshold=0)
    # reference: condense the massless DOFs out of K (static condensation is exact for them)
    i0, i1 = np.arange(n0), np.arange(n0, n)
    Kc = K[np.ix_(i1, i1)] - K[np.ix_(i1, i0)] @ np.linalg.solve(K[np.ix_(i0, i0)], K[np.ix_(i0, i1)])
    ref = eigh(Kc, np.diag(m[i1]), eigvals_only=True)[:6]
    np.testing.assert_allclose(lam, ref, rtol=1e-8)
