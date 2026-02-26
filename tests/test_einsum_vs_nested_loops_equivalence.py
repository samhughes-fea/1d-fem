"""
Test that einsum and nested-loop stiffness assembly give equivalent results.

Moved from pre_processing/element_library/euler_bernoulli/tests/ into fem_model/tests.
"""

import numpy as np
import pytest


def _compute_stiffness_einsum(weights, B_T_tensor, D, B_tensor, detJ):
    return np.einsum("i,ijk,kl,ilm->jm", weights, B_T_tensor, D, B_tensor) * detJ


def _compute_stiffness_loops(weights, B_T_tensor, D, B_tensor, detJ):
    num_gauss_points = weights.shape[0]
    size = B_T_tensor.shape[1]
    Ke = np.zeros((size, size))
    for i in range(num_gauss_points):
        for j in range(size):
            for m in range(size):
                sum_term = 0.0
                for k in range(2):
                    for l in range(2):
                        sum_term += B_T_tensor[i, j, k] * D[k, l] * B_tensor[i, l, m]
                Ke[j, m] += weights[i] * sum_term * detJ
    return Ke


def test_einsum_and_nested_loops_stiffness_agree():
    """Einsum and nested-loop stiffness computation must match within numerical tolerance."""
    num_gauss_points = 3
    dof_size = 6
    dim = 2
    np.random.seed(42)
    weights = np.random.rand(num_gauss_points)
    B_T_tensor = np.random.rand(num_gauss_points, dof_size, dim)
    D = np.random.rand(dim, dim)
    B_tensor = np.random.rand(num_gauss_points, dim, dof_size)
    detJ = np.random.rand()

    Ke_einsum = _compute_stiffness_einsum(weights, B_T_tensor, D, B_tensor, detJ)
    Ke_loops = _compute_stiffness_loops(weights, B_T_tensor, D, B_tensor, detJ)

    np.testing.assert_allclose(Ke_einsum, Ke_loops, atol=1e-10)
