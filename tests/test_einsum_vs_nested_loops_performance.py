"""
Benchmark: einsum vs nested-loop stiffness assembly.

Moved from pre_processing/element_library/euler_bernoulli/tests/ into fem_model/tests.
Runs a fixed number of trials and asserts both methods run; optional matplotlib plot is skipped in test.
"""

import numpy as np
import time


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


def test_einsum_vs_loops_benchmark_runs():
    """Run a short benchmark and assert einsum and loops both complete (no plot)."""
    num_gauss_points = 3
    dof_size = 6
    dim = 2
    num_trials = 100
    np.random.seed(43)

    einsum_times = []
    loop_times = []
    for _ in range(num_trials):
        weights = np.random.rand(num_gauss_points)
        B_T_tensor = np.random.rand(num_gauss_points, dof_size, dim)
        D = np.random.rand(dim, dim)
        B_tensor = np.random.rand(num_gauss_points, dim, dof_size)
        detJ = np.random.rand()

        t0 = time.perf_counter()
        _ = _compute_stiffness_einsum(weights, B_T_tensor, D, B_tensor, detJ)
        einsum_times.append(time.perf_counter() - t0)

        t0 = time.perf_counter()
        _ = _compute_stiffness_loops(weights, B_T_tensor, D, B_tensor, detJ)
        loop_times.append(time.perf_counter() - t0)

    avg_einsum = np.mean(einsum_times)
    avg_loop = np.mean(loop_times)
    assert avg_einsum >= 0 and avg_loop >= 0
    # Sanity: loops are typically slower; at least both ran
    assert len(einsum_times) == num_trials and len(loop_times) == num_trials
