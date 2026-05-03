"""
Optional micro-benchmark for §4 harmonic frequency sweep linear solves.

Not run in CI. Usage (from repo root)::

    python benchmarks/harmonic_sparse_microbench.py

Uses ``harmonic_linear_solver=splu`` vs ``spsolve`` on a medium sparse identity-like system.
**Profiling gate:** compare serial wall time across solvers; if the sweep dominates job time,
enable ``parallel_frequency_sweep`` in job settings. CSC data-buffer (``FEM_HARMONIC_SPLU_CSC_BUFFER``)
only reduces Python allocation overhead — full numeric refactor still runs each frequency.
Symbolic-only SuperLU reuse for changing **A(ω)** values is **not** in SciPy's high-level ``splu``;
a future optional backend (e.g. Pardiso) would require explicit dependency + parity tests.
"""

from __future__ import annotations

import os
import time

import numpy as np
from scipy.sparse import eye

from processing.harmonic.frequency_response import harmonic_damping_matrix, sweep_displacements


def main() -> None:
    n = 400
    K = eye(n, format="csr", dtype=np.float64) * 50.0
    M = eye(n, format="csr", dtype=np.float64) * 0.5
    zeta = 0.02
    omega_ref = 20.0
    C = harmonic_damping_matrix(M, K, zeta, omega_ref, 0.0, 0.0)
    F = np.zeros(n)
    F[0] = 1.0
    freqs = np.linspace(1.0, 80.0, 64)

    for solver in ("spsolve", "splu"):
        t0 = time.perf_counter()
        sweep_displacements(K, M, C, F, freqs, linear_solver=solver, parallel=False)
        dt = time.perf_counter() - t0
        print(f"{solver:8s}  serial  n={n}  n_freq={len(freqs)}  wall_s={dt:.4f}")

    t0 = time.perf_counter()
    sweep_displacements(
        K, M, C, F, freqs, linear_solver="splu", parallel=True, max_workers=4
    )
    dt = time.perf_counter() - t0
    print(f"{'splu':8s}  parallel(4)  n={n}  n_freq={len(freqs)}  wall_s={dt:.4f}")

    for buf in ("1", "0"):
        os.environ["FEM_HARMONIC_SPLU_CSC_BUFFER"] = buf
        t0 = time.perf_counter()
        sweep_displacements(K, M, C, F, freqs, linear_solver="splu", parallel=False)
        dt = time.perf_counter() - t0
        print(
            f"{'splu':8s}  serial  FEM_HARMONIC_SPLU_CSC_BUFFER={buf}  wall_s={dt:.4f}"
        )


if __name__ == "__main__":
    main()
