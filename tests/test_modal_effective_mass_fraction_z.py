"""Directional effective mass fraction helper (native eigen metrics)."""

from __future__ import annotations

import numpy as np
from scipy.sparse import eye

from processing.eigen.metrics.directional_effective_mass import modal_effective_mass_fraction_z


def test_modal_effective_mass_fraction_z_identity_modes():
    n = 6
    M = eye(n, format="csr", dtype=np.float64)
    Phi = np.eye(n, dtype=np.float64)
    out = modal_effective_mass_fraction_z(M, Phi, dof_per_node=6)
    assert out is not None
    assert out.shape == (n,)
    uz = np.zeros(n)
    uz[2] = 1.0
    r = uz / np.linalg.norm(uz)
    rMr = float(r @ M @ r)
    for j in range(n):
        v = Phi[:, j]
        mjj = float(v @ M @ v)
        lj = float(v @ M @ r)
        expect = (lj * lj) / (rMr * mjj)
        assert abs(out[j] - expect) < 1e-10


def test_modal_effective_mass_fraction_z_skips_non_six_dof():
    M = eye(7, format="csr", dtype=np.float64)
    Phi = np.eye(7, dtype=np.float64)
    assert modal_effective_mass_fraction_z(M, Phi, dof_per_node=7) is None
