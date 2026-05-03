"""Harmonic structural BCs use the same penalty ``fixed_dofs`` resolution as spectral."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, eye

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from processing.boundary_supports import resolve_penalty_fixed_dofs
from processing.harmonic.operations.modify_harmonic_structural import ModifyHarmonicStructuralMatrices
from processing.spectral.operations.modify_spectral_global import ModifySpectralGlobalSystem


def test_harmonic_and_spectral_bc_dofs_match_for_fixed_node():
    n = 12
    K = csr_matrix(eye(n, format="csr", dtype=np.float64))
    M = csr_matrix(eye(n, format="csr", dtype=np.float64))
    F = np.ones(n, dtype=np.float64)
    fixed = resolve_penalty_fixed_dofs(
        total_dof=n,
        dof_per_node=6,
        prescribed_displacement_dict=None,
        section_settings={"fixed_node_id": 1},
        grid_node_ids=np.arange(2, dtype=np.int32),
    )
    assert fixed is not None
    with tempfile.TemporaryDirectory() as td:
        _, _, bc_s = ModifySpectralGlobalSystem(
            job_results_dir=td,
            fixed_dofs=fixed,
            prescribed_displacements=None,
        ).run(K.copy(), M.copy())
        Kh, Mh, bc_h, Fh = ModifyHarmonicStructuralMatrices(
            prescribed_displacements=None,
            job_results_dir=td,
            fixed_dofs=fixed,
        ).run(K.copy(), M.copy(), F.copy())
    np.testing.assert_array_equal(np.sort(bc_s), np.sort(bc_h))
    assert Kh.shape == (n, n) and Mh.shape == (n, n)
    assert Fh.shape == (n,)
