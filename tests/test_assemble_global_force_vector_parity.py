"""Parity and edge cases for ``assemble_global_force_vector`` (shared static/transient scatter)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from processing.dynamic.assembly import assemble_global_force_vector


def test_none_force_vectors_returns_zero():
    F = assemble_global_force_vector([], None, total_dof=3, job_results_dir=None)
    assert F.shape == (3,) and F.dtype == np.float64


def test_empty_force_vector_list_returns_zero():
    F = assemble_global_force_vector([], [], total_dof=6, job_results_dir=None)
    assert F.shape == (6,) and F.dtype == np.float64 and np.all(F == 0.0)


def test_complex_accumulation_matches_reference():
    dof_maps = [
        np.array([0, 1], dtype=np.int32),
        np.array([2, 3], dtype=np.int32),
    ]
    Fe = [np.array([1.0 + 1.0j, 0.0]), np.array([0.0, 2.0])]
    F = assemble_global_force_vector(
        [],
        Fe,
        total_dof=4,
        local_global_dof_map=dof_maps,
    )
    ref = np.zeros(4, dtype=np.complex128)
    for f, d in zip(Fe, dof_maps):
        ref[d] += np.asarray(f, dtype=np.complex128).ravel()
    np.testing.assert_array_equal(F, ref)


def test_shape_mismatch_message():
    dof_maps = [np.array([0, 1], dtype=np.int32)]
    with pytest.raises(ValueError, match="Force vector 0 shape mismatch"):
        assemble_global_force_vector([], [np.array([1.0, 2.0, 3.0])], 4, local_global_dof_map=dof_maps)
