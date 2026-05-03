"""Unit tests for ``processing.boundary_supports`` penalty BC resolution."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from processing.boundary_supports import resolve_penalty_fixed_dofs


def test_resolve_legacy_first_six_without_prescribed():
    fd = resolve_penalty_fixed_dofs(
        total_dof=24,
        dof_per_node=6,
        prescribed_displacement_dict=None,
        section_settings={},
        grid_node_ids=np.arange(4),
    )
    assert fd is not None
    np.testing.assert_array_equal(fd, np.arange(6, dtype=np.int32))


def test_resolve_fixed_node_id():
    fd = resolve_penalty_fixed_dofs(
        total_dof=18,
        dof_per_node=6,
        prescribed_displacement_dict=None,
        section_settings={"fixed_node_id": 1},
        grid_node_ids=np.arange(3),
    )
    np.testing.assert_array_equal(fd, np.arange(6, 12, dtype=np.int32))


def test_resolve_fixed_node_id_unknown_raises():
    with pytest.raises(ValueError, match="not present"):
        resolve_penalty_fixed_dofs(
            total_dof=12,
            dof_per_node=6,
            prescribed_displacement_dict=None,
            section_settings={"fixed_node_id": 99},
            grid_node_ids=np.arange(2),
        )


def test_prescribed_only_no_extra_fixed():
    pd = {
        "global_dof": np.array([0, 1], dtype=np.int32),
        "value": np.array([0.0, 0.0], dtype=np.float64),
    }
    fd = resolve_penalty_fixed_dofs(
        total_dof=12,
        dof_per_node=6,
        prescribed_displacement_dict=pd,
        section_settings={},
        grid_node_ids=None,
    )
    assert fd is None


def test_prescribed_plus_fixed_node_merges_in_apply_bc():
    pd = {
        "global_dof": np.array([0], dtype=np.int32),
        "value": np.array([0.0], dtype=np.float64),
    }
    fd = resolve_penalty_fixed_dofs(
        total_dof=18,
        dof_per_node=6,
        prescribed_displacement_dict=pd,
        section_settings={"fixed_node_id": 1},
        grid_node_ids=np.arange(3),
    )
    np.testing.assert_array_equal(fd, np.arange(6, 12, dtype=np.int32))
