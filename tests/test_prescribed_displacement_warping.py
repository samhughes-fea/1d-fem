"""Tests for prescribed_displacement_parser with 7 DOF/node (warping mesh)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.parsing.prescribed_displacement_parser import parse_prescribed_displacement


def test_parse_prescribed_displacement_dof_per_node_7(tmp_path: Path) -> None:
    p = tmp_path / "prescribed_displacement.txt"
    p.write_text(
        "[Prescribed Displacement]\n"
        "[id]     [node_id]  [dof]   [value]     [type]          [comment]\n"
        "0        0          UX      0.0         displacement\n"
        "1        1          UY      0.0         displacement\n"
        "2        0          CHI     0.0         displacement\n"
        "3        1          W       0.01        displacement\n",
        encoding="utf-8",
    )

    d = parse_prescribed_displacement(str(p), dof_per_node=7)
    np.testing.assert_array_equal(d["global_dof"], np.array([0, 8, 6, 13]))
    np.testing.assert_array_equal(d["dof_index"], np.array([0, 1, 6, 6]))
    assert list(d["dof"]) == ["UX", "UY", "CHI", "W"]


def test_parse_prescribed_displacement_invalid_dof_per_node(tmp_path: Path) -> None:
    p = tmp_path / "prescribed_displacement.txt"
    p.write_text(
        "[Prescribed Displacement]\n"
        "0        0          UX      0.0         displacement\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="dof_per_node"):
        parse_prescribed_displacement(str(p), dof_per_node=5)
