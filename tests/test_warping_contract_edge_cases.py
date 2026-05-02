"""
Contract checks for ``beam_warping`` helpers: explicit ``[warping]`` vs legacy type inference,
and effective Γ when stiffness is toggled.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def test_mesh_warping_dof_false_when_warping_column_all_off():
    from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

    ed = {
        "ids": np.array([0]),
        "types": np.array(["LinearEulerBernoulliBeamElement3D"]),
        "warping": np.array([0], dtype=np.int8),
    }
    assert not mesh_uses_warping_dof(ed)


def test_mesh_warping_dof_true_when_explicit_column_on():
    from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

    ed = {
        "ids": np.array([0]),
        "types": np.array(["LinearEulerBernoulliBeamElement3D"]),
        "warping": np.array([1], dtype=np.int8),
    }
    assert mesh_uses_warping_dof(ed)


def test_mesh_warping_dof_legacy_type_containing_warping():
    from pre_processing.element_library.beam_warping import mesh_uses_warping_dof

    ed = {
        "ids": np.array([0]),
        "types": np.array(["LegacyWarpingBeamElement3D"]),
    }
    assert mesh_uses_warping_dof(ed)


def test_element_warping_stiffness_on_prefers_column_over_name():
    from pre_processing.element_library.beam_warping import element_warping_stiffness_on

    ed = {
        "ids": np.array([0]),
        "types": np.array(["LegacyWarpingBeamElement3D"]),
        "warping": np.array([0], dtype=np.int8),
    }
    assert not element_warping_stiffness_on(ed, 0, "LegacyWarpingBeamElement3D")


def test_effective_warping_gamma_zeros_when_stiffness_off():
    from pre_processing.element_library.beam_warping import effective_warping_gamma

    assert effective_warping_gamma(1.0e-6, False) == 0.0


def test_effective_warping_gamma_passes_gamma_when_stiffness_on():
    from pre_processing.element_library.beam_warping import effective_warping_gamma

    g = 2.5e-8
    assert effective_warping_gamma(g, True) == pytest.approx(g)
