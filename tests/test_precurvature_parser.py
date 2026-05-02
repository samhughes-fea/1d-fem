"""Tests for precurvature.txt parsing and Voigt reference strain helpers."""

import os
import tempfile

import numpy as np
import pytest

from pre_processing.parsing.precurvature_parser import (
    element_reference_strain_voigt,
    parse_precurvature,
    reference_strain_voigt,
    voigt_standard_beam_to_third_order_beam,
)


def test_parse_precurvature_missing_file_returns_zeros():
    ids = np.array([1, 2], dtype=np.int64)
    out = parse_precurvature("", ids)
    assert out.shape == (2, 3)
    assert np.allclose(out, 0.0)

    out2 = parse_precurvature("/nonexistent/path/precurvature.txt", ids)
    assert np.allclose(out2, 0.0)


def test_parse_precurvature_reads_section():
    ids = np.array([10, 20], dtype=np.int64)
    content = """# comment
[Precurvature]
[element_id] [k_x0] [k_y0] [k_z0]
10 0.01 0.02 0.03
20 0.0 0.0 0.0
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False, encoding="utf-8") as f:
        f.write(content)
        path = f.name
    try:
        out = parse_precurvature(path, ids)
        assert out[0, 0] == pytest.approx(0.01)
        assert out[0, 1] == pytest.approx(0.02)
        assert out[0, 2] == pytest.approx(0.03)
        assert np.allclose(out[1], 0.0)
    finally:
        os.unlink(path)


def test_reference_strain_voigt_ordering():
    e0 = reference_strain_voigt(np.array([0.1, 0.2, 0.3]))
    assert e0[0] == 0.0
    assert e0[1] == pytest.approx(0.2)
    assert e0[2] == pytest.approx(0.3)
    assert e0[3] == 0.0
    assert e0[4] == 0.0
    assert e0[5] == pytest.approx(0.1)


def test_voigt_standard_beam_to_third_order_beam_swaps_curvatures():
    e = reference_strain_voigt(np.array([0.1, 0.2, 0.3]))
    t = voigt_standard_beam_to_third_order_beam(e)
    assert t[1] == pytest.approx(e[2])
    assert t[2] == pytest.approx(e[1])
    assert t[0] == pytest.approx(e[0])
    assert t[5] == pytest.approx(e[5])


def test_element_reference_strain_voigt_from_dictionary():
    ed = {
        "ids": np.array([5, 7], dtype=np.int64),
        "precurvature_per_element": np.array([[0.0, 0.1, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64),
    }
    v5 = element_reference_strain_voigt(ed, 5)
    assert v5[1] == pytest.approx(0.1)
    v7 = element_reference_strain_voigt(ed, 7)
    assert np.allclose(v7, 0.0)

    ed2 = {"ids": np.array([0], dtype=np.int64)}
    v0 = element_reference_strain_voigt(ed2, 0)
    assert np.allclose(v0, 0.0)
