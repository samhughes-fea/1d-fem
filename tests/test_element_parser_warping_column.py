"""Tests for optional ``[warping]`` and combined ``[curvature]`` columns in ``element.txt``."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from pre_processing.parsing.element_parser import ElementParser


def _write(tmp: Path, body: str) -> str:
    p = tmp / "element.txt"
    p.write_text(body, encoding="utf-8")
    return str(p)


def test_parse_base_no_optional_columns(tmp_path):
    txt = """[Element]
[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order]
0 0 1 LinearEulerBernoulliBeamElement3D 3 3 3 2 2 3 2
"""
    fp = _write(tmp_path, txt)
    out = ElementParser(fp, str(tmp_path)).parse()
    ed = out["element_dictionary"]
    assert ed["warping"].shape == (1,)
    assert ed["warping"][0] == 0
    assert ed["curvature"][0] == 0.0


def test_parse_warping_only_column(tmp_path):
    txt = """[Element]
[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]
0 0 1 LinearEulerBernoulliBeamElement3D 3 3 3 2 2 3 2 1
"""
    fp = _write(tmp_path, txt)
    out = ElementParser(fp, str(tmp_path)).parse()
    ed = out["element_dictionary"]
    assert ed["warping"][0] == 1
    assert ed["curvature"][0] == 0.0


def test_parse_curvature_and_warping(tmp_path):
    txt = """[Element]
[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [curvature] [warping]
0 0 1 LinearEulerBernoulliBeamElement3D 3 3 3 2 2 3 2 0.01 0
"""
    fp = _write(tmp_path, txt)
    out = ElementParser(fp, str(tmp_path)).parse()
    ed = out["element_dictionary"]
    assert np.isclose(ed["curvature"][0], 0.01)
    assert ed["warping"][0] == 0


def test_parse_warping_true_false_strings(tmp_path):
    txt = """[Element]
[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]
0 0 1 LinearEulerBernoulliBeamElement3D 3 3 3 2 2 3 2 false
1 1 2 LinearEulerBernoulliBeamElement3D 3 3 3 2 2 3 2 True
"""
    fp = _write(tmp_path, txt)
    out = ElementParser(fp, str(tmp_path)).parse()
    ed = out["element_dictionary"]
    assert ed["warping"][0] == 0 and ed["warping"][1] == 1


def test_parse_warping_levinson_and_reddy(tmp_path):
    txt = """[Element]
[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order] [warping]
0 0 1 LinearLevinsonBeamElement3D 3 3 3 3 3 3 2 1
1 1 2 LinearReddyBeamElement3D 3 3 3 3 3 3 2 0
"""
    fp = _write(tmp_path, txt)
    out = ElementParser(fp, str(tmp_path)).parse()
    ed = out["element_dictionary"]
    assert ed["warping"][0] == 1 and ed["warping"][1] == 0
