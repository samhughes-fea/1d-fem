"""IEA 15 MW section validation (post_processing/validation_visualisers/iea15mw_sections)."""

import importlib.util
import os

import numpy as np
import pytest

pytest.importorskip("implicit_geometry_kernel")

ROOT = os.path.join(os.path.dirname(__file__), "..")
RUN = os.path.join(ROOT, "post_processing", "validation_visualisers", "iea15mw_sections", "run.py")


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def iea():
    return _load(RUN, "iea15mw_sections_run")


def test_root_tube_all_models_agree(iea, tmp_path):
    s = iea.main(stations=(0,), plot=False, out_dir=str(tmp_path))["stations"]["0"]
    for model in ("solid_2d", "shell", "prism"):
        for q in iea.QUANTITIES:
            assert abs(s["ratios"][model][q]) < 1e-3, (model, q)
    # PreComp measures wall areas on the outer surface: 2 % high for t/R = 4 %
    for q in iea.QUANTITIES:
        assert 0.015 < s["ratios"]["precomp"][q] < 0.025
    assert s["prism_vs_2d_6x6"]["max_normalised_difference"] < 1e-3


def test_seventy_percent_span_shell_matches_thin_walled(iea, tmp_path):
    s = iea.main(stations=(20,), plot=False, out_dir=str(tmp_path))["stations"]["20"]
    r = s["ratios"]["shell"]
    assert abs(r["EA"]) < 0.01 and abs(r["EI_flap"]) < 0.01 and abs(r["EI_edge"]) < 0.01
    assert abs(r["GJ"]) < 0.03
    assert s["centre_offsets"]["shell_sc"] < 0.01 and s["centre_offsets"]["shell_tc"] < 0.01
    # the 2D caps reach the webs' outer faces; with their carbon matched EA and EI agree
    c = s["carbon"]
    assert 0.05 < 1.0 / c["ratio"] - 1.0 < 0.08
    assert abs(c["matched_ratios"]["EA"]) < 0.01
    assert abs(c["matched_ratios"]["EI_flap"]) < 0.015


def test_summary_conventions(iea):
    K6 = np.diag([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    out = iea.summary_1dfem(K6, (0.0, 0.0), (0.0, 0.0), 1.0)
    assert (out["EA"], out["EI_flap"], out["EI_edge"], out["GJ"]) == pytest.approx((1.0, 2.0, 3.0, 6.0))


def test_seventy_percent_span_prism_end_restraint():
    """The longer prism keeps EA and EI and loses most of its extra shear and torsion stiffness."""
    import json

    ref = os.path.join(ROOT, "post_processing", "validation_visualisers", "iea15mw_sections", "reference")
    short, long = (np.linalg.inv(np.array(json.load(open(os.path.join(ref, f)))["K"]))
                   for f in ("prism_20.json", "prism_20_4chords.json"))
    ratio = np.diag(short) / np.diag(long)       # long-prism over short-prism effective stiffness
    assert np.abs(ratio[:3] - 1.0).max() < 1e-3   # N, M_y, M_z
    assert 0.80 < ratio[3] < 0.88                 # edgewise shear (V_y, y chordwise)
    assert 0.83 < ratio[5] < 0.89                 # torsion
