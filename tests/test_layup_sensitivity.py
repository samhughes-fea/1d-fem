"""Sensitivities to ply angles and thicknesses: section derivatives and job adjoints."""

from __future__ import annotations

import csv
import logging
import os
import shutil
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "examples", "thin_walled"))

import common  # noqa: E402

from pre_processing.parsing.layup_sensitivity_parser import (  # noqa: E402
    Response,
    format_layup_sensitivity_file,
    parse_layup_sensitivity_file,
)
from pre_processing.section_integration.section_analysis import SectionAnalysis  # noqa: E402
from pre_processing.section_integration.thin_walled import ThinWalledSection, Wall  # noqa: E402
from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial, WallLaminate  # noqa: E402
from pre_processing.section_integration.thin_walled.layup_sensitivity import (  # noqa: E402
    LayupParameter,
    section_derivative,
)
from pre_processing.section_integration.thin_walled.section_analysis import ThinWalledSectionAnalysis  # noqa: E402

CF = PlyMaterial(135e9, 10e9, 5e9, 0.3, 5e9, 3.8e9, 1550.0)
GF = PlyMaterial(40e9, 10e9, 4e9, 0.28, 4e9, 3.5e9, 1900.0)
B, H = 0.6, 0.3
BOX = ThinWalledSection({1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)},
                        [Wall(1, 1, 2, 1), Wall(2, 2, 3, 2), Wall(3, 3, 4, 1), Wall(4, 4, 1, 2)])
LAYUP = {1: [(CF, 45, 5e-4), (CF, -45, 5e-4), (CF, 0, 5e-4), (CF, 20, 5e-4), (GF, 0, 8e-4), (CF, -45, 5e-4),
             (CF, 45, 5e-4)],
         2: [(GF, 45, 5e-4), (GF, -30, 5e-4), (GF, -45, 7e-4), (GF, 45, 5e-4)]}
_WORK = tempfile.mkdtemp()


def teardown_module(module):
    shutil.rmtree(_WORK, ignore_errors=True)


def _box(change=None, hoop="free"):
    lams = {}
    for lid, plies in LAYUP.items():
        pl = [list(p) for p in plies]
        if change is not None and change[0].laminate_id == lid:
            pl[change[0].ply][1 if change[0].kind == "angle" else 2] += change[1]
        lams[lid] = WallLaminate([Ply(*p) for p in pl])
    return ThinWalledSectionAnalysis(BOX, lams, hoop=hoop)


PARAMS = [(LayupParameter(1, 3, "angle"), 1e-3), (LayupParameter(2, 1, "angle"), 1e-3),
          (LayupParameter(1, 4, "thickness"), 1e-8), (LayupParameter(2, 2, "thickness"), 1e-8)]


@pytest.mark.parametrize("hoop", ["free", "restrained"])
@pytest.mark.parametrize("prm, h", PARAMS, ids=[p.label for p, _ in PARAMS])
def test_section_stiffness_derivative_matches_central_differences(prm, h, hoop):
    d = section_derivative(_box(hoop=hoop), prm)
    fd = (_box((prm, h), hoop).stiffness - _box((prm, -h), hoop).stiffness) / (2 * h)
    assert np.abs(d.stiffness - fd).max() <= 1e-6 * np.abs(fd).max()


@pytest.mark.parametrize("prm", [LayupParameter(1, 4, "thickness"), LayupParameter(2, 0, "thickness"),
                                 LayupParameter(1, 3, "angle")], ids=lambda p: p.label)
def test_section_mass_derivative_matches_central_differences(prm):
    h = 1e-7
    d = section_derivative(_box(), prm)
    fd = (SectionAnalysis.from_thin_walled(_box((prm, h))).mass_matrix()
          - SectionAnalysis.from_thin_walled(_box((prm, -h))).mass_matrix()) / (2 * h)
    np.testing.assert_allclose(d.mass_matrix, fd, atol=1e-6 * np.abs(fd).max() + 1e-12)


def test_angle_change_leaves_mass_unchanged():
    d = section_derivative(_box(), LayupParameter(1, 3, "angle"))
    np.testing.assert_array_equal(d.mass_matrix, 0.0)


def test_parser_round_trip(tmp_path):
    params = [LayupParameter(1, 4, "angle"), LayupParameter(3, 0, "thickness")]
    responses = [Response("displacement", 12, "UZ"), Response("compliance"), Response("frequency", mode=2)]
    p = tmp_path / "layup_sensitivity.txt"
    p.write_text(format_layup_sensitivity_file(params, responses))
    assert parse_layup_sensitivity_file(str(p)) == (params, responses)
    assert [r.label for r in responses] == ["UZ@12", "compliance", "f2"]
    assert params[0].label == "L1.p5.angle"
    assert parse_layup_sensitivity_file(str(tmp_path / "absent.txt")) == ([], [])


def test_parser_rejects_bad_rows(tmp_path):
    p = tmp_path / "layup_sensitivity.txt"
    p.write_text("[Layup_Parameter]\n1 0 angle\n")
    with pytest.raises(ValueError, match="count from 1"):
        parse_layup_sensitivity_file(str(p))
    with pytest.raises(ValueError):
        Response("frequency")
    with pytest.raises(ValueError):
        Response("displacement", node_id=3)


# ------------------------------------------------------------------ jobs

CFJ = common.Material(135e9, 10e9, 5e9, 0.3, 5e9, 3.8e9, 1550.0)
N_EL = 4
CAP = [45.0, -45.0, 0.0, 0.0, 10.0, 0.0, -45.0, 45.0]
WEB = [45.0, -45.0, -45.0, 45.0]
JOB_PARAMS = [LayupParameter(1, 4, "angle"), LayupParameter(1, 2, "thickness")]
STATIC = [Response("displacement", N_EL, "UZ"), Response("displacement", N_EL, "RX"), Response("compliance")]
EIGEN = [Response("frequency", mode=1), Response("frequency", mode=3)]
BUCKLING = [Response("buckling", mode=1), Response("buckling", mode=2)]
RESPONSES = {"static": STATIC, "eigen": EIGEN, "buckling": BUCKLING}
#: tip load of the buckling jobs: compression with shear, torque and bending, so every
#: second-order term of the prestress is present
BUCKLING_TIP = [-2e5, 4e3, 6e3, 3e3, 2e3, -1e3]


def _spec(change=None):
    cap, t = list(CAP), [0.5e-3] * len(CAP)
    if change is not None:
        (cap if change[0].kind == "angle" else t)[change[0].ply] += change[1]
    return common.SectionSpec({1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)},
                              [(1, 1, 2, 1), (2, 2, 3, 2), (3, 3, 4, 1), (4, 4, 1, 2)], {1: CFJ},
                              {1: [(1, a, tt) for a, tt in zip(cap, t)], 2: [(1, a, 0.5e-3) for a in WEB]})


def _taper_file(sec) -> str:
    """``section_taper.txt`` for every element: the box's depth growing by 10 % per metre and
    its width shrinking by 5 % per metre, cut either side of the section."""
    from pre_processing.parsing.section_taper_parser import format_section_taper_file
    from pre_processing.section_integration.thin_walled.taper import section_taper

    an = common.build_analysis(sec)
    d = 0.01

    def cut(s):
        pts = {k: (y * (1 - 0.05 * s), z * (1 + 0.10 * s)) for k, (y, z) in sec.points.items()}
        return ThinWalledSection(pts, [Wall(*w) for w in sec.walls])

    t = section_taper(cut(-d), cut(d), an.laminates, 2 * d)
    return format_section_taper_file({e: t for e in range(N_EL)})


def _run(name, sec, kind, element="LinearTimoshenkoBeamElement3D", precurvature=None, taper=False):
    from workflow_orchestrator.run_job import process_job

    tip = BUCKLING_TIP if kind == "buckling" else [0.0, 3e3, 20e3, 5e3, 0.0, 0.0]
    spec = common.CantileverSpec(4.0, N_EL, element, sec, tip, warping=int("Benscoter" in element))
    job = os.path.join(_WORK, "jobs", name)
    common.write_job(job, spec)
    common.write_supports(job, spec)
    if taper:
        with open(os.path.join(job, "section_taper.txt"), "w", encoding="utf-8") as f:
            f.write(_taper_file(sec))
    if precurvature is not None:
        with open(os.path.join(job, "precurvature.txt"), "w", encoding="utf-8") as f:
            f.write("[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n")
            f.writelines(f"{e} {precurvature[0]} {precurvature[1]} {precurvature[2]}\n" for e in range(N_EL))
    if kind == "eigen":
        with open(os.path.join(job, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write("[Simulation]\n[Type]\neigen\n\n[Solver]\ntype = direct\n\n[Eigen]\nnum_modes = 4\n")
    if kind == "buckling":
        with open(os.path.join(job, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write("[Simulation]\n[Type]\nbuckling\n\n[Solver]\ntype = direct\n\n[Buckling]\nenabled = true\n"
                    "num_modes = 4\nbuckling_prestress = linear_static\nbuckling_load_factor = 1.0\n")
    with open(os.path.join(job, "layup_sensitivity.txt"), "w", encoding="utf-8") as f:
        f.write(format_layup_sensitivity_file(JOB_PARAMS, RESPONSES[kind]))
    res = common.results_dir(os.path.join(_WORK, "results", name))
    logging.disable(logging.WARNING)
    try:
        process_job(job, res, {}, {}, force_serial=True, max_processes_per_job=1)
    finally:
        logging.disable(logging.NOTSET)
    with open(os.path.join(res, "secondary_results", "sensitivity", "layup_sensitivity.csv"), newline="",
              encoding="utf-8") as f:
        return {(r["parameter"], r["response"]): r for r in csv.DictReader(f)}


@pytest.mark.parametrize("prm, h", PARAMS, ids=[p.label for p, _ in PARAMS])
def test_benscoter_derivative_matches_central_differences(prm, h):
    """The 8 x 8 Benscoter stiffness and the warping inertia, both hoop options."""
    for hoop in ("free", "restrained"):
        d = section_derivative(_box(hoop=hoop), prm, benscoter=True)
        p, m = _box((prm, h), hoop), _box((prm, -h), hoop)
        fd = (p.benscoter_stiffness - m.benscoter_stiffness) / (2 * h)
        assert np.abs(d.benscoter_stiffness - fd).max() <= 1e-6 * np.abs(fd).max()
        assert d.warping_inertia == pytest.approx((p.warping_inertia - m.warping_inertia) / (2 * h), rel=1e-6)
        np.testing.assert_allclose(d.benscoter_stiffness[:6, :6], d.stiffness)


@pytest.mark.slow
@pytest.mark.parametrize("element", ["LinearTimoshenkoBeamElement3D", "LinearBenscoterBeamElement3D"])
@pytest.mark.parametrize("kind", ["static", "eigen"])
def test_job_sensitivities_match_central_differences_of_re_solved_jobs(kind, element):
    tag = f"{kind}_{element[6:10]}"
    base = _run(tag, _spec(), kind, element)
    assert len(base) == len(JOB_PARAMS) * len(RESPONSES[kind])
    for prm in JOB_PARAMS:
        h = 1e-3 if prm.kind == "angle" else 1e-8
        rp = _run(f"{tag}_p", _spec((prm, h)), kind, element)
        rm = _run(f"{tag}_m", _spec((prm, -h)), kind, element)
        for (label, resp), row in base.items():
            if label != prm.label:
                continue
            fd = (float(rp[label, resp]["response_value"]) - float(rm[label, resp]["response_value"])) / (2 * h)
            assert float(row["derivative"]) == pytest.approx(fd, rel=2e-5), (label, resp)
            if prm.kind == "thickness":
                rel = 0.5e-3 / float(row["response_value"]) * float(row["derivative"])
                assert float(row["relative_derivative"]) == pytest.approx(rel, rel=1e-12)


@pytest.mark.slow
@pytest.mark.parametrize("kind", ["static", "eigen"])
def test_pretwisted_curved_job_sensitivities_match_central_differences(kind):
    """A pretwisted, prebent cantilever: curved elements condensed from straight sub-elements,
    with the helical-fibre coupling of the twist rate."""
    curv = (0.15, 0.0, 0.05)
    tag = f"curved_{kind}"
    base = _run(tag, _spec(), kind, precurvature=curv)
    for prm in JOB_PARAMS:
        h = 1e-3 if prm.kind == "angle" else 1e-8
        rp = _run(f"{tag}_p", _spec((prm, h)), kind, precurvature=curv)
        rm = _run(f"{tag}_m", _spec((prm, -h)), kind, precurvature=curv)
        for (label, resp), row in base.items():
            if label != prm.label:
                continue
            fd = (float(rp[label, resp]["response_value"]) - float(rm[label, resp]["response_value"])) / (2 * h)
            assert float(row["derivative"]) == pytest.approx(fd, rel=2e-5), (label, resp)


@pytest.mark.slow
@pytest.mark.parametrize("kind", ["static", "eigen"])
def test_tapered_pretwisted_job_sensitivities_match_central_differences(kind):
    """Tapered (Resal coupling split by laminate) and pretwisted: the blade case."""
    curv = (0.15, 0.0, 0.0)
    tag = f"tapered_{kind}"
    base = _run(tag, _spec(), kind, precurvature=curv, taper=True)
    for prm in JOB_PARAMS:
        h = 1e-3 if prm.kind == "angle" else 1e-8
        rp = _run(f"{tag}_p", _spec((prm, h)), kind, precurvature=curv, taper=True)
        rm = _run(f"{tag}_m", _spec((prm, -h)), kind, precurvature=curv, taper=True)
        for (label, resp), row in base.items():
            if label != prm.label:
                continue
            fd = (float(rp[label, resp]["response_value"]) - float(rm[label, resp]["response_value"])) / (2 * h)
            assert float(row["derivative"]) == pytest.approx(fd, rel=2e-5), (label, resp)


def test_taper_split_by_laminate_sums_to_the_coupling(tmp_path):
    from pre_processing.parsing.section_taper_parser import parse_section_taper_file
    from pre_processing.section_integration.thin_walled.wall_stiffness import beam_wall_stiffness

    sec = _spec()
    p = tmp_path / "section_taper.txt"
    p.write_text(_taper_file(sec))
    t = parse_section_taper_file(str(p))[0]
    lams = common.build_analysis(sec).laminates
    total = sum(beam_wall_stiffness(lams[lid]).matrix[0, 0] * row for lid, row in t.parts.items())
    np.testing.assert_allclose(total, t.as_row(), rtol=1e-12, atol=1e-9 * np.abs(t.as_row()).max())


@pytest.mark.slow
def test_buckling_load_sensitivities_match_central_differences_of_re_solved_jobs():
    """Load factors of a cantilever under compression, shear, torque and bending."""
    base = _run("buckling", _spec(), "buckling")
    assert len(base) == len(JOB_PARAMS) * len(BUCKLING)
    for prm in JOB_PARAMS:
        h = 1e-3 if prm.kind == "angle" else 1e-8
        rp = _run("buckling_p", _spec((prm, h)), "buckling")
        rm = _run("buckling_m", _spec((prm, -h)), "buckling")
        for (label, resp), row in base.items():
            if label != prm.label:
                continue
            fd = (float(rp[label, resp]["response_value"]) - float(rm[label, resp]["response_value"])) / (2 * h)
            assert float(row["derivative"]) == pytest.approx(fd, rel=2e-5), (label, resp)


def test_unsupported_elements_are_refused():
    from post_processing.sensitivity.layup_sensitivity import element_derivatives

    class LinearEulerBernoulliBeamElement3D:
        element_id = 0

    class LinearTimoshenkoBeamElement3D:
        element_id = 1

        class _assigned_section_analysis:
            tapered = True

    with pytest.raises(NotImplementedError, match="LinearTimoshenkoBeamElement3D"):
        element_derivatives(LinearEulerBernoulliBeamElement3D(), None)
    with pytest.raises(NotImplementedError, match="tapered"):
        element_derivatives(LinearTimoshenkoBeamElement3D(), None)
