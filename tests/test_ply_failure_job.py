"""Ply failure at the Gauss points of a solved job (``ply_strength.txt`` -> ``ply_failure.csv``)."""

from __future__ import annotations

import csv
import logging
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from examples.thin_walled import common  # noqa: E402
from onedfem.pre_processing.parsing.ply_strength_parser import format_ply_strength_file, parse_ply_strength_file  # noqa: E402
from onedfem.pre_processing.section_integration.ply_failure import CRITERIA, PlyStrength, exposures  # noqa: E402

CFRP = common.Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
STRENGTH = PlyStrength(Xt=1500e6, Xc=1200e6, Yt=50e6, Yc=200e6, S12=70e6, p_TL_t=0.35, p_TL_c=0.30)
B, H, L, N_EL, P = 0.6, 0.3, 6.0, 2, 20e3


def _section():
    cap = [(1, a, 0.5e-3) for a in (45, -45, 0, 0, 0, 0, -45, 45)]
    web = [(1, a, 0.5e-3) for a in (45, -45, -45, 45)]
    return common.SectionSpec(
        points={1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 2), (3, 3, 4, 1), (4, 4, 1, 2)], materials={1: CFRP},
        laminates={1: cap, 2: web})


def _run(tmp_path, name, load, kind="thin_walled"):
    from onedfem.workflow_orchestrator.run_job import process_job

    spec = common.CantileverSpec(L, N_EL, "LinearTimoshenkoBeamElement3D", _section(), load)
    job = os.path.join(str(tmp_path), "jobs", name)
    common.write_job(job, spec)
    common.write_supports(job, spec)
    if kind == "2d":
        from onedfem.interfaces.blade.solid_sections import write_section_mesh_2d
        from onedfem.pre_processing.section_integration.solid.mesh import SectionMesh

        an = common.build_analysis(spec.section)
        mesh = SectionMesh.from_thin_walled(an.section, an.laminates, elements_along=6)
        write_section_mesh_2d(os.path.join(job, "section_mesh_2d.txt"), {1: mesh})
        os.remove(os.path.join(job, "section_geometry.txt"))
    with open(os.path.join(job, "ply_strength.txt"), "w", encoding="utf-8") as f:
        f.write(format_ply_strength_file({1: STRENGTH}))
    res = common.results_dir(os.path.join(str(tmp_path), "results", name))
    logging.getLogger().setLevel(logging.WARNING)
    process_job(job, res, {}, {}, force_serial=True, max_processes_per_job=1)
    out = os.path.join(res, "secondary_results", "ply_failure")
    with open(os.path.join(out, "ply_failure.csv"), newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    with open(os.path.join(out, "ply_failure_summary.csv"), newline="", encoding="utf-8") as f:
        summary = list(csv.DictReader(f))
    return rows, summary


def test_strength_file_round_trip(tmp_path):
    p = tmp_path / "ply_strength.txt"
    p.write_text(format_ply_strength_file({1: STRENGTH, 4: PlyStrength(1, 2, 3, 4, 5, S23=6.0)}), encoding="utf-8")
    back = parse_ply_strength_file(str(p))
    assert back[1] == STRENGTH and back[4].S23 == 6.0 and back[4].Zt is None
    p.write_text("[Ply_Strength]\n[material_id] [Xt] [Xc] [Yt] [Yc] [S12] [p_TT]\n2 1 1 1 1 1 -\n", encoding="utf-8")
    assert parse_ply_strength_file(str(p))[2].p_TT == 0.25
    p.write_text("[Ply_Strength]\n[material_id] [Xt] [Xc] [Yt] [Yc] [S12] [bad]\n2 1 1 1 1 1 1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unknown columns"):
        parse_ply_strength_file(str(p))
    assert parse_ply_strength_file(str(tmp_path / "missing.txt")) == {}


@pytest.mark.slow
def test_thin_walled_ply_failure_at_gauss_points(tmp_path):
    rows, summary = _run(tmp_path, "p1", [0, 0, P, 0, 0, 0])
    # every element is checked at each Gauss point of its stiffness integration
    gp = {(r["element_id"], r["gauss_point"]) for r in rows}
    n_g = len(gp) // N_EL
    assert len(gp) == N_EL * n_g and n_g >= 2
    xis = sorted({float(r["xi"]) for r in rows})
    np.testing.assert_allclose(xis, np.polynomial.legendre.leggauss(n_g)[0], atol=1e-12)
    # the written exposures are the criteria of the written stresses
    S = np.array([[float(r[k]) for k in ("sigma_11", "sigma_22", "sigma_33", "tau_23", "tau_13", "tau_12")]
                  for r in rows])
    e = exposures(S, STRENGTH)
    for c in CRITERIA:
        np.testing.assert_allclose([float(r[c]) for r in rows], e[c], rtol=1e-12, atol=1e-15)
    f = np.array([float(r["exposure"]) for r in rows])
    np.testing.assert_allclose([float(r["margin"]) for r in rows], 1 / f - 1, rtol=1e-12)
    # the 0-degree cap plies at the root carry the bending: tension in the bottom cap,
    # compression in the top cap, about E_11 kappa H / 2 with the section's bending compliance.
    # The two-node element's curvature is constant: that of the moment at its mid-point.
    an = common.build_analysis(_section())
    for wall, sign in (("1", 1.0), ("3", -1.0)):
        root = [r for r in rows if r["element_id"] == "0" and r["gauss_point"] == "0" and r["wall_id"] == wall
                and r["ply"] in ("2", "3", "4", "5")]
        s11 = np.array([float(r["sigma_11"]) for r in root])
        assert np.all(sign * s11 > 0)
        kappa = np.linalg.solve(an.classical_stiffness, [0.0, P * (L - L / (2 * N_EL)), 0.0, 0.0])[1]
        assert np.median(np.abs(s11)) == pytest.approx(CFRP.E1 * abs(kappa) * H / 2, rel=0.02)
    # the summary holds each element's worst point per criterion
    worst = {(s["element_id"], s["criterion"]): float(s["exposure"]) for s in summary}
    for el in ("0", "1"):
        for c in CRITERIA:
            assert worst[(el, c)] == pytest.approx(max(float(r[c]) for r in rows if r["element_id"] == el), rel=1e-12)
    assert worst[("0", "max_stress")] > worst[("1", "max_stress")]     # the root is worse
    # exposures scale with the load
    rows2, _ = _run(tmp_path, "p2", [0, 0, 2 * P, 0, 0, 0])
    np.testing.assert_allclose([float(r["exposure"]) for r in rows2], 2 * f, rtol=1e-9)


@pytest.mark.slow
def test_2d_section_ply_failure_matches_thin_walled(tmp_path):
    """The same box as a 2D mesh: the root's worst fibre exposure within 5 % of the thin-walled one."""
    rows_t, sum_t = _run(tmp_path, "tw", [0, 0, P, 0, 0, 0])
    rows_s, sum_s = _run(tmp_path, "solid", [0, 0, P, 0, 0, 0], kind="2d")
    assert all(r["section_element"] != "" and r["material_id"] == "1" for r in rows_s)
    S = np.array([[float(r[k]) for k in ("sigma_11", "sigma_22", "sigma_33", "tau_23", "tau_13", "tau_12")]
                  for r in rows_s])
    np.testing.assert_allclose([float(r["puck"]) for r in rows_s], exposures(S, STRENGTH)["puck"], rtol=1e-12)
    worst = lambda summ: {s["criterion"]: float(s["exposure"]) for s in summ if s["element_id"] == "0"}  # noqa: E731
    wt, ws = worst(sum_t), worst(sum_s)
    assert ws["max_stress"] == pytest.approx(wt["max_stress"], rel=0.05)
