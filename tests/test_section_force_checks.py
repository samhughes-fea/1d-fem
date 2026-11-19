"""Ply checks under section forces given from outside the job (``section_forces.txt``)."""

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

pytest.importorskip("anisotropic_section_model.section_strength")

from examples.thin_walled import common  # noqa: E402
from onedfem.pre_processing.parsing.ply_strength_parser import format_ply_strength_file  # noqa: E402
from onedfem.pre_processing.parsing.section_forces_parser import (  # noqa: E402
    format_section_forces_file,
    parse_section_forces_file,
    to_engine_forces,
)
from anisotropic_section_model import axes  # noqa: E402
from anisotropic_section_model.ply_failure import PlyStrength  # noqa: E402

CFRP = common.Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.8e9, rho=1550.0)
STRENGTH = PlyStrength(Xt=1500e6, Xc=1200e6, Yt=50e6, Yc=200e6, S12=70e6, p_TL_t=0.35, p_TL_c=0.30)
B, H, L, N_EL = 0.6, 0.3, 6.0, 2
TIP = np.array([2e3, 3e3, 20e3, 500.0, 200.0, -300.0])     # [F_x, F_y, F_z, M_x, M_y, M_z] at x = L


def _section():
    cap = [(1, a, 0.5e-3) for a in (45, -45, 0, 0, 0, 0, -45, 45)]
    web = [(1, a, 0.5e-3) for a in (45, -30, -45, 60)]
    return common.SectionSpec(
        points={1: (-B / 2, -H / 2), 2: (B / 2, -H / 2), 3: (B / 2, H / 2), 4: (-B / 2, H / 2)},
        walls=[(1, 1, 2, 1), (2, 2, 3, 2), (3, 3, 4, 1), (4, 4, 1, 2)], materials={1: CFRP},
        laminates={1: cap, 2: web})


def cut_forces(x):
    """Right-handed ``[N, M_y, M_z, V_y, V_z, T]`` at ``x`` of the cantilever under ``TIP``."""
    Fx, Fy, Fz, Mx, My, Mz = TIP
    a = L - x
    return np.array([Fx, My - a * Fz, Mz + a * Fy, Fy, Fz, Mx])


def _read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _run(tmp_path, name, force_rows, load=TIP, kind="thin_walled"):
    from onedfem.workflow_orchestrator.run_job import process_job

    spec = common.CantileverSpec(L, N_EL, "LinearTimoshenkoBeamElement3D", _section(), list(load))
    job = os.path.join(str(tmp_path), "jobs", name)
    common.write_job(job, spec)
    common.write_supports(job, spec)
    if kind == "2d":
        from onedfem.interfaces.blade.solid_sections import write_section_mesh_2d
        from anisotropic_section_model.solid.mesh import SectionMesh

        an = common.build_analysis(spec.section)
        mesh = SectionMesh.from_thin_walled(an.section, an.laminates, elements_along=6)
        write_section_mesh_2d(os.path.join(job, "section_mesh_2d.txt"), {1: mesh})
        os.remove(os.path.join(job, "section_geometry.txt"))
    with open(os.path.join(job, "ply_strength.txt"), "w", encoding="utf-8") as f:
        f.write(format_ply_strength_file({1: STRENGTH}))
    with open(os.path.join(job, "section_forces.txt"), "w", encoding="utf-8") as f:
        f.write(format_section_forces_file(force_rows))
    res = common.results_dir(os.path.join(str(tmp_path), "results", name))
    logging.getLogger().setLevel(logging.WARNING)
    process_job(job, res, {}, {}, force_serial=True, max_processes_per_job=1)
    return res


def test_file_round_trip_and_frames(tmp_path):
    f = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    rows = [("a", 0, 0.5, "engine", f), ("a", 1, 1.5, "public", axes.from_engine(f)),
            ("b", 1, 1.5, "beamdyn", axes.MIRROR @ axes.from_engine(f))]
    p = tmp_path / "section_forces.txt"
    p.write_text(format_section_forces_file(rows), encoding="utf-8")
    back = parse_section_forces_file(str(p))
    assert [r[:4] for r in back] == [r[:4] for r in rows]
    for case, eid, x, frame, g in back:
        np.testing.assert_allclose(to_engine_forces(g, frame), f)
    p.write_text(format_section_forces_file(rows + [rows[0]]), encoding="utf-8")
    with pytest.raises(ValueError, match="given twice"):
        parse_section_forces_file(str(p))


def _solved_gauss_forces(res, eid):
    """Right-handed section forces at the Gauss points of a solved element (the tertiary output
    carries M_y with the element sign)."""
    path = os.path.join(res, "tertiary_results", "gaussian", "section_forces", f"section_forces_elem_{eid:06d}.csv")
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(ln for ln in f if not ln.startswith("#")))
    return [np.array([float(r["N"]), -float(r["My"]), float(r["Mz"]), float(r["Vy"]), float(r["Vz"]), float(r["T"])])
            for r in rows]


@pytest.mark.slow
@pytest.mark.parametrize("kind", ["thin_walled", "2d"])
def test_external_forces_give_the_solved_ply_failure(tmp_path, kind):
    """The solved cantilever's own section forces, given back as a table in each frame, reproduce
    its ``ply_failure.csv`` at every Gauss point."""
    first = _run(tmp_path, "probe", [("probe", 0, 0.0, "engine", np.zeros(6))], kind=kind)
    gauss = _read(os.path.join(first, "secondary_results", "ply_failure", "ply_failure.csv"))
    solved: dict = {}
    xs: dict = {}
    for r in gauss:
        eid, g = int(r["element_id"]), int(r["gauss_point"])
        xs[(eid, g)] = float(r["x"])
        solved[(eid, g)] = max(solved.get((eid, g), 0.0), float(r["exposure"]))
    rows = []
    for k, (eid, g) in enumerate(sorted(solved)):
        f = _solved_gauss_forces(first, eid)[g]
        frame = ("engine", "public", "beamdyn")[k % 3]
        f = {"engine": f, "public": axes.from_engine(f), "beamdyn": axes.MIRROR @ axes.from_engine(f)}[frame]
        rows.append(("solved", eid, xs[(eid, g)], frame, f))
    res = _run(tmp_path, "cut", rows, kind=kind)
    checks = _read(os.path.join(res, "secondary_results", "section_force_checks", "section_force_checks.csv"))
    got = {(int(r["element_id"]), float(r["x"])): float(r["exposure"]) for r in checks if r["criterion"] == "exposure"}
    assert len(got) == len(solved)
    for key, f in solved.items():
        assert got[(key[0], xs[key])] == pytest.approx(f, rel=1e-6), key


def test_envelope_over_load_cases(tmp_path):
    f = cut_forces(1.5)
    rows = [("light", 0, 1.5, "engine", f), ("heavy", 0, 1.5, "engine", 2.0 * f),
            ("reversed", 0, 1.5, "engine", -1.5 * f), ("light", 1, 4.5, "engine", cut_forces(4.5))]
    res = _run(tmp_path, "env", rows)
    out = os.path.join(res, "secondary_results", "section_force_checks")
    checks = _read(os.path.join(out, "section_force_checks.csv"))
    env = _read(os.path.join(out, "section_force_envelope.csv"))
    assert len(checks) == 4 * 5 and len(env) == 2 * 5
    ms = {r["load_case"]: float(r["exposure"]) for r in checks if r["criterion"] == "max_stress"
          and r["element_id"] == "0"}
    assert ms["heavy"] == pytest.approx(2.0 * ms["light"], rel=1e-10)   # max stress is linear in the load
    for r in env:
        if r["element_id"] != "0":
            assert r["load_case"] == "light"
            continue
        cases = [float(c["exposure"]) for c in checks if c["element_id"] == "0" and c["criterion"] == r["criterion"]]
        assert float(r["exposure"]) == max(cases)
        assert r["load_case"] in ("heavy", "reversed")
        assert float(r["margin"]) == pytest.approx(1.0 / float(r["exposure"]) - 1.0)


def test_unknown_elements_are_refused(tmp_path):
    from types import SimpleNamespace

    from onedfem.post_processing.failure.section_force_checks import write_section_force_checks

    (tmp_path / "section_forces.txt").write_text(format_section_forces_file([("a", 99, 0.0, "engine", np.ones(6))]),
                                                 encoding="utf-8")
    elements = [SimpleNamespace(element_id=0)]
    with pytest.raises(ValueError, match="needs ply_strength.txt"):
        write_section_force_checks(str(tmp_path), str(tmp_path), elements, None, None)
    (tmp_path / "ply_strength.txt").write_text(format_ply_strength_file({1: STRENGTH}), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown elements"):
        write_section_force_checks(str(tmp_path), str(tmp_path), elements, None, None)
    (tmp_path / "section_forces.txt").write_text(format_section_forces_file([("a", 0, 0.0, "engine", np.ones(6))]),
                                                 encoding="utf-8")
    with pytest.raises(ValueError, match="no thin-walled or 2D section"):
        write_section_force_checks(str(tmp_path), str(tmp_path), elements, None, None)
