"""Blade beam jobs written from an implicit-geometry-kernel loft (``blade_beam_job``).

Shell reference: composite-laminate-model 0.1, ``BladeBuilder.from_geometry`` on the same
lofts (curved Quad8 FSDT shells on the wall mid-surfaces), root clamped, tip loads through
a rigid tip diaphragm at the pitch axis. ``SHELL_TIP`` holds the tip reference-node
displacements, blade frame.
"""

from __future__ import annotations

import csv
import logging
import os

import numpy as np
import pytest

pytest.importorskip("implicit_geometry_kernel")

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402

from pre_processing.mesh_library.blade_beam_job import blade_beam_job, blade_to_job, job_to_blade  # noqa: E402
from pre_processing.section_integration.thin_walled import Ply, PlyMaterial, WallLaminate  # noqa: E402

GLASS = PlyMaterial(E1=38.6e9, E2=8.27e9, G12=4.14e9, nu12=0.26, G13=4.14e9, G23=3.063e9, rho=1800.0)
CARBON = PlyMaterial(E1=181e9, E2=10.3e9, G12=7.17e9, nu12=0.28, G13=7.17e9, G23=3.679e9, rho=1600.0)


def lam(mat, angles, t=0.5e-3):
    return WallLaminate([Ply(mat, a, t) for a in angles])


LAMS = {"skin": lam(GLASS, [0, 45, -45, -45, 45, 0]), "web": lam(GLASS, [45, -45, -45, 45]), "cap": lam(CARBON, [0] * 16)}
AF = Airfoil.naca4("0024")


def geometry(stations):
    return BladeGeometry(stations, "2D-F", 0.003, 0.008, 0.002, web_positions=(0.18, 0.48))


def run(job, tmp_path):
    from workflow_orchestrator.run_job import process_job

    res = str(tmp_path / (os.path.basename(job.job_dir) + "_results"))
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(os.path.join(res, sub), exist_ok=True)
    logging.getLogger().setLevel(logging.WARNING)
    process_job(job.job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return U.reshape(job.nodes.shape[0], -1)


def test_axes():
    v = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    np.testing.assert_array_equal(blade_to_job(v), [3, 1, 2, 6, 4, 5])
    np.testing.assert_array_equal(job_to_blade(blade_to_job(v)), v)
    # right handed: blade x cross y = z maps to job Y cross Z = X
    np.testing.assert_array_equal(np.cross(blade_to_job([1, 0, 0]), blade_to_job([0, 1, 0])), blade_to_job([0, 0, 1]))


def test_straight_loft_matches_a_hand_written_job(tmp_path):
    """A prismatic, untwisted loft gives the job one would write by hand from its section."""
    from pre_processing.parsing.laminate_parser import parse_laminate_file
    from pre_processing.parsing.section_geometry_parser import parse_section_geometry_file
    from pre_processing.section_integration.thin_walled.from_geometry import thin_walled_from_cut
    from pre_processing.section_integration.thin_walled.section_analysis import ThinWalledSectionAnalysis

    geo = geometry([Station(0.0, 1.0, AF, pitch_axis=0.3), Station(4.0, 1.0, AF, pitch_axis=0.3)])
    F = 1000.0
    job = blade_beam_job(str(tmp_path / "blade"), geo, LAMS, 4, loads={4.0: blade_to_job([0.0, F, 0, 0, 0, 0])})
    assert not os.path.exists(os.path.join(job.job_dir, "precurvature.txt"))
    np.testing.assert_allclose(job.nodes, np.column_stack([np.linspace(0, 4, 5), np.zeros(5), np.zeros(5)]))
    # the files hold the cut's section exactly
    direct = thin_walled_from_cut(geo.section_at(0.5), LAMS)
    lams = parse_laminate_file(os.path.join(job.job_dir, "laminate.txt"))
    secs = parse_section_geometry_file(os.path.join(job.job_dir, "section_geometry.txt"), lams.keys())
    assert len(secs) == 4
    K_file = ThinWalledSectionAnalysis(secs[1], lams).classical_stiffness
    np.testing.assert_allclose(K_file, direct.analysis().classical_stiffness, rtol=1e-12, atol=1e-6)
    # the hand-written job: examples/thin_walled/common.py's cantilever with the same section
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples", "thin_walled"))
    import common

    mats, mat_id, plies = {}, {}, {}
    for lid, lam_ in direct.laminates.items():
        rows = []
        for p in lam_.plies:
            if p.material not in mat_id:
                mat_id[p.material] = len(mat_id) + 1
                m = p.material
                mats[mat_id[p.material]] = common.Material(m.E1, m.E2, m.G12, m.nu12, m.G13, m.G23, m.rho)
            rows.append((mat_id[p.material], p.angle_deg, p.thickness))
        plies[lid] = rows
    sec = direct.section
    spec = common.CantileverSpec(4.0, 4, "LinearTimoshenkoBeamElement3D", common.SectionSpec(
        {pid: tuple(p) for pid, p in sec.points.items()}, [(w.wall_id, w.start, w.end, w.laminate_id) for w in sec.walls],
        mats, plies), tip_load=[0, 0, F, 0, 0, 0])
    hand = common.run(spec, str(tmp_path / "hand"), "hand").U
    U = run(job, tmp_path)
    # common.py writes 9 significant digits
    np.testing.assert_allclose(U, hand, rtol=1e-7, atol=1e-12 * np.abs(hand).max())
    assert U[-1, 2] > 0.0                         # flap load lifts the tip


def test_prebent_swept_nodes_follow_the_offsets(tmp_path):
    st = [Station(0.0, 1.0, AF, pitch_axis=0.3), Station(2.0, 1.0, AF, pitch_axis=0.3, offset=(0.05, 0.02)),
          Station(4.0, 1.0, AF, pitch_axis=0.3, offset=(0.15, 0.10))]
    job = blade_beam_job(str(tmp_path / "bent"), geometry(st), LAMS, 4)
    np.testing.assert_allclose(job.nodes[:, 0], [0, 1, 2, 3, 4])
    np.testing.assert_allclose(job.nodes[:, 1], [0, 0.025, 0.05, 0.10, 0.15], atol=1e-12)   # sweep: blade x -> Y
    np.testing.assert_allclose(job.nodes[:, 2], [0, 0.01, 0.02, 0.06, 0.10], atol=1e-12)    # prebend: blade y -> Z
    with open(os.path.join(job.job_dir, "grid.txt"), encoding="utf-8") as f:
        rows = [ln.split() for ln in f.read().splitlines()[2:]]
    np.testing.assert_allclose(np.array(rows, float)[:, 1:], job.nodes, rtol=0, atol=0)
    # sections ride on the pitch axis: the same section as a straight blade's
    straight = blade_beam_job(str(tmp_path / "straight"), geometry([Station(0.0, 1.0, AF, pitch_axis=0.3),
                                                                     Station(4.0, 1.0, AF, pitch_axis=0.3)]), LAMS, 4)
    for sid in (1, 4):
        np.testing.assert_allclose(job.sections[sid].analysis().classical_stiffness,
                                   straight.sections[sid].analysis().classical_stiffness, rtol=1e-9, atol=1e-3)


def test_pretwist_rate_from_the_twist_gradient(tmp_path):
    st = [Station(0.0, 1.0, AF, twist=10.0, pitch_axis=0.3), Station(5.0, 1.0, AF, twist=0.0, pitch_axis=0.3)]
    job = blade_beam_job(str(tmp_path / "twist"), geometry(st), LAMS, 5)
    # nose-up twist falls 2 degrees per metre: the section turns +2 degrees per metre about the span
    np.testing.assert_allclose(job.pretwist, np.radians(2.0), rtol=1e-12)
    from pre_processing.parsing.precurvature_parser import parse_precurvature

    k = parse_precurvature(os.path.join(job.job_dir, "precurvature.txt"), np.arange(5))
    np.testing.assert_allclose(k[:, 0], np.radians(2.0), rtol=1e-12)
    # the arc's frame at an element's end turns the section as the geometry does
    from pre_processing.element_library.element_frame import arc_geometry

    arc = arc_geometry(job.nodes[0], job.nodes[1], k[0])
    R = arc.R_mid.T @ arc.R1
    turn = np.degrees(np.arctan2(R[2, 1], R[1, 1]))
    assert turn == pytest.approx(-(job.twist[1] - (job.twist[0] + job.twist[1]) / 2), abs=1e-9)


# composite-laminate-model shell (12 x 1260 curved Quad8; 24 spans moves these < 0.1 %),
# blade frame (x, y, z, rx, ry, rz) of the tip reference node under 1 kN or 1 kN m at the tip.
SHELL_TIP = {
    # prismatic, untwisted (chord 1 m, 6 m)
    "prismatic": {"flap": [6.838e-9, 6.764595e-3, 0, 0, 0, -1.750554e-5],
                  "edge": [3.453640e-3, 6.837e-9, 4.222867e-5, 0, 0, 0],
                  "torque": [0, -1.750553e-5, 0, 0, 0, 3.124483e-3]},
    # prismatic, twisted 8 degrees nose up at the root to 0 at the tip
    "twisted": {"flap": [6.648177e-4, 6.701384e-3, 0, 0, 0, 1.724417e-5],
                "edge": [3.582891e-3, 6.648177e-4, 4.241060e-5, 0, 0, 6.534208e-5]},
    # tapered (chord 1.2 m to 0.8 m) and twisted 8 degrees to 0
    "tapered": {"flap": [4.693385e-4, 5.112116e-3, 0, 0, 0, -8.171013e-6],
                "edge": [2.787707e-3, 4.693385e-4, 3.812701e-5, 0, 0, 5.692878e-5]},
}
WINGS = {
    "prismatic": [Station(0.0, 1.0, AF, pitch_axis=0.3), Station(6.0, 1.0, AF, pitch_axis=0.3)],
    "twisted": [Station(0.0, 1.0, AF, twist=8.0, pitch_axis=0.3), Station(6.0, 1.0, AF, pitch_axis=0.3)],
    "tapered": [Station(0.0, 1.2, AF, twist=8.0, pitch_axis=0.3), Station(6.0, 0.8, AF, pitch_axis=0.3)],
}
LOADS = {"flap": [0, 1000.0, 0, 0, 0, 0], "edge": [1000.0, 0, 0, 0, 0, 0], "torque": [0, 0, 0, 0, 0, 1000.0]}


def _tip(wing, case, tmp_path, n=8, **kw):
    job = blade_beam_job(str(tmp_path / f"{wing}_{case}"), geometry(WINGS[wing]), LAMS, n,
                         loads={6.0: blade_to_job(LOADS[case])}, **kw)
    return job_to_blade(run(job, tmp_path)[-1, :6])


@pytest.mark.slow
@pytest.mark.parametrize("case", ["flap", "edge", "torque"])
def test_prismatic_wing_matches_the_shell(case, tmp_path):
    tip, ref = _tip("prismatic", case, tmp_path), SHELL_TIP["prismatic"][case]
    k = {"flap": 1, "edge": 0, "torque": 5}[case]
    assert tip[k] == pytest.approx(ref[k], rel=0.01)


@pytest.mark.slow
def test_twisted_wing_matches_the_shell(tmp_path):
    """The pretwist rate gives the edge-twist coupling of the twisted wing; without it there is none."""
    flap, edge = _tip("twisted", "flap", tmp_path), _tip("twisted", "edge", tmp_path)
    ref = SHELL_TIP["twisted"]
    assert flap[1] == pytest.approx(ref["flap"][1], rel=0.01)
    assert flap[0] == pytest.approx(ref["flap"][0], rel=0.06)          # edge motion under flap load
    assert edge[0] == pytest.approx(ref["edge"][0], rel=0.01)
    assert edge[5] == pytest.approx(ref["edge"][5], rel=0.03)          # twist under edge load
    plain = _tip("twisted", "edge", tmp_path / "plain", pretwist=False)
    assert abs(plain[5]) < 0.05 * abs(ref["edge"][5])


@pytest.mark.slow
def test_tapered_twisted_wing_against_the_shell(tmp_path):
    """Edgewise within 1 %. Flapwise 3.7 % soft: the beam's shear flexibility ignores the
    share of the shear the tapered (inclined) caps carry, which the shell has."""
    flap, edge = _tip("tapered", "flap", tmp_path), _tip("tapered", "edge", tmp_path)
    ref = SHELL_TIP["tapered"]
    assert edge[0] == pytest.approx(ref["edge"][0], rel=0.01)
    assert edge[5] == pytest.approx(ref["edge"][5], rel=0.06)
    assert 1.0 < flap[1] / ref["flap"][1] < 1.045


def test_line_load_on_a_bent_blade(tmp_path):
    """A flapwise line load: a blade bent by a micrometre (rotated elements, line loads spread
    along them) gives the straight blade's tip deflection."""
    straight = [Station(0.0, 1.0, AF, pitch_axis=0.3), Station(4.0, 1.0, AF, pitch_axis=0.3)]
    bent = [Station(0.0, 1.0, AF, pitch_axis=0.3), Station(4.0, 1.0, AF, pitch_axis=0.3, offset=(0.0, 1e-6))]
    q = {0.0: blade_to_job([0, 400.0, 0, 0, 0, 0]), 4.0: blade_to_job([0, 100.0, 0, 0, 0, 0])}
    tips = []
    for name, st in (("straight", straight), ("bent", bent)):
        job = blade_beam_job(str(tmp_path / name), geometry(st), LAMS, 4, line_loads=q)
        tips.append(run(job, tmp_path)[-1, :6])
    assert tips[0][2] > 0.0
    np.testing.assert_allclose(tips[1], tips[0], rtol=1e-4, atol=1e-4 * abs(tips[0][2]))
    with pytest.raises(ValueError, match="cover the blade"):
        blade_beam_job(str(tmp_path / "short"), geometry(straight), LAMS, 4, line_loads={0.0: [0] * 6, 3.0: [0] * 6})


def test_two_d_sections(tmp_path):
    """``kind="2d"`` writes layered 2D section meshes; the stiffness is the thin-walled one's within 1.5 %."""
    from pre_processing.section_integration.solid.mesh import parse_section_mesh_2d

    geo = geometry([Station(0.0, 1.0, AF, pitch_axis=0.3), Station(2.0, 1.0, AF, pitch_axis=0.3)])
    job = blade_beam_job(str(tmp_path / "solid"), geo, LAMS, 1, kind="2d", mesh_options={"element_size": 0.02})
    meshes = parse_section_mesh_2d(os.path.join(job.job_dir, "section_mesh_2d.txt"), {})
    assert set(meshes) == {1}
    from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis

    K2 = SolidSectionAnalysis(meshes[1]).classical_stiffness
    Kt = blade_beam_job(str(tmp_path / "tw"), geo, LAMS, 1).sections[1].analysis().classical_stiffness
    for i in (0, 1, 2):
        assert K2[i, i] == pytest.approx(Kt[i, i], rel=0.015)


def test_rejects_types_without_frames(tmp_path):
    geo = geometry([Station(0.0, 1.0, AF), Station(1.0, 1.0, AF)])
    with pytest.raises(ValueError, match="element_type"):
        blade_beam_job(str(tmp_path / "x"), geo, LAMS, 1, element_type="LinearGBTBeamElement3D")
    with pytest.raises(ValueError, match="no node"):
        blade_beam_job(str(tmp_path / "x"), geo, LAMS, 2, loads={0.3: [0, 0, 1.0, 0, 0, 0]})
