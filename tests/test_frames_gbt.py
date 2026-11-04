"""GBT beams off the global x axis: a turned GBT cantilever is the GBT cantilever, turned.

Its mode amplitudes (the extra nodal DOFs) belong to the section, so they map to
themselves; translations and rotations turn with the beam.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples", "thin_walled"))
import common  # noqa: E402

from pre_processing.element_library.element_frame import straight_frame  # noqa: E402

GBT = "LinearGBTBeamElement3D"


def _box_spec(load):
    b, h, t = 0.3, 0.12, 2e-3
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [(1, 1, 2, 1), (2, 2, 3, 1), (3, 3, 4, 1), (4, 4, 1, 1)]
    mat = {1: common.Material(70e9, 70e9, 26.9e9, 0.3, 26.9e9, 26.9e9, 2700.0)}
    sec = common.SectionSpec(pts, walls, mat, {1: [(1, 0.0, t)]}, n_sub=2)
    return common.CantileverSpec(1.5, 6, GBT, sec, tip_load=load, modes=2, warping=1)


def _turn_job(job_dir, d, twist, R, spec):
    """Rewrite the grid along ``d``, add the ``[twist]`` column and turn the tip load."""
    n, L = spec.n_el, spec.length
    p0 = np.array([0.4, -0.3, 0.2])
    pts = p0 + np.outer(np.linspace(0.0, L, n + 1), d)
    with open(os.path.join(job_dir, "grid.txt"), "w", encoding="utf-8") as f:
        f.write("[Grid]\n[node_id] [x] [y] [z]\n" + "".join(f"{i} {float(x)!r} {float(y)!r} {float(z)!r}\n" for i, (x, y, z) in enumerate(pts)))
    path = os.path.join(job_dir, "element.txt")
    lines = open(path, encoding="utf-8").read().splitlines()
    lines[1] += " [twist]"
    lines[2:] = [ln + f" {float(twist)!r}" for ln in lines[2:]]
    open(path, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    load = np.asarray(spec.tip_load, float)
    F = np.concatenate([R @ load[:3], R @ load[3:]])
    with open(os.path.join(job_dir, "point_load.txt"), "w", encoding="utf-8") as f:
        f.write("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                + " ".join(repr(float(v)) for v in (*pts[-1], *F)) + "\n")


def _run(job_dir, res_dir, n_nodes):
    import csv
    import logging

    from workflow_orchestrator.run_job import process_job

    res = common.results_dir(res_dir)
    logging.getLogger().setLevel(logging.WARNING)
    process_job(job_dir, res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return U.reshape(n_nodes, -1)


@pytest.mark.parametrize("d,twist", [([0.0, 1.0, 0.0], 0.0), ([1.0, 2.0, -2.0], 35.0)])
def test_a_turned_gbt_beam_is_the_beam_turned(tmp_path, d, twist):
    spec = _box_spec([200.0, -500.0, 800.0, 60.0, -40.0, 90.0])
    ref_dir, tur_dir = str(tmp_path / "ref"), str(tmp_path / "turned")
    for job in (ref_dir, tur_dir):
        common.write_job(job, spec)
        common.write_supports(job, spec)
    d = np.asarray(d, float) / np.linalg.norm(d)
    R = straight_frame([0, 0, 0], d, twist)
    _turn_job(tur_dir, d, twist, R, spec)
    U0 = _run(ref_dir, str(tmp_path / "ref_res"), spec.n_el + 1)
    U1 = _run(tur_dir, str(tmp_path / "tur_res"), spec.n_el + 1)
    assert U0.shape[1] > 6
    scale = np.abs(U0).max()
    np.testing.assert_allclose(U1[:, :3], U0[:, :3] @ R.T, atol=1e-9 * scale)
    np.testing.assert_allclose(U1[:, 3:6], U0[:, 3:6] @ R.T, atol=1e-9 * scale)
    np.testing.assert_allclose(U1[:, 6:], U0[:, 6:], atol=1e-9 * scale)       # mode amplitudes



def _bend_job(job_dir, spec, Rc, curved):
    """Make the cantilever a quarter circle of radius ``Rc`` in the x-y plane: curved
    elements (``k_z0 = 1/Rc``), or straight chords through points on the arc."""
    n = spec.n_el
    th = np.linspace(0.0, np.pi / 2, n + 1)
    pts = np.column_stack([Rc * np.sin(th), Rc * (1 - np.cos(th)), np.zeros(n + 1)])
    with open(os.path.join(job_dir, "grid.txt"), "w", encoding="utf-8") as f:
        f.write("[Grid]\n[node_id] [x] [y] [z]\n"
                + "".join(f"{i} {float(x)!r} {float(y)!r} {float(z)!r}\n" for i, (x, y, z) in enumerate(pts)))
    if curved:
        with open(os.path.join(job_dir, "precurvature.txt"), "w", encoding="utf-8") as f:
            f.write("[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n"
                    + "".join(f"{e} 0.0 0.0 {1.0 / Rc!r}\n" for e in range(n)))
    with open(os.path.join(job_dir, "point_load.txt"), "w", encoding="utf-8") as f:
        f.write("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                + " ".join(repr(float(v)) for v in (*pts[-1], *spec.tip_load)) + "\n")


def test_a_curved_gbt_box_is_its_chain_of_straight_elements(tmp_path):
    """Quarter-circle GBT box cantilever under an out-of-plane tip force and torque.

    A curved element is condensed from 8 straight sub-elements, so 8 curved elements equal
    64 straight chords through the same arc points, mode amplitudes included. The run
    also goes through post-processing, which must cope with elements that keep no Gauss data.
    """
    import dataclasses

    load = [0.0, 0.0, 800.0, 60.0, 0.0, 0.0]
    runs = {}
    for name, n, curved in (("arc", 8, True), ("fine", 64, False)):
        spec = dataclasses.replace(_box_spec(load), n_el=n)
        job, res = str(tmp_path / name), str(tmp_path / f"{name}_res")
        common.write_job(job, spec)
        common.write_supports(job, spec)
        _bend_job(job, spec, 1.0, curved)
        runs[name] = _run(job, res, n + 1)[-1]
        summary = [os.path.join(r, f) for r, _, fs in os.walk(res) for f in fs if f == "tertiary_summary.csv"]
        assert summary, f"{name}: the run did not finish post-processing"
    tip_a, tip_f = runs["arc"], runs["fine"]
    assert tip_a.size > 6 and abs(tip_f[3]) > 0
    np.testing.assert_allclose(tip_a, tip_f, rtol=1e-7, atol=1e-9 * np.abs(tip_f).max())
