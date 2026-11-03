"""GBT on 2D section meshes: the same thin box, as a 2D mesh and as a thin-walled section.

With the cross-section modes held at zero both are the same Benscoter-type beam; with
modes free, the 2D section's distortional modes (in-plane eigenmodes of the meshed
walls) catch the local flange modes of the box, which the thin-walled GBT basis
(a few nodes per wall) is too coarse to show. The reference is a CalculiX C3D20R
model of the same box (``post_processing/validation_visualisers/gbt_2d_box``).
"""

from __future__ import annotations

import csv
import json
import logging
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "examples", "thin_walled"))
import common  # noqa: E402

GBT = "LinearGBTBeamElement3D"
N_EL = 16
REFERENCE = os.path.join(HERE, "..", "post_processing", "validation_visualisers", "gbt_2d_box", "reference.json")


def _box_section():
    b, h, t = 0.3, 0.12, 2e-3
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [(1, 1, 2, 1), (2, 2, 3, 1), (3, 3, 4, 1), (4, 4, 1, 1)]
    mat = {1: common.Material(70e9, 70e9, 26.9e9, 0.3, 26.9e9, 26.9e9, 2700.0)}
    return common.SectionSpec(pts, walls, mat, {1: [(1, 0.0, t)]}, n_sub=4)


def box_mesh():
    """The box as a 2D section mesh (Q8 elements, 12 along each wall)."""
    from pre_processing.section_integration.solid.mesh import SectionMesh

    tw = common.build_analysis(_box_section())
    return SectionMesh.from_thin_walled(tw.section, tw.laminates, elements_along=12)


def _run(tmp_path, kind, modes, eigen=False):
    from pre_processing.section_integration.solid.from_geometry import write_section_mesh_2d
    from workflow_orchestrator.run_job import process_job

    load = [100.0, 300.0, 500.0, 200.0, 50.0, 80.0]
    spec = common.CantileverSpec(1.5, N_EL, GBT, _box_section(), tip_load=load, modes=modes, warping=1)
    name = f"{kind}_{modes}_{'eigen' if eigen else 'static'}"
    jd = str(tmp_path / "jobs" / name)
    common.write_job(jd, spec)
    common.write_supports(jd, spec)
    if kind == "2d":
        os.remove(os.path.join(jd, "section_geometry.txt"))
        os.remove(os.path.join(jd, "laminate.txt"))
        write_section_mesh_2d(os.path.join(jd, "section_mesh_2d.txt"), {1: box_mesh()})
    if eigen:
        with open(os.path.join(jd, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write("[Simulation]\n[Type]\neigen\n\n[Eigen]\nnum_modes = 8\n")
        with open(os.path.join(jd, "material.txt"), "w", encoding="utf-8") as f:
            f.write("[Material]\n[element_id] [E] [G] [nu] [rho]\n"
                    + "".join(f"{e} 7e10 2.7e10 0.3 2700\n" for e in range(N_EL)))
    res = common.results_dir(str(tmp_path / "res" / name))
    logging.getLogger().setLevel(logging.WARNING)
    process_job(jd, res, {}, {}, force_serial=True, max_processes_per_job=1)
    if eigen:
        path = os.path.join(res, "primary_results", "modal_results", f"{name}_frequencies.txt")
        return np.loadtxt(path)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return U.reshape(N_EL + 1, -1)


def test_section_mesh_2d_file_keeps_junction_ties(tmp_path):
    from pre_processing.section_integration.solid.from_geometry import write_section_mesh_2d
    from pre_processing.section_integration.solid.mesh import parse_section_mesh_2d

    mesh = box_mesh()
    assert mesh.ties
    write_section_mesh_2d(str(tmp_path / "m.txt"), {1: mesh})
    back = parse_section_mesh_2d(str(tmp_path / "m.txt"), {})[1]
    assert back.ties.keys() == mesh.ties.keys()
    for n, row in mesh.ties.items():
        np.testing.assert_allclose([c for _, c in back.ties[n]], [c for _, c in row])


@pytest.mark.parametrize("modes", [0, 6])
def test_2d_gbt_statics_match_thin_walled_gbt(tmp_path, modes):
    """Tip translations and rotations of the two section models agree within 3 %
    (the 2D walls overlap at the tied corners; the twist differs most)."""
    tw = _run(tmp_path, "tw", modes)[-1, :6]
    d2 = _run(tmp_path, "2d", modes)[-1, :6]
    np.testing.assert_allclose(d2, tw, rtol=3e-2, atol=1e-2 * np.abs(tw).max())


def test_2d_gbt_beam_frequencies_match_calculix_box(tmp_path):
    """Without distortional modes, the first flap and edge frequencies are within 3 % and
    5 % of the CalculiX box (64.0 Hz flap, 135.9 Hz edge; its second mode, 90.4 Hz, is a
    torsion-distortion mode)."""
    with open(REFERENCE, encoding="utf-8") as f:
        ref = json.load(f)
    f2d = _run(tmp_path, "2d", 0, eigen=True)
    assert f2d[0] == pytest.approx(ref["flap_hz"], rel=0.03)
    assert f2d[1] == pytest.approx(ref["edge_hz"], rel=0.05)
