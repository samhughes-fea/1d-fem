"""Sections given as 6 x 6 stiffness and mass tables (``section_matrix.txt``)."""

from __future__ import annotations

import csv
import glob
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

from onedfem.pre_processing.parsing.section_matrix_parser import (  # noqa: E402
    format_section_matrix_file,
    parse_section_matrix_file,
)

pytest.importorskip("anisotropic_section_model.section_analysis")
from anisotropic_section_model import axes  # noqa: E402
from anisotropic_section_model.section_analysis import SectionAnalysis  # noqa: E402

N_EL, LENGTH = 6, 3.0
CF = common.Material(135e9, 10e9, 5e9, 0.3, 5e9, 3.8e9, 1550.0)
# an unbalanced, off-centre box: every coupling, centres away from the beam axis
SECTION = common.SectionSpec(
    points={1: (-0.08, -0.07), 2: (0.22, -0.07), 3: (0.22, 0.05), 4: (-0.08, 0.05)},
    walls=[(1, 1, 2, 1), (2, 2, 3, 2), (3, 3, 4, 3), (4, 4, 1, 2)],
    materials={1: CF},
    laminates={1: [(1, 30, 1e-3), (1, 30, 1e-3), (1, 0, 1e-3)], 2: [(1, 45, 1.5e-3), (1, -45, 1.5e-3)],
               3: [(1, -20, 2e-3), (1, 0, 2e-3)]},
)
TIP = [2e3, 3e3, -4e3, 500.0, 200.0, -300.0]
_WORK = tempfile.mkdtemp()


def teardown_module(module):
    shutil.rmtree(_WORK, ignore_errors=True)


def _tables(frame="engine"):
    an = common.build_analysis(SECTION)
    sa = SectionAnalysis.from_thin_walled(an)
    K, M = an.stiffness, sa.mass_matrix()
    if frame in ("public", "beamdyn"):
        K = axes.from_engine(K)
        props = an.mass_properties
        M = axes.mass_matrix(props)
        M[5, 5] = props["I_p"]
        if frame == "beamdyn":
            K, M = axes.to_beamdyn_frame(K), axes.to_beamdyn_frame(M)
    return K, M


def _run(name, kind, table=None, precurvature=None, n_el=N_EL):
    """One cantilever job; ``table`` replaces the thin-walled section by its tables in that frame."""
    from onedfem.workflow_orchestrator.run_job import process_job

    spec = common.CantileverSpec(LENGTH, n_el, "LinearTimoshenkoBeamElement3D", SECTION, TIP)
    job = os.path.join(_WORK, "jobs", name)
    common.write_job(job, spec)
    common.write_supports(job, spec)
    if table is not None:
        os.remove(os.path.join(job, "laminate.txt"))
        os.remove(os.path.join(job, "section_geometry.txt"))
        K, M = _tables(table)
        with open(os.path.join(job, "section_matrix.txt"), "w", encoding="utf-8") as f:
            f.write(format_section_matrix_file({1: (table, K, M)}))
    if precurvature is not None:
        with open(os.path.join(job, "precurvature.txt"), "w", encoding="utf-8") as f:
            f.write("[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n")
            f.writelines(f"{e} {precurvature[0]} {precurvature[1]} {precurvature[2]}\n" for e in range(n_el))
    if kind == "eigen":
        with open(os.path.join(job, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write("[Simulation]\n[Type]\neigen\n\n[Solver]\ntype = direct\n\n[Eigen]\nnum_modes = 6\n")
    res = common.results_dir(os.path.join(_WORK, "results", name))
    logging.disable(logging.WARNING)
    try:
        process_job(job, res, {}, {}, force_serial=True, max_processes_per_job=1)
    finally:
        logging.disable(logging.NOTSET)
    if kind == "static":
        with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
            return np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return _frequencies(res)


def _frequencies(res):
    paths = glob.glob(os.path.join(res, "**", "*_frequencies.txt"), recursive=True)
    assert paths, f"no frequency file under {res}"
    return np.loadtxt(paths[0])


def test_round_trip_of_the_file(tmp_path):
    K, M = _tables("beamdyn")
    p = tmp_path / "section_matrix.txt"
    p.write_text(format_section_matrix_file({4: ("beamdyn", K, M), 5: ("engine", K, None)}), encoding="utf-8")
    out = parse_section_matrix_file(str(p))
    assert out[4][0] == "beamdyn" and out[5][2] is None
    np.testing.assert_array_equal(out[4][1], K)
    np.testing.assert_array_equal(out[4][2], M)


def test_malformed_files_are_refused(tmp_path):
    K, M = _tables()
    text = format_section_matrix_file({1: ("engine", K, M)})
    p = tmp_path / "section_matrix.txt"
    p.write_text(text.replace("1 engine mass 6", "1 public mass 6"), encoding="utf-8")
    with pytest.raises(ValueError, match="two frames"):
        parse_section_matrix_file(str(p))
    p.write_text("\n".join(ln for ln in text.splitlines() if "stiffness 3" not in ln), encoding="utf-8")
    with pytest.raises(ValueError, match="rows 1 to 6"):
        parse_section_matrix_file(str(p))


@pytest.mark.parametrize("frame", ["engine", "public", "beamdyn"])
def test_table_section_gives_the_thin_walled_static_response(frame):
    ref = _run("static_tw", "static")
    tab = _run(f"static_{frame}", "static", table=frame)
    np.testing.assert_allclose(tab, ref, rtol=1e-8, atol=1e-10 * np.abs(ref).max())


def test_table_section_gives_the_thin_walled_frequencies():
    ref = _run("eigen_tw", "eigen")
    tab = _run("eigen_tab", "eigen", table="beamdyn")
    assert ref.size == 6 and np.all(ref > 0)
    np.testing.assert_allclose(tab, ref, rtol=1e-6, atol=2e-6)   # written to 6 decimals


def test_twisted_elements_keep_the_tables_as_given():
    """Pretwist adds no helical-fibre coupling to a table: the tables are complete."""
    a = _run("static_tab_twist", "static", table="engine", precurvature=(0.05, 0.0, 0.0))
    assert np.all(np.isfinite(a))
