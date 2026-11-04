"""The shell validation models: mesh, laminate, loads, Abaqus deck and result reduction.

The shell runs themselves need Abaqus (or CalculiX); these tests check that the models
mean what the 1d-fem beams mean, without a solver.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from post_processing.validation_visualisers.abaqus.shell_validation import cases as C
from post_processing.validation_visualisers.abaqus.shell_validation.results import (
    ShellFrames,
    read_abaqus_csv,
    read_calculix_dat,
    section_motion,
)
from post_processing.validation_visualisers.abaqus.shell_validation.shell_model import write_inp


@pytest.fixture(params=["S4R", "S8R"])
def element(request):
    return request.param


def _coords(mesh):
    out = {}
    for k, x in enumerate(mesh.x):
        for j in mesh.section_nodes(k):
            out[mesh.label(k, j)] = np.array([x, *mesh.section_yz[j]])
    return out


@pytest.mark.parametrize("name", list(C.CASES))
def test_shell_normals_are_the_wall_normals(name, element):
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import section_analysis

    case = C.build_case(name, fine=0.25)
    case.model.element_type = element
    mesh = case.model.mesh()
    X = _coords(mesh)
    geo = section_analysis(case.model.section).section
    n_corner = 0
    for wid, *_ in case.model.section.walls:
        n_wall = geo.normal(geo.wall(wid))
        for conn in mesh.elements[wid]:
            a, b, c, d = (X[i] for i in conn[:4])
            n = np.cross(b - a, d - a)
            n /= np.linalg.norm(n)
            np.testing.assert_allclose(n, [0.0, *n_wall], atol=1e-12)
            n_corner += 1
            if element == "S8R":   # midside nodes halfway along the edges
                for m, (p, q) in zip(conn[4:], ((a, b), (b, c), (c, d), (d, a))):
                    np.testing.assert_allclose(X[m], 0.5 * (p + q), atol=1e-12)
    assert n_corner == case.model.n_x * sum(case.model.section.wall_divisions().values())


@pytest.mark.parametrize("make", [C.i_section, C.box_section, C.channel_section])
def test_laminates_match_the_section_analysis(make):
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import section_analysis

    sec = make()
    an = section_analysis(sec)
    for lid, lam in sec.laminates.items():
        np.testing.assert_allclose(lam.membrane_stiffness(), an.laminates[lid].A(), rtol=1e-12, atol=1e-6)


def test_section_weights_add_up_to_the_axial_stiffness(element):
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import section_analysis

    case = C.build_case("ibeam_large_twist", fine=0.5)
    case.model.element_type = element
    mesh = case.model.mesh()
    EA = 1.0 / section_analysis(case.model.section).flexibility[0, 0]
    assert mesh.section_weight.sum() == pytest.approx(EA, rel=2e-3)   # 1d-fem adds the plate shear terms


@pytest.mark.parametrize("load", [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (-3.0, 0.2, -0.7)])
def test_plane_section_forces_have_the_requested_resultant(load, element):
    case = C.build_case("channel_ft", fine=0.5)
    case.model.element_type = element
    mesh = case.model.mesh()
    f = mesh.plane_section_forces(*load)
    y, z = mesh.section_yz.T
    np.testing.assert_allclose([f.sum(), f @ z, -(f @ y)], load, atol=1e-12)


def test_compression_at_the_tension_centre_is_a_uniform_stress():
    from post_processing.validation_visualisers.abaqus.shell_validation.onedfem import tension_centre

    case = C.build_case("channel_ft", fine=0.5)
    mesh = case.model.mesh()
    y_t, z_t = tension_centre(case.model.section)
    f = mesh.plane_section_forces(1.0, z_t, -y_t)
    np.testing.assert_allclose(f / mesh.section_weight, 1.0 / mesh.section_weight.sum(), rtol=1e-6)


def _deck(tmp_path, name, dialect="abaqus", element="S4R"):
    case = C.build_case(name, fine=0.25)
    case.model.element_type = element
    path = tmp_path / f"{name}.inp"
    mesh = write_inp(case.model, str(path), dialect=dialect)
    return case, mesh, path.read_text()


def _cloads(text):
    block = text.split("*Cload\n")[1].split("*")[0]
    return [(int(a), int(b), float(c)) for a, b, c in (ln.split(",") for ln in block.strip().splitlines())]


@pytest.mark.parametrize("name", list(C.CASES))
def test_deck_loads_have_the_beam_resultants(tmp_path, name):
    case, mesh, text = _deck(tmp_path, name)
    X = _coords(mesh)
    n_nodes = mesh.x.size * mesh.n_sec
    refs = {n_nodes + 1: (0.0, *case.model.root.reference), n_nodes + 2: (case.model.length, *case.model.tip.reference)}
    F = np.zeros(3)
    M = np.zeros(3)
    for label, dof, v in _cloads(text):
        p = X.get(label, None)
        if p is None:
            p = np.array(refs[label])
        if dof <= 3:
            f = np.zeros(3)
            f[dof - 1] = v
            F += f
            M += np.cross(p, f)
        else:
            M[dof - 4] += v
    Fb = np.zeros(3)
    Mb = np.zeros(3)
    for x, vals in case.beam["loads"]:
        Fb += vals[:3]
        Mb += np.asarray(vals[3:]) + np.cross([x, 0.0, 0.0], vals[:3])
    scale = max(np.abs(Fb).max(), np.abs(Mb).max() / case.model.length, 1.0)
    np.testing.assert_allclose(F, Fb, atol=1e-9 * scale)
    np.testing.assert_allclose(M, Mb, atol=1e-9 * scale * case.model.length)


def test_deck_cards(tmp_path):
    case, mesh, text = _deck(tmp_path, "ibeam_ltb_path")
    assert text.count("*Coupling") == 2 and text.count("*Kinematic\n2, 3\n") == 2
    assert "*Step, name=LOAD, nlgeom=YES" in text
    assert "*Output, field, time points=LEVELS" in text
    n_el = sum(len(v) for v in mesh.elements.values())
    assert len(re.findall(r"^\d+, \d+, \d+, \d+, \d+$", text, flags=re.M)) == n_el
    # orientations: orthonormal, fibre in the wall plane
    for a in re.findall(r"\*Orientation, name=\S+, system=RECTANGULAR\n(.+)\n", text):
        e = np.array([float(v) for v in a.split(",")]).reshape(2, 3)
        np.testing.assert_allclose(e @ e.T, np.eye(2), atol=1e-9)
    _, _, buckle = _deck(tmp_path, "channel_ft")
    assert "*Step, name=BUCKLE, perturbation\n*Buckle\n12," in buckle
    assert "RP_ROOT, 2, 4" in buckle and "RP_TIP, 2, 4" in buckle
    _, _, box = _deck(tmp_path, "composite_box")
    assert "*Kinematic\n1, 6\n" in box and "ROOT, 1, 6" in box
    assert box.count("*Shell Section, elset=WALL") == 4 and ", 3, CFRP, ORI_W1_P20" in box


def test_calculix_deck_differs_only_in_its_cards(tmp_path):
    _, _, ccx = _deck(tmp_path, "ibeam_large_twist", dialect="calculix", element="S8R")
    assert "*Step, nlgeom, inc=10000" in ccx and "*Node Print, nset=STATIONS, time points=LEVELS" in ccx
    assert "*Shell Section, elset=WALL1, material=STEEL, orientation=" in ccx


# ------------------------------------------------------------------ results
def _rigid_frames(mesh, k, u, v, w, phi):
    c, s = np.cos(phi), np.sin(phi)
    disp = {}
    for j in mesh.section_nodes(k):
        y, z = mesh.section_yz[j]
        disp[mesh.label(k, j)] = np.array([u, c * y - s * z - y + v, s * y + c * z - z + w])
    return disp


@pytest.mark.parametrize("phi", [1e-3, 0.4, 1.2])
def test_section_motion_recovers_a_rigid_motion(phi, element):
    case = C.build_case("channel_ft", fine=0.5)
    case.model.element_type = element
    mesh = case.model.mesh()
    k = mesh.station_index(0.5)
    m = section_motion(mesh, k, _rigid_frames(mesh, k, 0.01, -0.02, 0.03, phi))
    np.testing.assert_allclose(m.as_list(), [0.01, -0.02, 0.03, phi], atol=1e-12)
    assert m.residual < 1e-10


def test_linear_section_motion():
    case = C.build_case("ibeam_large_twist", fine=0.5)
    mesh = case.model.mesh()
    k = mesh.last
    disp = {}
    for j in mesh.section_nodes(k):
        y, z = mesh.section_yz[j]
        disp[mesh.label(k, j)] = np.array([0.0, 0.1 - 0.9 * z, -0.2 + 0.9 * y])
    m = section_motion(mesh, k, disp, linear=True)
    np.testing.assert_allclose(m.as_list(), [0.0, 0.1, -0.2, 0.9], atol=1e-12)


def test_readers(tmp_path):
    csv_path = tmp_path / "r.csv"
    csv_path.write_text("kind,frame,value,label,u1,u2,u3\nbuckle,1,12.5,3,0,1,2\nbuckle,2,20,3,0,3,4\n")
    fr = read_abaqus_csv(str(csv_path))
    assert fr.kind == "buckle" and fr.values == [12.5, 20.0]
    np.testing.assert_allclose(fr.displacements[1][3], [0, 3, 4])
    dat = tmp_path / "r.dat"
    dat.write_text(
        "\n displacements (vx,vy,vz) for set STATIONS and time  0.0000000E+00\n\n"
        "         3  1.0E+00  0.0E+00  0.0E+00\n\n"
        "     B U C K L I N G   F A C T O R   O U T P U T\n\n MODE NO       BUCKLING\n                FACTOR\n\n"
        "      1   0.1115019E+06\n      2   0.1463216E+06\n\n"
        "                    E I G E N V A L U E    N U M B E R     1\n\n\n"
        " displacements (vx,vy,vz) for set STATIONS and time  0.0000000E+00\n\n"
        "         3  1.000000E-01  2.000000E-01  3.000000E-01\n\n"
        "                    E I G E N V A L U E    N U M B E R     2\n\n\n"
        " displacements (vx,vy,vz) for set STATIONS and time  0.0000000E+00\n\n"
        "         3  4.000000E-01  5.000000E-01  6.000000E-01\n")
    fr = read_calculix_dat(str(dat), "buckle")
    assert fr.values == [111501.9, 146321.6]
    np.testing.assert_allclose(fr.displacements[1][3], [0.4, 0.5, 0.6])


def test_path_comparison_of_matching_results():
    from post_processing.validation_visualisers.abaqus.shell_validation.run import compare_path

    case = C.build_case("ibeam_large_twist", fine=0.25)
    mesh = case.model.mesh()
    levels = list(case.beam["levels"])
    beam = {"levels": levels, "stations": {"1": [[0.0, 0.0, 0.0, 0.7 * lv] for lv in levels]}}
    frames = ShellFrames("static", [0.0] + levels,
                         [_rigid_frames(mesh, mesh.last, 0, 0, 0, 0.0)]
                         + [{**_rigid_frames(mesh, mesh.last, 0, 0, 0, 0.7 * lv),
                             **_rigid_frames(mesh, mesh.station_index(0.5), 0, 0, 0, 0.3 * lv)} for lv in levels])
    case.model.stations = (1.0,)
    out = compare_path(case, mesh, frames, beam)
    assert out["passed"] and out["max_error"] < 1e-12
