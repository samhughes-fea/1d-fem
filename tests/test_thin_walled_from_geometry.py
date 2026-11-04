"""Thin-walled sections cut from an implicit-geometry-kernel blade.

Reference: anisotropic-section-model 0.1 ``ThinWalledSection.from_geometry`` on the same
cut (its example 06 section: NACA 0024, 1 m chord, layout 2D-F with webs at 18 and
48 % chord, glass skins and webs, carbon caps). That model is an ``E t`` / ``G t``
mid-line model, so it places each wall's whole stiffness on the stack's mid-surface.
1d-fem keeps the cap-on-skin stack's membrane-bending coupling ``B``: the carbon sits
inside the glass, nearer the neutral axis, which makes the flapwise ``EI`` 2.1 % lower.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("implicit_geometry_kernel")

from implicit_geometry_kernel import Airfoil, BladeGeometry, Station  # noqa: E402

from pre_processing.parsing.region_laminate_parser import (  # noqa: E402
    format_region_laminate,
    parse_region_laminate_file,
)
from pre_processing.section_integration.thin_walled import Ply, PlyMaterial, WallLaminate  # noqa: E402
from pre_processing.section_integration.thin_walled.from_geometry import thin_walled_from_cut  # noqa: E402

GLASS = PlyMaterial(E1=38.6e9, E2=8.27e9, G12=4.14e9, nu12=0.26, G13=4.14e9, G23=3.063e9, rho=1800.0)
CARBON = PlyMaterial(E1=181e9, E2=10.3e9, G12=7.17e9, nu12=0.28, G13=7.17e9, G23=3.679e9, rho=1600.0)

# anisotropic-section-model, wall_bending="restrained", same cut (z = 3 m)
ASM = {"EA": 1.0151621576e9, "EI_flap": 1.1727352785e7, "EI_edge": 2.1333255973e7, "GJ": 1.9053912317e6,
       "shear_centre_x": 0.0061910508, "centroid_x": 0.0499841074}


def lam(mat, angles, t=0.5e-3):
    return WallLaminate([Ply(mat, a, t) for a in angles])


LAMS = {"skin": lam(GLASS, [0, 45, -45, -45, 45, 0]), "web": lam(GLASS, [45, -45, -45, 45]), "cap": lam(CARBON, [0] * 16)}


@pytest.fixture(scope="module")
def cut():
    af = Airfoil.naca4("0024")
    st = [Station(0.0, 1.0, af, pitch_axis=0.3), Station(6.0, 1.0, af, pitch_axis=0.3)]
    return BladeGeometry(st, "2D-F", 0.003, 0.008, 0.002, web_positions=(0.18, 0.48)).section_at(3.0)


def _numbers(a):
    K = a.classical_stiffness
    yt, zt = a.tension_centre
    EA = K[0, 0]
    return {"EA": EA, "EI_flap": K[1, 1] - EA * zt**2, "EI_edge": K[2, 2] - EA * yt**2,
            "GJ": 1.0 / a.flexibility[5, 5], "shear_centre_x": a.shear_centre[0], "centroid_x": yt}


def test_matches_anisotropic_section_model(cut):
    g = thin_walled_from_cut(cut, LAMS)
    got = _numbers(g.analysis(hoop="restrained"))
    for key in ("EA", "EI_edge", "GJ"):
        assert got[key] == pytest.approx(ASM[key], rel=5e-3), key
    # the cap-on-skin coupling B (see the module docstring)
    assert got["EI_flap"] == pytest.approx(ASM["EI_flap"], rel=0.03)
    assert got["EI_flap"] < ASM["EI_flap"]
    for key in ("shear_centre_x", "centroid_x"):
        assert got[key] == pytest.approx(ASM[key], abs=1e-3), key   # 0.1 % of the chord


def test_walls_stacks_and_labels(cut):
    g = thin_walled_from_cut(cut, LAMS)
    names = set(g.stack_names.values())
    assert "spar_cap_suction_surface+skin_suction_surface:bay_0" in names
    assert {"shear_web_0", "shear_web_1"} <= names
    # a cap stack lists the cap plies first (inside) then the skin plies
    sid = next(k for k, v in g.stack_names.items() if v.startswith("spar_cap_suction"))
    plies = g.laminates[sid].plies
    assert len(plies) == 16 + 6 and plies[0].material is CARBON and plies[-1].material is GLASS
    # three cells: LE, box, TE
    assert g.section.n_cells == 3
    # skin normals point out of the section
    sec = g.section
    for w in sec.walls:
        if g.wall_names[w.wall_id].startswith("skin_suction"):
            assert sec.normal(w)[1] > 0.0


def test_coordinates_from_the_pitch_axis():
    af = Airfoil.naca4("0018")
    st = [Station(0.0, 1.0, af, pitch_axis=0.3, offset=(0.2, -0.5)), Station(2.0, 1.0, af, pitch_axis=0.3, offset=(0.2, -0.5))]
    c = BladeGeometry(st, "2D-F", 0.003, 0.008, 0.002).section_at(1.0)
    a = thin_walled_from_cut(c, LAMS)
    b = thin_walled_from_cut(c, LAMS, origin="blade")
    shift = np.array([0.2, -0.5])
    for pid, p in a.section.points.items():
        np.testing.assert_allclose(b.section.points[pid], p + shift, atol=1e-12)
    np.testing.assert_allclose(_numbers(a.analysis())["EA"], _numbers(b.analysis())["EA"], rtol=1e-12)


def test_missing_laminate_names_the_wall(cut):
    with pytest.raises(KeyError, match="shear_web_0"):
        thin_walled_from_cut(cut, {"skin": LAMS["skin"], "cap": LAMS["cap"]})


def test_thickness_mismatch(cut):
    with pytest.raises(ValueError, match="mm thick in the geometry"):
        thin_walled_from_cut(cut, {**LAMS, "skin": lam(GLASS, [0, 45, -45, 0])})


def test_region_laminate_file(tmp_path, cut):
    from pre_processing.section_integration.thin_walled.section_analysis import ThinWalledSectionAnalysis

    path = tmp_path / "region_laminate.txt"
    path.write_text(format_region_laminate({"skin": 1, "cap": 2, "web": 3, "shear_web_1": 3}))
    mapping = parse_region_laminate_file(str(path), laminate_ids=[1, 2, 3])
    assert mapping["shear_web_1"] == 3
    library = {1: LAMS["skin"], 2: LAMS["cap"], 3: LAMS["web"]}
    g = thin_walled_from_cut(cut, mapping, library)
    a = ThinWalledSectionAnalysis(g.section, g.laminates)
    assert a.hoop == "free"
    assert _numbers(a)["EA"] == pytest.approx(_numbers(thin_walled_from_cut(cut, LAMS).analysis())["EA"], rel=1e-12)
    with pytest.raises(ValueError, match="unknown laminate_id"):
        parse_region_laminate_file(str(path), laminate_ids=[1, 2])


def test_hoop_free_is_the_default():
    from pre_processing.section_integration.thin_walled.section_analysis import wall_law

    np.testing.assert_allclose(wall_law(LAMS["skin"]), wall_law(LAMS["skin"], "free"))
