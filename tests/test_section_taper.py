"""Taper (Resal) coupling of thin-walled sections: ``thin_walled.taper`` and ``section_taper.txt``."""

from __future__ import annotations

import os

import numpy as np
import pytest

from pre_processing.section_integration.section_analysis import ELEMENT_SIGN, SectionAnalysis
from pre_processing.section_integration.thin_walled.laminate import WallLaminate
from pre_processing.section_integration.thin_walled.section_analysis import ThinWalledSectionAnalysis
from pre_processing.section_integration.thin_walled.taper import SectionTaper, section_taper
from pre_processing.section_integration.thin_walled.thin_walled_geometry import ThinWalledSection, Wall

E, NU, T = 70e9, 0.3, 0.004
LAM = {1: WallLaminate.isotropic(E, NU, T)}


def box(b: float, h: float, turn: float = 0.0) -> ThinWalledSection:
    pts = np.array([[-b / 2, -h / 2], [b / 2, -h / 2], [b / 2, h / 2], [-b / 2, h / 2]])
    c, s = np.cos(turn), np.sin(turn)
    pts = pts @ np.array([[c, s], [-s, c]])
    walls = [Wall(i + 1, i + 1, (i + 1) % 4 + 1, 1) for i in range(4)]
    return ThinWalledSection({i + 1: tuple(p) for i, p in enumerate(pts)}, walls, [1])


def membrane_inertia(b, h):
    """int A z^2 ds and int A y^2 ds of the box walls (membrane only)."""
    A = E * T / (1 - NU**2) * (1 - NU**2)   # hoop-free axial membrane stiffness = E t
    Iy = A * (2 * b * (h / 2) ** 2 + 2 * h**3 / 12)
    Iz = A * (2 * h * (b / 2) ** 2 + 2 * b**3 / 12)
    return Iy, Iz


def test_depth_taper_of_a_box():
    b, h0, h1, d = 0.3, 0.2, 0.19, 0.1
    tp = section_taper(box(b, h0), box(b, h1), LAM, d)
    h = 0.5 * (h0 + h1)
    beta = (h1 - h0) / d / h
    Iy, _ = membrane_inertia(b, h)
    C = np.zeros((3, 2))
    C[1, 1] = beta * Iy
    assert np.allclose(tp.coupling, C, atol=1e-9 * Iy)
    assert tp.shear[1, 1] == pytest.approx(beta**2 * Iy, rel=1e-10)
    assert abs(tp.shear[0, 0]) < 1e-9 * Iy and abs(tp.shear[0, 1]) < 1e-9 * Iy


def test_web_shear_follows_resal():
    """gamma_xz = (V_z - M_y h'/h) / GA with the taper coupling added to the section."""
    b, h0, h1, d = 0.3, 0.2, 0.19, 0.1
    h = 0.5 * (h0 + h1)
    beta = (h1 - h0) / d / h
    K = ThinWalledSectionAnalysis(box(b, h), LAM, hoop="free").stiffness
    Kt = K + section_taper(box(b, h0), box(b, h1), LAM, d).matrix()
    # a cantilever under an upward tip force: hogging M_y (top fibres in compression) with V_z > 0,
    # and the depth falling towards the tip, so the inclined flanges take part of V_z
    F = np.array([0.0, -1.0e3, 0.0, 0.0, 1.0e4, 0.0])          # N, M_y, M_z, V_y, V_z, T
    gamma_tapered = np.linalg.solve(Kt, F)[4]
    F_resal = F.copy()
    F_resal[4] = F[4] - beta * F[1]
    gamma_resal = np.linalg.solve(K, F_resal)[4]
    assert gamma_tapered == pytest.approx(gamma_resal, rel=2e-3)
    # the webs carry 5 % less shear than in the prismatic section
    ratio = gamma_tapered / np.linalg.solve(K, F)[4]
    assert ratio == pytest.approx(1.0 - beta * F[1] / F[4], rel=2e-3)
    assert ratio < 0.96


def test_prismatic_pair_has_no_taper():
    tp = section_taper(box(0.3, 0.2), box(0.3, 0.2), LAM, 0.5)
    assert not np.any(tp.matrix())


def test_turning_section_gives_the_helical_fibre_coupling():
    k, d = 0.02, 0.01
    tp = section_taper(box(0.3, 0.2, -0.5 * k * d), box(0.3, 0.2, 0.5 * k * d), LAM, d)
    Iy, Iz = membrane_inertia(0.3, 0.2)
    # r' = k (-z, y): kappa_y couples with gamma_xy by -k Iy, kappa_z with gamma_xz by -k Iz
    assert tp.coupling[1, 0] == pytest.approx(-k * Iy, rel=1e-4)
    assert tp.coupling[2, 1] == pytest.approx(-k * Iz, rel=1e-4)
    assert abs(tp.coupling[0]).max() < 1e-6 * k * Iy


def test_with_taper_in_element_signs():
    sa = SectionAnalysis.from_thin_walled(ThinWalledSectionAnalysis(box(0.3, 0.2), LAM), 1)
    tp = section_taper(box(0.3, 0.2), box(0.3, 0.19), LAM, 0.1)
    out = sa.with_taper(tp)
    assert np.allclose(out.stiffness, sa.stiffness + ELEMENT_SIGN @ tp.matrix() @ ELEMENT_SIGN)
    assert np.allclose(out.classical_stiffness, sa.classical_stiffness)
    assert np.allclose(out.benscoter_stiffness[:6, :6], out.stiffness)
    assert sa.with_taper(SectionTaper.zero()) is sa


def test_section_taper_file_round_trip(tmp_path):
    from pre_processing.parsing.section_taper_parser import format_section_taper_file, parse_section_taper_file

    tp = section_taper(box(0.3, 0.2), box(0.28, 0.19), LAM, 0.1)
    path = tmp_path / "section_taper.txt"
    path.write_text(format_section_taper_file({3: tp}))
    back = parse_section_taper_file(str(path))
    assert set(back) == {3}
    assert np.allclose(back[3].matrix(), tp.matrix(), rtol=1e-15)


def test_attach_applies_the_taper(tmp_path):
    from pre_processing.parsing.laminate_parser import format_laminate_file
    from pre_processing.parsing.section_geometry_parser import format_section_geometry
    from pre_processing.parsing.section_taper_parser import format_section_taper_file
    from pre_processing.section_integration.thin_walled.job_input import attach_section_analyses

    sec = box(0.3, 0.2)
    sec.section_id = 1
    (tmp_path / "laminate.txt").write_text(format_laminate_file(LAM))
    (tmp_path / "section_geometry.txt").write_text(format_section_geometry({1: sec}))
    (tmp_path / "element_section.txt").write_text("[Element_Section]\n[element_id] [section_id]\n0 1\n1 1\n")
    tp = section_taper(box(0.3, 0.2), box(0.3, 0.19), LAM, 0.1)
    (tmp_path / "section_taper.txt").write_text(format_section_taper_file({1: tp}))
    ed = {"ids": np.array([0, 1]), "types": np.array(["LinearTimoshenkoBeamElement3D"] * 2)}
    out = attach_section_analyses(str(tmp_path), ed)["section_analysis_per_element"]
    assert np.allclose(out[1].stiffness - out[0].stiffness, ELEMENT_SIGN @ tp.matrix() @ ELEMENT_SIGN)
    (tmp_path / "section_taper.txt").write_text(format_section_taper_file({5: tp}))
    with pytest.raises(ValueError, match="without a section"):
        attach_section_analyses(str(tmp_path), ed)
