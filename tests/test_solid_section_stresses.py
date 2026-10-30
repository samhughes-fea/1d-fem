"""3D stresses over 2D-meshed sections and their restrained-warping (Benscoter) stiffness."""

from __future__ import annotations

import os

import numpy as np
import pytest

from pre_processing.section_integration.section_analysis import SectionAnalysis
from pre_processing.section_integration.solid.material import SolidMaterial, bond_matrix, material_axes
from pre_processing.section_integration.solid.mesh import SectionMesh
from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis
from pre_processing.section_integration.thin_walled import ThinWalledSection, Wall, WallLaminate
from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial
from pre_processing.section_integration.thin_walled.section_analysis import ThinWalledSectionAnalysis
from tests.test_benscoter_beam_element import BS, fe_cantilever
from tests.test_solid_section_analysis import CFRP, ORTHO, _i_section

E_AL, NU_AL = 70e9, 0.3
G_AL = E_AL / (2 * (1 + NU_AL))
ALU = PlyMaterial(E1=E_AL, E2=E_AL, G12=G_AL, nu12=NU_AL, G13=G_AL, G23=G_AL, rho=2700.0)


# ------------------------------------------------------------------ stresses
def test_rectangle_bending_and_shear_stresses():
    b, h = 0.02, 0.04
    a = SolidSectionAnalysis(SectionMesh.rectangle(b, h, 4, 8, SolidMaterial.isotropic(1e9, 0.0)))
    I = b * h**3 / 12
    for p in a.point_stresses(forces=[0, 5.0, 0, 0, 0, 0]):          # M_y = 5
        assert p.stress[0] == pytest.approx(5.0 * p.z / I, rel=1e-9, abs=1e-6)
        np.testing.assert_allclose(p.stress[1:], 0.0, atol=1e-6)
    # V_z = 3 with nu = 0: the parabola of elementary theory (cubic warping, quadratic elements)
    V = 3.0
    for p in a.point_stresses(forces=[0, 0, 0, 0, V, 0]):
        ref = 1.5 * V / (b * h) * (1 - 4 * p.z**2 / h**2)
        assert p.stress[4] == pytest.approx(ref, abs=1e-2 * 1.5 * V / (b * h))
        assert abs(p.stress[5]) < 1e-3 * V / (b * h)


def test_stresses_carry_the_section_forces():
    sec, lams = _i_section()
    a = SolidSectionAnalysis(SectionMesh.from_thin_walled(sec, lams, elements_along=6))
    F = np.array([1e3, 2e2, -3e2, 4e2, -5e2, 60.0])
    out = np.zeros(6)
    for q, p in zip(a._points, a.point_stresses(forces=F)):
        s, y, z = p.stress, p.y, p.z
        out += q.weight * np.array([s[0], z * s[0], -y * s[0], s[5], s[4], y * s[4] - z * s[5]])
    np.testing.assert_allclose(out, F, rtol=1e-8, atol=1e-8 * np.abs(F).max())


def test_material_and_wall_axes():
    a = SolidSectionAnalysis(SectionMesh.rectangle(0.02, 0.02, 2, 2, ORTHO, theta=30.0, alpha=60.0))
    p = a.point_stresses(forces=[1e3, 10.0, -20.0, 50.0, 30.0, 5.0])[7]
    np.testing.assert_allclose(bond_matrix(material_axes(30.0, 60.0)) @ p.stress_material, p.stress, rtol=1e-12)
    np.testing.assert_allclose(bond_matrix(material_axes(0.0, 60.0)) @ p.stress_wall, p.stress, rtol=1e-12)
    assert p.interlaminar == (p.stress_wall[2], p.stress_wall[3], p.stress_wall[4])
    with pytest.raises(ValueError, match="either"):
        a.point_stresses()


def test_interlaminar_shear_of_a_cross_ply_strip():
    """A [0/90/0/90] strip (wall along y, plies stacked along z) under V_z. With nu = 0 the
    flexure stress is that of a layered beam, tau_xz(z) = V int_z^top E (z' - z_n) dz' / EI,
    continuous across the ply interfaces."""
    mat = SolidMaterial(135e9, 10e9, 10e9, 5e9, 5e9, 3.8e9, 0.0, 0.0, 0.0)
    ply = PlyMaterial(E1=135e9, E2=10e9, G12=5e9, nu12=0.0, G13=5e9, G23=3.8e9)
    t, b = 0.5e-3, 0.02
    angles = (0, 90, 0, 90)
    sec = ThinWalledSection({1: (-b / 2, 0.0), 2: (b / 2, 0.0)}, [Wall(1, 1, 2, 1)])
    lam = {1: WallLaminate([Ply(ply, a, t) for a in angles])}
    mesh = SectionMesh.from_thin_walled(sec, lam, elements_along=4, rows_per_ply=6, materials={ply: mat})
    a = SolidSectionAnalysis(mesh)
    # E_x(z): 0-degree plies have E1 along x, 90-degree plies E2
    bounds = -2 * t + t * np.arange(5)
    Ez = np.array([135e9 if ang == 0 else 10e9 for ang in angles])
    EA = sum(Ez * t)
    z_n = sum(Ez * t * 0.5 * (bounds[:-1] + bounds[1:])) / EA
    EI = sum(Ez * ((bounds[1:] - z_n) ** 3 - (bounds[:-1] - z_n) ** 3) / 3)

    def tau(z):
        k = np.searchsorted(bounds, z, side="right") - 1
        top = sum(Ez[j] * ((bounds[j + 1] - z_n) ** 2 - (max(bounds[j], z) - z_n) ** 2) / 2
                  for j in range(k, 4))
        return V * top / (EI * b)

    V = 10.0
    F = [0, 0, 0, 0, V, 0]
    for k in range(1, 4):                                # the three ply interfaces
        z_i = bounds[k]
        # the row below the interface ends at eta = +1, the row above starts at eta = -1
        below = [e for e in range(len(mesh.elements)) if abs(mesh.nodes[mesh.elements[e][3], 1] - z_i) < 1e-12]
        above = [e for e in range(len(mesh.elements)) if abs(mesh.nodes[mesh.elements[e][0], 1] - z_i) < 1e-12]
        e_lo, e_hi = below[len(below) // 2], above[len(above) // 2]   # mid-width
        t_lo = a.stress_at(e_lo, 0.0, 1.0, forces=F).interlaminar[2]
        t_hi = a.stress_at(e_hi, 0.0, -1.0, forces=F).interlaminar[2]
        # quadratic elements evaluated on their edges: within 2 % with six rows per ply,
        # converging as (rows)^-2
        assert t_lo == pytest.approx(tau(z_i), rel=2e-2)
        assert t_hi == pytest.approx(tau(z_i), rel=2e-2)


def test_core_stresses_of_a_sandwich():
    """Aluminium faces on a soft core: the core carries the shear, the faces the bending."""
    face, core = SolidMaterial.isotropic(70e9, 0.3), SolidMaterial.isotropic(0.1e9, 0.3)
    t_f, t_c, b = 1e-3, 20e-3, 0.05
    parts = [SectionMesh.rectangle(b, t_f, 4, 1, face, centre=(0, -(t_c + t_f) / 2)),
             SectionMesh.rectangle(b, t_c, 4, 4, core, centre=(0, 0)),
             SectionMesh.rectangle(b, t_f, 4, 1, face, centre=(0, (t_c + t_f) / 2))]
    nodes = np.vstack([p.nodes for p in parts])
    elems, off = [], 0
    for p in parts:
        elems += [c + off for c in p.elements]
        off += p.nodes.shape[0]
    from pre_processing.section_integration.solid.mesh import _merge_nodes

    nodes, elems = _merge_nodes(nodes, elems)
    mids = [1] * 4 + [2] * 16 + [1] * 4
    mesh = SectionMesh(nodes, elems, ["Q8"] * len(elems), mids, [0.0] * len(elems), [0.0] * len(elems),
                       {1: face, 2: core})
    a = SolidSectionAnalysis(mesh)
    V = 100.0
    pts = a.point_stresses(forces=[0, 0, 0, 0, V, 0])
    core_tau = [p.stress[4] for p in pts if mesh.material_ids[p.element] == 2]
    # sandwich theory: the core shear is nearly uniform, V / (b d) with d the face-centre distance
    assert np.mean(core_tau) == pytest.approx(V / (b * (t_c + t_f)), rel=0.03)
    assert (max(core_tau) - min(core_tau)) < 0.05 * np.mean(core_tau)
    M = 10.0
    face_sig = [p.stress[0] for p in a.point_stresses(forces=[0, M, 0, 0, 0, 0]) if mesh.material_ids[p.element] == 1]
    assert max(face_sig) == pytest.approx(M / (b * t_f * (t_c + t_f)), rel=0.06)


# ------------------------------------------------------------------ restrained warping
def test_benscoter_stiffness_of_an_isotropic_i_beam():
    sec, _ = _i_section()
    lams = {1: WallLaminate([Ply(ALU, 0.0, 4e-3)]), 2: WallLaminate([Ply(ALU, 0.0, 3e-3)])}
    tw = ThinWalledSectionAnalysis(sec, lams, hoop="free")
    so = SolidSectionAnalysis(SectionMesh.from_thin_walled(sec, lams, elements_along=8, rows_per_ply=2))
    Ks, Kt = so.benscoter_stiffness, tw.benscoter_stiffness
    # Vlasov: E t_f b^3 h^2 / 24
    assert Ks[7, 7] == pytest.approx(E_AL * 4e-3 * 0.15**3 * 0.3**2 / 24, rel=1e-3)
    assert so.warping_constant == Ks[7, 7]
    np.testing.assert_allclose(Ks[:6, :6], so.stiffness)
    # the warping shear d = psi - phi: the 2D section also shears its walls through their
    # thickness (the through-thickness warping follows psi), which Kirchhoff walls cannot
    assert 1.0 < Ks[6, 6] / Kt[6, 6] < 1.5
    assert so.warping_inertia == pytest.approx(tw.warping_inertia, rel=2e-3)
    np.testing.assert_allclose(Ks[:6, 6:], 0.0, atol=1e-9 * np.sqrt(Ks[0, 0] * Ks[7, 7]))


def test_benscoter_stiffness_of_a_composite_i_beam():
    sec, lams = _i_section()
    tw = ThinWalledSectionAnalysis(sec, lams, hoop="free")
    so = SolidSectionAnalysis(SectionMesh.from_thin_walled(sec, lams, elements_along=8))
    Ks, Kt = so.benscoter_stiffness, tw.benscoter_stiffness
    assert Ks[7, 7] == pytest.approx(Kt[7, 7], rel=2e-3)
    # the +-45 flanges couple the warping gradient with flange shear
    assert Ks[3, 7] == pytest.approx(Kt[3, 7], rel=0.1)
    X = so.benscoter_geometric_terms
    assert X.shape == (7, 2)


def test_restrained_torsion_cantilever_with_a_2d_section():
    """Benscoter cantilever, warping fixed at the root, tip torque: 2D and thin-walled sections."""
    sec, _ = _i_section()
    lams = {1: WallLaminate([Ply(ALU, 0.0, 4e-3)]), 2: WallLaminate([Ply(ALU, 0.0, 3e-3)])}
    tw = SectionAnalysis.from_thin_walled(ThinWalledSectionAnalysis(sec, lams, hoop="free"))
    so_an = SolidSectionAnalysis(SectionMesh.from_thin_walled(sec, lams, elements_along=8, rows_per_ply=2))
    so = SectionAnalysis.from_solid(so_an)
    L, T = 3.0, 100.0
    tip_t, psi_t, _ = fe_cantilever(BS, tw, L, [0, 0, 0, T, 0, 0], n_el=16)
    tip_s, psi_s, els = fe_cantilever(BS, so, L, [0, 0, 0, T, 0, 0], n_el=16)
    # GJ (2D vs Kirchhoff) differs by about 1 %; the restraint effect is the same
    assert tip_s[3] == pytest.approx(tip_t[3], rel=0.02)
    free = T * L / so.stiffness[5, 5]
    assert tip_s[3] < 0.9 * free                                    # warping restraint stiffens it
    # a Benscoter element's state: twist rate phi and a warping amplitude c (d = c - phi)
    from pre_processing.section_integration.solid.element_recovery import element_section_state

    el = els[0]
    phi, c = 2e-3, 5e-3
    U_e = np.zeros(14)
    U_e[9] = phi * el.L                                    # theta_x at the second node
    U_e[12] = U_e[13] = c                                  # psi at both nodes
    forces, e8 = element_section_state(el, U_e, 0.0, so_an)
    assert forces is None
    np.testing.assert_allclose(e8, [0, 0, 0, 0, 0, phi, c - phi, 0], atol=1e-12)
    pts = so_an.point_stresses(benscoter=e8)
    # d shears the section: torque of the warping deviation is K[5, 6] d + K[6, 6]-type, no axial stress
    assert max(abs(p.stress[0]) for p in pts) < 1e-6 * max(abs(p.stress[5]) for p in pts)


def test_section_stress_csv(tmp_path):
    from post_processing.solid.section_stresses import SECTION_STRESS_COLUMNS, write_section_stresses

    class _El:
        element_id, x_start, x_end, L = 3, 0.0, 1.0, 1.0

        def assemble_global_dof_indices(self):
            return np.arange(12)

    # an element type without recovery is skipped
    path = write_section_stresses(str(tmp_path), [_El()], [object()], np.zeros(12))
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines == [",".join(SECTION_STRESS_COLUMNS)]
