"""Third-order shear (Reddy, Levinson) elements on isotropic and thin-walled sections."""

from __future__ import annotations

import importlib
import os
import shutil
import tempfile

import numpy as np
import pytest

from pre_processing.element_library.dof_layout import mesh_dof_per_node, row_extra_dofs
from pre_processing.element_library.element_factory import ElementFactory
from pre_processing.section_integration.section_analysis import SectionAnalysis
from pre_processing.section_integration.thin_walled import (
    ThinWalledSection,
    ThinWalledSectionAnalysis,
    Wall,
    WallLaminate,
)
from pre_processing.section_integration.thin_walled.laminate import Ply, PlyMaterial

REDDY, LEVINSON, TIMO = "LinearReddyBeamElement3D", "LinearLevinsonBeamElement3D", "LinearTimoshenkoBeamElement3D"
E, NU = 2.0e11, 0.3
G = E / (2 * (1 + NU))
B_SEC, H_SEC = 0.1, 0.2  # solid rectangle: width (y), depth (z)


def _ed(etype, n_el):
    o = {"axial": 2, "bending_y": 2, "bending_z": 2, "shear_y": 1, "shear_z": 1, "torsion": 2, "load": 2}
    return {
        "ids": np.arange(n_el),
        "connectivity": np.column_stack([np.arange(n_el), np.arange(1, n_el + 1)]),
        "types": np.array([etype] * n_el),
        "warping": np.zeros(n_el, dtype=np.int8),
        "integration_orders": {k: np.full(n_el, v) for k, v in o.items()},
    }


def solve(etype, L, n_el, tip=None, q=None, sa=None, clamp_slopes=True, return_K=False):
    """Cantilever along x by direct assembly, clamped at x = 0.

    ``tip`` maps local DOF -> force at the tip node; ``q`` is a uniform distributed load
    ``[q_x, q_y, q_z, m_x, m_y, m_z]``. ``sa`` is a thin-walled SectionAnalysis (else the
    solid rectangle). Returns nodal U (n_nodes, stride).
    """
    Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[etype]), etype)
    ed = _ed(etype, n_el)
    if sa is not None:
        ed["section_analysis_per_element"] = tuple([sa] * n_el)
    x = np.linspace(0.0, L, n_el + 1)
    coords = np.column_stack([x, 0 * x, 0 * x])
    A, Iy, Iz = B_SEC * H_SEC, B_SEC * H_SEC**3 / 12, H_SEC * B_SEC**3 / 12
    mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()}
    sec = {k: np.full(n_el, v) for k, v in dict(A=A, I_x=Iy + Iz, I_y=Iy, I_z=Iz, J_t=0.2 * Iy, kappa=5 / 6,
                                                alpha=0.0).items()}
    dist = np.empty((0, 9))
    if q is not None:
        xs = np.linspace(0.0, L, 7)
        dist = np.column_stack([xs, 0 * xs, 0 * xs] + [np.full_like(xs, v) for v in q])
    jr = os.path.join(tempfile.mkdtemp(), "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    try:
        dpn = mesh_dof_per_node(ed)
        n = dpn * (n_el + 1)
        K = np.zeros((n, n))
        F = np.zeros(n)
        for e in range(n_el):
            el = Cls(element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
                     material_dictionary=mat, section_dictionary=sec, point_load_array=np.empty((0, 9)),
                     distributed_load_array=dist, job_results_dir=jr)
            el.logger_operator = None
            g = el.assemble_global_dof_indices()
            K[np.ix_(g, g)] += el.element_stiffness_matrix().K_e
            if q is not None:
                F[g] += el.element_force_vector().F_e
        for local, val in (tip or {}).items():
            F[dpn * n_el + local] += val
        fixed = set(range(6 if not clamp_slopes else dpn))
        free = np.array(sorted(set(range(n)) - fixed))
        U = np.zeros(n)
        U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        U = U.reshape(n_el + 1, dpn)
        return (U, K) if return_K else U
    finally:
        shutil.rmtree(os.path.dirname(jr), ignore_errors=True)


def _plane(I, A=B_SEC * H_SEC, kappa=5 / 6):
    EI, kGA = E * I, kappa * G * A
    lam = np.sqrt(kGA / (EI / 84.0))
    return EI, kGA, lam


# --------------------------------------------------------------------- #
# Section terms
# --------------------------------------------------------------------- #
def _iso_mat():
    return PlyMaterial(E, E, G, NU, G, G, 0.0)


def i_beam(h=0.2, tw=0.01, bf=0.002, tf=0.0005, n_sub=4):
    m = _iso_mat()
    pts = {1: (-bf / 2, -h / 2), 2: (0.0, -h / 2), 3: (bf / 2, -h / 2), 4: (-bf / 2, h / 2), 5: (0.0, h / 2),
           6: (bf / 2, h / 2)}
    walls = [Wall(1, 1, 2, 2), Wall(2, 2, 3, 2), Wall(3, 4, 5, 2), Wall(4, 5, 6, 2), Wall(5, 2, 5, 1)]
    lams = {1: WallLaminate([Ply(m, 0.0, tw)]), 2: WallLaminate([Ply(m, 0.0, tf)])}
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), lams, n_sub=n_sub)


def composite_box(angles_cap=(0, 0, 30, 30), angles_web=(45, -45), b=0.4, h=0.2):
    cf = PlyMaterial(135e9, 10e9, 5e9, 0.3, 5e9, 3.8e9, 1550.0)
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 2), Wall(3, 3, 4, 1), Wall(4, 4, 1, 2)]
    lams = {1: WallLaminate([Ply(cf, a, 1e-3) for a in angles_cap]),
            2: WallLaminate([Ply(cf, a, 1e-3) for a in angles_web])}
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), lams, n_sub=2)


def test_rectangle_web_gives_reddy_section_term():
    """The slaved shear warping of a thin web is Reddy's cubic: chi_z stiffness E I / 84."""
    an = i_beam()
    K8 = an.third_order_stiffness
    I_web = 0.01 * 0.2**3 / 12
    assert K8[7, 7] == pytest.approx(E * I_web / 84, rel=5e-3)
    np.testing.assert_allclose(K8[:6, :6], an.stiffness, rtol=1e-12)
    # the slaved warping carries no N, M_y, M_z
    assert np.abs(K8[6:, :3]).max() < 1e-9 * K8[7, 7]
    np.testing.assert_allclose(an.third_order_rate, 0.0, atol=1e-12)


def test_isotropic_section_analysis_third_order_terms():
    A, Iy, Iz = 0.02, 6e-5, 1.5e-5
    sa = SectionAnalysis.from_isotropic(E=E, G=G, A=A, I_y=Iy, I_z=Iz, J_t=1e-5)
    K8 = sa.third_order_stiffness
    np.testing.assert_allclose(K8[:6, :6], sa.stiffness)
    assert K8[6, 6] == pytest.approx(E * Iz / 84) and K8[7, 7] == pytest.approx(E * Iy / 84)
    np.testing.assert_allclose(sa.third_order_rate, 0.0)


def test_thin_walled_section_analysis_is_element_signed():
    an = composite_box()
    sa = SectionAnalysis.from_thin_walled(an)
    S8 = np.diag([1.0, -1, 1, 1, 1, 1, 1, 1])
    np.testing.assert_allclose(sa.third_order_stiffness, S8 @ an.third_order_stiffness @ S8)
    np.testing.assert_allclose(sa.third_order_rate, an.third_order_rate @ S8[:6, :6])


def test_reddy_and_levinson_rows_carry_two_slope_dofs():
    ed = _ed(REDDY, 2)
    assert row_extra_dofs(ed, 0) == 2 and mesh_dof_per_node(ed) == 8
    ed["types"] = np.array([LEVINSON, TIMO])
    assert row_extra_dofs(ed, 1) == 0 and mesh_dof_per_node(ed) == 8


# --------------------------------------------------------------------- #
# Isotropic cantilevers against closed forms
# --------------------------------------------------------------------- #
@pytest.mark.parametrize("etype", [REDDY, LEVINSON])
def test_free_slopes_tip_load_is_timoshenko(etype):
    """Shear warping free at the root, constant shear: the Timoshenko beam, node for node."""
    L, P = 0.6, 1.0e4
    for n in (2, 8):
        U = solve(etype, L, n, tip={2: P, 1: 0.5 * P}, clamp_slopes=False)
        T = solve(TIMO, L, n, tip={2: P, 1: 0.5 * P})
        np.testing.assert_allclose(U[:, :6], T[:, :6], rtol=1e-10, atol=1e-16)


@pytest.mark.parametrize("plane", ["xz", "xy"])
def test_reddy_clamped_tip_load_matches_closed_form(plane):
    """gamma(0) = 0 at the clamp: w = P L^3 / 3EI + (P / kGA)(L - tanh(lam L) / lam), lam^2 = 84 kGA / EI."""
    L, P = 0.4, 1.0e4
    I = B_SEC * H_SEC**3 / 12 if plane == "xz" else H_SEC * B_SEC**3 / 12
    EI, kGA, lam = _plane(I)
    exact = P * L**3 / (3 * EI) + P / kGA * (L - np.tanh(lam * L) / lam)
    dof = 2 if plane == "xz" else 1
    errs = [abs(solve(REDDY, L, n, tip={dof: P})[-1, dof] / exact - 1) for n in (32, 64)]
    assert errs[1] < 2e-4 and errs[1] < errs[0]
    timo = P * L**3 / (3 * EI) + P * L / kGA
    assert exact < timo  # the restrained shear warping stiffens the beam


def test_uniform_load_free_warping_matches_timoshenko_tip():
    """With free shear warping the boundary layers leave the tip at q L^4 / 8EI + q L^2 / 2kGA."""
    L, qz = 0.5, 2.0e5
    EI, kGA, _ = _plane(B_SEC * H_SEC**3 / 12)
    exact = qz * L**4 / (8 * EI) + qz * L**2 / (2 * kGA)
    for et in (REDDY, LEVINSON):
        w = solve(et, L, 64, q=[0, 0, qz, 0, 0, 0], clamp_slopes=False)[-1, 2]
        assert w == pytest.approx(exact, rel=3e-4)


def test_uniform_load_clamped_reddy_matches_closed_form():
    """gamma(0) = 0, gamma'(L) = 0 with kGA gamma - H gamma'' = q (L - x)."""
    L, qz = 0.5, 2.0e5
    EI, kGA, lam = _plane(B_SEC * H_SEC**3 / 12)
    # gamma = q (L - x) / kGA + a cosh(lam x) + b sinh(lam x), a = -q L / kGA, b from gamma'(L) = 0
    a, c = -qz * L / kGA, qz / (kGA * lam)
    int_gamma = qz * L**2 / (2 * kGA) + (a * np.tanh(lam * L) + c * (1 - 1 / np.cosh(lam * L))) / lam
    exact = qz * L**4 / (8 * EI) + int_gamma
    w = solve(REDDY, L, 96, q=[0, 0, qz, 0, 0, 0])[-1, 2]
    assert w == pytest.approx(exact, rel=3e-4)


def test_thin_beam_limit_is_euler_bernoulli():
    """No shear locking: L/h = 100 converges to P L^3 / 3EI."""
    L, P = 20.0, 1.0
    EI = E * B_SEC * H_SEC**3 / 12
    for et in (REDDY, LEVINSON):
        w = solve(et, L, 16, tip={2: P})[-1, 2]
        assert w == pytest.approx(P * L**3 / (3 * EI), rel=2e-3)


# --------------------------------------------------------------------- #
# Thin-walled sections
# --------------------------------------------------------------------- #
def test_thin_walled_box_free_slopes_is_saint_venant_and_clamp_stiffens():
    an = composite_box(angles_cap=(0, 0, 0, 0), angles_web=(45, -45))
    sa = SectionAnalysis.from_thin_walled(an)
    L, P = 2.0, 1.0e4
    T = solve(TIMO, L, 8, tip={2: P}, sa=sa)
    for et in (REDDY, LEVINSON):
        U = solve(et, L, 8, tip={2: P}, sa=sa, clamp_slopes=False)
        np.testing.assert_allclose(U[:, :6], T[:, :6], rtol=1e-9, atol=1e-15)
    Uc = solve(REDDY, L, 8, tip={2: P}, sa=sa)
    assert 0.95 * T[-1, 2] < Uc[-1, 2] < T[-1, 2]


def test_coupled_section_levinson_is_unsymmetric_reddy_symmetric():
    """Unbalanced caps couple the slaved warping to the classical strains."""
    an = composite_box(angles_cap=(30, 30, 30, 30), angles_web=(45, 45))
    sa = SectionAnalysis.from_thin_walled(an)
    K8 = sa.third_order_stiffness
    assert np.abs(K8[6:, 3:6]).max() > 1e-6 * np.abs(K8[6:, 6:]).max()
    _, Kr = solve(REDDY, 1.0, 2, tip={2: 1.0}, sa=sa, return_K=True)
    _, Kl = solve(LEVINSON, 1.0, 2, tip={2: 1.0}, sa=sa, return_K=True)
    np.testing.assert_allclose(Kr, Kr.T, atol=1e-6 * np.abs(Kr).max())
    assert np.abs(Kl - Kl.T).max() > 1e-9 * np.abs(Kl).max()


def test_wall_stress_recovery_reddy_matches_saint_venant():
    from pre_processing.section_integration.thin_walled.element_recovery import element_wall_stress
    from pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

    an = composite_box(angles_cap=(0, 0, 0, 0), angles_web=(45, -45))
    sa = SectionAnalysis.from_thin_walled(an)
    L, P, n_el = 2.0, 1.0e4, 4
    U = solve(REDDY, L, n_el, tip={2: P}, sa=sa, clamp_slopes=False)
    Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[REDDY]), REDDY)
    ed = _ed(REDDY, n_el)
    ed["section_analysis_per_element"] = tuple([sa] * n_el)
    x = np.linspace(0.0, L, n_el + 1)
    jr = tempfile.mkdtemp()
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    try:
        mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()}
        sec = {k: np.full(n_el, v) for k, v in dict(A=1e-3, I_x=2e-8, I_y=1e-8, I_z=1e-8, J_t=1e-9).items()}
        el = Cls(element_id=1, element_dictionary=ed, grid_dictionary={"coordinates": np.column_stack([x, 0 * x, 0 * x])},
                 material_dictionary=mat, section_dictionary=sec, point_load_array=np.empty((0, 9)),
                 distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
        Ue = np.r_[U[1, :6], U[2, :6], U[1, 6:], U[2, 6:]]
        xm = 0.5 * (x[1] + x[2])
        got = element_wall_stress(el, Ue, 0.0, an, 3, 0.2)
        ref = SaintVenantRecovery(an).at(3, 0.2, [0, P * (xm - L), 0, 0, P, 0])
        np.testing.assert_allclose(got.resultants, ref.resultants, rtol=1e-6, atol=1e-6 * np.abs(ref.resultants).max())
    finally:
        shutil.rmtree(jr, ignore_errors=True)


def test_warping_column_is_refused():
    ed = _ed(REDDY, 1)
    ed["warping"] = np.array([1], dtype=np.int8)
    Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[REDDY]), REDDY)
    jr = tempfile.mkdtemp()
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    try:
        with pytest.raises(NotImplementedError, match="Benscoter"):
            Cls(element_id=0, element_dictionary=ed, grid_dictionary={"coordinates": np.array([[0, 0, 0], [1.0, 0, 0]])},
                material_dictionary={k: np.array([v]) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()},
                section_dictionary={k: np.array([v]) for k, v in dict(A=1e-3, I_x=2e-8, I_y=1e-8, I_z=1e-8, J_t=1e-9,
                                                                       kappa=5 / 6, alpha=0.0, y_sc=0.0, z_sc=0.0,
                                                                       Gamma=1e-9).items()},
                point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
    finally:
        shutil.rmtree(jr, ignore_errors=True)


# --------------------------------------------------------------------- #
# Trend against 2D plane-stress elasticity
# --------------------------------------------------------------------- #
def _plane_stress_q8_tip(L, h, P, nx, ny, t=B_SEC):
    """Tip deflection (neutral axis) of a fully clamped plane-stress cantilever, serendipity Q8,
    parabolic end shear ``P``."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spl

    D = E / (1 - NU**2) * np.array([[1, NU, 0], [NU, 1, 0], [0, 0, (1 - NU) / 2]])
    X, Y = np.linspace(0, L, 2 * nx + 1), np.linspace(-h / 2, h / 2, 2 * ny + 1)
    idx = -np.ones((2 * nx + 1, 2 * ny + 1), int)
    k = 0
    for i in range(2 * nx + 1):
        for j in range(2 * ny + 1):
            if not (i % 2 and j % 2):
                idx[i, j] = k
                k += 1
    coords = np.zeros((k, 2))
    for i in range(2 * nx + 1):
        for j in range(2 * ny + 1):
            if idx[i, j] >= 0:
                coords[idx[i, j]] = X[i], Y[j]
    gp, gw = np.polynomial.legendre.leggauss(3)
    xi_n = np.array([-1, 1, 1, -1, 0, 1, 0, -1.0])
    et_n = np.array([-1, -1, 1, 1, -1, 0, 1, 0.0])
    rows, cols, vals = [], [], []
    for ex in range(nx):
        for ey in range(ny):
            i, j = 2 * ex, 2 * ey
            nodes = [idx[i, j], idx[i + 2, j], idx[i + 2, j + 2], idx[i, j + 2],
                     idx[i + 1, j], idx[i + 2, j + 1], idx[i + 1, j + 2], idx[i, j + 1]]
            xy = coords[nodes]
            Ke = np.zeros((16, 16))
            for a, wa in zip(gp, gw):
                for b, wb in zip(gp, gw):
                    dN = np.zeros((8, 2))
                    for n in range(4):
                        dN[n] = (0.25 * xi_n[n] * (1 + et_n[n] * b) * (2 * xi_n[n] * a + et_n[n] * b),
                                 0.25 * et_n[n] * (1 + xi_n[n] * a) * (xi_n[n] * a + 2 * et_n[n] * b))
                    for n in (4, 6):
                        dN[n] = (-a * (1 + et_n[n] * b), 0.5 * (1 - a * a) * et_n[n])
                    for n in (5, 7):
                        dN[n] = (0.5 * xi_n[n] * (1 - b * b), -b * (1 + xi_n[n] * a))
                    J = dN.T @ xy
                    dNx = dN @ np.linalg.inv(J).T
                    Bm = np.zeros((3, 16))
                    Bm[0, 0::2], Bm[1, 1::2] = dNx[:, 0], dNx[:, 1]
                    Bm[2, 0::2], Bm[2, 1::2] = dNx[:, 1], dNx[:, 0]
                    Ke += Bm.T @ D @ Bm * np.linalg.det(J) * wa * wb * t
            dofs = np.ravel([[2 * n, 2 * n + 1] for n in nodes])
            rows += list(np.repeat(dofs, 16))
            cols += list(np.tile(dofs, 16))
            vals += list(Ke.ravel())
    K = sp.csr_matrix((vals, (rows, cols)), shape=(2 * k, 2 * k))
    F = np.zeros(2 * k)
    for ey in range(ny):
        j = 2 * ey
        nodes = np.array([idx[2 * nx, j], idx[2 * nx, j + 1], idx[2 * nx, j + 2]])
        y0, y1 = Y[j], Y[j + 2]
        for s, ws in zip(gp, gw):
            N = np.array([0.5 * s * (s - 1), 1 - s * s, 0.5 * s * (s + 1)])
            y = 0.5 * (y0 + y1) + 0.5 * (y1 - y0) * s
            F[2 * nodes + 1] += N * 1.5 * P / h * (1 - 4 * y * y / h**2) * 0.5 * (y1 - y0) * ws
    fixed = np.ravel([[2 * idx[0, j], 2 * idx[0, j] + 1] for j in range(2 * ny + 1)])
    free = np.setdiff1d(np.arange(2 * k), fixed)
    U = np.zeros(2 * k)
    U[free] = spl.spsolve(K[free][:, free].tocsc(), F[free])
    return U[2 * idx[2 * nx, ny] + 1]


def test_thick_clamped_cantilever_reddy_is_closer_to_plane_stress_than_timoshenko():
    """L/h = 2: the fully clamped 2D solid restrains warping; the clamped Reddy beam captures part of it."""
    L, P = 0.4, 1.0e4
    w2d = _plane_stress_q8_tip(L, H_SEC, P, 20, 10)
    w_t = solve(TIMO, L, 64, tip={2: P})[-1, 2]
    w_r = solve(REDDY, L, 64, tip={2: P})[-1, 2]
    assert w2d < w_r < w_t
    assert abs(w_r / w2d - 1) < 0.01 and abs(w_t / w2d - 1) > abs(w_r / w2d - 1)
