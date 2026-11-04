"""Geometric stiffness of the GBT element: local, distortional and global buckling."""

from __future__ import annotations

import importlib
import os
import shutil
import sys
import tempfile

import numpy as np
import pytest
import scipy.linalg as sla
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "examples", "thin_walled", "gbt_buckling"))

from test_gbt_cross_section_modes import E, G, NU, _ed, box, iso  # noqa: E402

from finite_strip import column_stress, critical_stress  # noqa: E402
from pre_processing.element_library.dof_layout import mesh_dof_per_node  # noqa: E402
from pre_processing.section_integration.section_analysis import SectionAnalysis  # noqa: E402
from pre_processing.section_integration.thin_walled import ThinWalledSection, ThinWalledSectionAnalysis, Wall  # noqa: E402

GBT = "LinearGBTBeamElement3D"
_WORK = tempfile.mkdtemp()
TS = "LinearTimoshenkoBeamElement3D"


def _elements(etype, sa, L, n_el, modes):
    from pre_processing.element_library.element_factory import ElementFactory

    Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[etype]), etype)
    ed = _ed(etype, n_el, modes, sa)
    ed["section_analysis_per_element"] = tuple([sa] * n_el)
    x = np.linspace(0.0, L, n_el + 1)
    jr = os.path.join(_WORK, "job")
    for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
        os.makedirs(os.path.join(jr, sub), exist_ok=True)
    mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=0.0).items()}
    sec = {k: np.full(n_el, v) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9, kappa=5 / 6,
                                                alpha=0.0, y_sc=0.0, z_sc=0.0, Gamma=1e-9).items()}
    els = []
    for e in range(n_el):
        el = Cls(element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": np.column_stack([x, 0 * x, 0 * x])},
                 material_dictionary=mat, section_dictionary=sec, point_load_array=np.empty((0, 9)),
                 distributed_load_array=np.empty((0, 9)), job_results_dir=jr)
        el.logger_operator = None
        els.append(el)
    return els, mesh_dof_per_node(ed)


def _pinned_column(etype, sa, L, n_el, modes, n_lowest=3):
    """Lowest buckling loads of a column on forks and diaphragms, compressed by a unit end force."""
    els, dpn = _elements(etype, sa, L, n_el, modes)
    p = sa.gbt_sv_modes if etype == GBT else 0
    n = dpn * (n_el + 1)
    K = np.zeros((n, n))
    dofs = []
    for el in els:
        g = el.assemble_global_dof_indices()
        dofs.append(g)
        K[np.ix_(g, g)] += el.element_stiffness_matrix().K_e
    diaphragm = [6 + p + 2 * k for k in range(modes)] if etype == GBT else []
    fixed = {0, 1, 2, 3, *diaphragm} | {dpn * n_el + j for j in [1, 2, 3, *diaphragm]}
    free = np.array(sorted(set(range(n)) - fixed))
    F = np.zeros(n)
    F[dpn * n_el] = -1.0
    U = np.zeros(n)
    U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    Kg = np.zeros((n, n))
    for el, g in zip(els, dofs):
        Kg[np.ix_(g, g)] += el.linear_geometric_stiffness_matrix(U[g])
    mu = sla.eigh(-Kg[np.ix_(free, free)], K[np.ix_(free, free)], eigvals_only=True)
    lam = np.sort(1.0 / mu[mu > 1e-12 * mu.max()])
    return lam[:n_lowest], K[np.ix_(free, free)], Kg[np.ix_(free, free)], els, U, dofs


def teardown_module(module):
    shutil.rmtree(_WORK, ignore_errors=True)


def test_euler_load_of_a_tube_without_modes_is_the_timoshenko_one():
    tw = box(0.3, 0.3, 2e-3)
    lam_t = _pinned_column(TS, SectionAnalysis.from_thin_walled(tw), 6.0, 12, 0)[0]
    lam_g = _pinned_column(GBT, SectionAnalysis.from_thin_walled(tw, gbt=True), 6.0, 12, 0)[0]
    np.testing.assert_allclose(lam_g[:2], lam_t[:2], rtol=2e-5)


def test_rigid_section_terms_are_those_of_the_beam():
    """Under compression the GBT stresses give the beam's N and Wagner term ``W``."""
    tw = box(0.3, 0.2, 2e-3)
    _, _, _, els, U, dofs = _pinned_column(GBT, SectionAnalysis.from_thin_walled(tw, gbt_modes=2), 2.0, 4, 2)
    c = els[1].section_geometric_terms(U[dofs[1]])
    _, _, _, els_t, U_t, dofs_t = _pinned_column(TS, SectionAnalysis.from_thin_walled(tw), 2.0, 4, 0)
    ct = els_t[1].section_geometric_terms(U_t[dofs_t[1]])
    assert c[0] == pytest.approx(-1.0, rel=1e-9)
    assert c[3] == pytest.approx(ct[3], rel=1e-6)                       # polar term of N
    np.testing.assert_allclose(c[[1, 2, 4, 5, 6]], 0.0, atol=1e-9)


def test_square_tube_local_buckling_is_the_plate_value():
    """Four simply supported plates, four half-waves: k = 4."""
    b, t, M = 0.3, 2e-3, 8
    sa = SectionAnalysis.from_thin_walled(box(b, b, t, n_sub=4), gbt_modes=M)
    lam = _pinned_column(GBT, sa, 4 * b, 16, M)[0]
    plate = 4 * np.pi**2 * E / (12 * (1 - NU**2)) * (t / b) ** 2 * 4 * b * t
    assert lam[0] == pytest.approx(plate, rel=1e-3)


def test_finite_strip_reference_is_the_plate_value():
    b, t = 0.3, 2e-3
    pts = [(-b / 2, -b / 2), (b / 2, -b / 2), (b / 2, b / 2), (-b / 2, b / 2)]
    s = critical_stress([(pts[i], pts[(i + 1) % 4]) for i in range(4)], t, E, NU, b, n_sub=8)
    assert s == pytest.approx(4 * np.pi**2 * E / (12 * (1 - NU**2)) * (t / b) ** 2, rel=2e-4)


def _lipped_channel(n_sub):
    h, bf, c, t = 0.1, 0.06, 0.01, 2e-3
    y_c = (2 * bf * bf / 2 + 2 * c * bf) / (h + 2 * bf + 2 * c)
    P = {1: (bf, h / 2 - c), 2: (bf, h / 2), 3: (0.0, h / 2), 4: (0.0, -h / 2), 5: (bf, -h / 2), 6: (bf, -h / 2 + c)}
    P = {k: (y - y_c, z) for k, (y, z) in P.items()}
    W = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]
    tw = ThinWalledSectionAnalysis(ThinWalledSection(P, [Wall(i + 1, a, b, 1) for i, (a, b) in enumerate(W)]),
                                   {1: iso(t)}, n_sub=n_sub)
    return tw, [(P[a], P[b]) for a, b in W], t * (h + 2 * bf + 2 * c)


@pytest.mark.parametrize("L, n_el, kind", [(0.08, 8, "local"), (0.3, 30, "distortional")])
def test_lipped_channel_matches_the_finite_strip_column(L, n_el, kind):
    M = 10
    tw, walls, A = _lipped_channel(2)
    sa = SectionAnalysis.from_thin_walled(tw, gbt_modes=M)
    lam = _pinned_column(GBT, sa, L, n_el, M)[0]
    s_fsm, _ = column_stress(walls, 2e-3, E, NU, L)
    # GBT's Vlasov modes are slightly stiffer than full plate theory
    assert 1.0 <= lam[0] / A / s_fsm < 1.03


def test_sparse_buckling_solver_matches_dense():
    from processing.buckling.linear_buckling import _sparse_buckling

    tw, _, _ = _lipped_channel(2)
    sa = SectionAnalysis.from_thin_walled(tw, gbt_modes=6)
    lam, K, Kg, *_ = _pinned_column(GBT, sa, 0.3, 24, 6, n_lowest=4)
    out = _sparse_buckling(csr_matrix(K), csr_matrix(Kg), 4)
    assert out is not None
    np.testing.assert_allclose(out[0], lam, rtol=1e-8)
    A = csr_matrix(K)
    A[0, 1] += 1.0                                       # a nonsymmetric K goes to the dense solver
    assert _sparse_buckling(A, csr_matrix(Kg), 4) is None


@pytest.mark.slow
def test_2d_section_geometric_terms_match_thin_walled():
    """The same box as a 2D mesh: Euler load as the Timoshenko beam's, and the geometric
    term of the first local mode as the thin-walled section's."""
    from test_gbt_2d_sections import box_mesh

    from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis

    an = SolidSectionAnalysis(box_mesh())
    lam_g = _pinned_column(GBT, SectionAnalysis.from_solid(an, gbt=True), 6.0, 12, 0)[0]
    lam_t = _pinned_column(TS, SectionAnalysis.from_solid(an), 6.0, 12, 0)[0]
    np.testing.assert_allclose(lam_g[:2], lam_t[:2], rtol=1e-3)    # shear warping vs shear factor

    def first_mode_terms(sa, m):
        Kd, _ = sa.gbt_matrices(m)
        TX, _ = sa.gbt_geometric_terms(m)
        E0 = np.zeros(Kd.shape[0])
        E0[0] = 1.0
        X = TX @ E0
        oa = 6 + 2 * sa.gbt_sv_modes
        kk = np.pi / 0.3                                     # a 0.3 m half-wave: the local mode is lowest

        def single(j):
            C, D, B = (Kd[oa + i * m + j, oa + i * m + j] for i in range(3))
            CB = Kd[oa + j, oa + 2 * m + j]
            return (C * kk**4 + D * kk**2 + B - 2 * CB * kk**2) / (X[3 + j, 3 + j] * kk**2)

        k = int(np.argmin([single(j) for j in range(m)]))
        return X[3 + k, 3 + k] / X[0, 0]

    x2 = first_mode_terms(SectionAnalysis.from_solid(an, gbt_modes=4), 4)
    xt = first_mode_terms(SectionAnalysis.from_thin_walled(box(0.3, 0.12, 2e-3, n_sub=4), gbt_modes=8), 8)
    assert x2 == pytest.approx(xt, rel=0.01)
