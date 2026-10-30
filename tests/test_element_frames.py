"""Beams off the global x axis, section twist, and curved and pretwisted elements.

References:

- A straight beam turned to any direction and twisted about its axis gives the same
  results as the unturned beam, turned.
- Quarter-circle cantilever of radius ``R`` with an out-of-plane tip force ``P``
  (Roark, table 9.3, case 1): ``delta = P R^3 (pi / (4 EI) + (3 pi / 4 - 2) / GJ)``.
- Helical spring of ``n`` turns, coil radius ``R`` and helix angle ``alpha``, under an
  axial force on its axis (Timoshenko, Strength of Materials II, sec. 54; Wahl 1963):
  ``delta = 2 pi n P R^3 / cos(alpha) (cos^2 alpha / GJ + sin^2 alpha / EI)``.
- Pretwisted strip under tension (Rosen 1980): with ``I_p = int r^2 dA`` and
  ``I_rr = int r^4 dA``, the twist rate is ``-E k I_p eps / (GJ + E k^2 I_rr)`` with
  ``eps`` from ``EA eps + E k I_p phi' = N``.
"""

from __future__ import annotations

import importlib
import os
import shutil
import tempfile

import numpy as np
import pytest

from pre_processing.element_library.dof_layout import mesh_dof_per_node
from pre_processing.element_library.element_factory import ElementFactory
from pre_processing.element_library.element_frame import arc_geometry, rotation_exp, straight_frame
from pre_processing.section_integration.section_analysis import SectionAnalysis
from tests.test_second_order import E, G, NU, channel

EB, TS, BS = "LinearEulerBernoulliBeamElement3D", "LinearTimoshenkoBeamElement3D", "LinearBenscoterBeamElement3D"
RD, LV = "LinearReddyBeamElement3D", "LinearLevinsonBeamElement3D"


class Mesh:
    """Elements of one type on arbitrary nodes, assembled directly."""

    def __init__(self, etype, coords, conn=None, twist=None, curvature=None, sa=None, props=None, dl=None, wall=None):
        coords = np.asarray(coords, dtype=float)
        n_el = coords.shape[0] - 1 if conn is None else len(conn)
        conn = np.array([[i, i + 1] for i in range(n_el)] if conn is None else conn, dtype=np.int64)
        warping = "Benscoter" in etype
        ed = {
            "ids": np.arange(n_el), "connectivity": conn, "types": np.array([etype] * n_el),
            "integration_orders": {k: np.full(n_el, 2) for k in
                                   ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")},
            "warping": np.full(n_el, 1 if warping else 0, dtype=np.int8),
            "twist": np.zeros(n_el) if twist is None else np.broadcast_to(np.asarray(twist, float), (n_el,)).copy(),
        }
        if curvature is not None:
            ed["precurvature_per_element"] = np.broadcast_to(np.asarray(curvature, float), (n_el, 3)).copy()
        if sa is not None:
            ed["section_analysis_per_element"] = tuple([SectionAnalysis.from_thin_walled(sa)] * n_el)
        if wall is not None:      # one generalised wall load per length for every element
            ed["wall_load_per_element"] = tuple([np.asarray(wall, float)] * n_el)
        p = dict(A=1e-3, I_x=2e-7, I_y=1e-7, I_z=1e-7, J_t=2e-7, kappa=0.9, alpha=0.0, y_sc=0.0, z_sc=0.0,
                 Gamma=1e-12)
        p.update(props or {})
        if not warping:     # the warping columns only go with a warping mesh
            for key in ("y_sc", "z_sc", "Gamma"):
                p.pop(key)
        sec = {k: np.broadcast_to(np.asarray(v, float), (n_el,)).copy() for k, v in p.items()}
        mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=G, nu=NU, rho=7850.0).items()}
        Cls = getattr(importlib.import_module(ElementFactory.ELEMENT_CLASS_MAP[etype]), etype)
        tmp = tempfile.mkdtemp()
        jr = os.path.join(tmp, "job")
        for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
            os.makedirs(os.path.join(jr, sub), exist_ok=True)
        self.dpn = mesh_dof_per_node(ed)
        self.n = self.dpn * coords.shape[0]
        self.K = np.zeros((self.n, self.n))
        self.M = np.zeros((self.n, self.n))
        self.elements, self.maps = [], []
        self.F = np.zeros(self.n)
        for e in range(n_el):
            el = Cls(element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
                     material_dictionary=mat, section_dictionary=sec, point_load_array=np.empty((0, 9)),
                     distributed_load_array=np.empty((0, 9)) if dl is None else np.asarray(dl, float),
                     job_results_dir=jr)
            el.logger_operator = None
            g = np.asarray(el.assemble_global_dof_indices())
            self.K[np.ix_(g, g)] += el.element_stiffness_matrix().K_e
            self.M[np.ix_(g, g)] += el.element_mass_matrix().M_e
            if dl is not None or wall is not None:
                self.F[g] += np.real(el.element_force_vector().F_e)
            self.elements.append(el)
            self.maps.append(g)
        shutil.rmtree(tmp, ignore_errors=True)

    def dof(self, node, local):
        return self.dpn * node + local

    def solve(self, F, fixed):
        free = np.setdiff1d(np.arange(self.n), fixed)
        U = np.zeros(self.n)
        U[free] = np.linalg.solve(self.K[np.ix_(free, free)], F[free])
        return U

    def K_sigma(self, U):
        Kg = np.zeros((self.n, self.n))
        for el, g in zip(self.elements, self.maps):
            Kg[np.ix_(g, g)] += el.linear_geometric_stiffness_matrix(U[g])
        return Kg

    def first_buckling(self, F, fixed):
        from scipy import linalg

        U = self.solve(F, fixed)
        free = np.setdiff1d(np.arange(self.n), fixed)
        ev = linalg.eigvals(self.K[np.ix_(free, free)], -self.K_sigma(U)[np.ix_(free, free)])
        ev = ev[np.isfinite(ev) & (np.abs(ev.imag) < 1e-6 * np.abs(ev.real)) & (ev.real > 0)].real
        return float(ev.min())


def _line(p0, d, L, n):
    return np.asarray(p0, float) + np.outer(np.linspace(0.0, L, n + 1), d)


# --------------------------------------------------------------------- frames
def test_straight_frame_defaults():
    np.testing.assert_allclose(straight_frame([0, 0, 0], [2, 0, 0]), np.eye(3), atol=1e-15)
    R = straight_frame([0, 0, 0], [0, 3, 0])                # along +Y: y = -X, z = Z
    np.testing.assert_allclose(R, [[0, -1, 0], [1, 0, 0], [0, 0, 1]], atol=1e-15)
    R = straight_frame([0, 0, 0], [0, 0, 1])                # vertical: y = Y
    np.testing.assert_allclose(R[:, 1], [0, 1, 0], atol=1e-15)
    R = straight_frame([0, 0, 0], [1, 0, 0], 90.0)          # twist turns y to Z
    np.testing.assert_allclose(R[:, 1], [0, 0, 1], atol=1e-15)


def test_arc_geometry_of_a_circle_and_a_helix():
    Rc = 2.0
    # quarter circle in the x-y plane, bending about the local z axis
    arc = arc_geometry([0, 0, 0], [Rc, Rc, 0], [0, 0, 1 / Rc])
    assert arc.length == pytest.approx(np.pi * Rc / 2, rel=1e-12)
    np.testing.assert_allclose(arc.R0[:, 0], [1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(arc.R1[:, 0], [0, 1, 0], atol=1e-12)
    # a helix turn split in eight arcs: each arc's end frame is the next one's start frame
    R, alpha = 0.05, np.deg2rad(10.0)
    pts, k = _helix(R, alpha, 1.0, 8)
    arcs = [arc_geometry(pts[i], pts[i + 1], k) for i in range(8)]
    for a, b in zip(arcs[:-1], arcs[1:]):
        np.testing.assert_allclose(a.R1, b.R0, atol=1e-10)
    np.testing.assert_allclose(arcs[0].R0[:, 0], _helix_tangent(alpha, 0.0), atol=1e-10)


def _helix(R, alpha, turns, n):
    t = np.linspace(0.0, 2 * np.pi * turns, n + 1)
    pts = np.column_stack([R * np.cos(t), R * np.sin(t), R * np.tan(alpha) * t])
    k = np.array([np.sin(alpha) * np.cos(alpha) / R, 0.0, np.cos(alpha) ** 2 / R])   # [torsion, 0, curvature]
    return pts, k


def _helix_tangent(alpha, t):
    return np.array([-np.sin(t) * np.cos(alpha), np.cos(t) * np.cos(alpha), np.sin(alpha)])


# --------------------------------------------------------- straight, rotated
def _turned_cantilever(etype, sa, d, twist, n=8, L=1.5):
    ref = Mesh(etype, _line([0, 0, 0], [1, 0, 0], L, n), sa=sa)
    turned = Mesh(etype, _line([0.3, -0.2, 0.7], d, L, n), twist=twist, sa=sa)
    R = straight_frame([0, 0, 0], d, twist)
    load = np.array([30.0, -50.0, 80.0, 6.0, -4.0, 9.0])
    fixed = list(range(ref.dpn))
    F = np.zeros(ref.n)
    F[ref.dof(n, 0):ref.dof(n, 6)] = load
    U_ref = ref.solve(F, fixed)
    Fg = np.zeros(turned.n)
    Fg[turned.dof(n, 0):turned.dof(n, 3)] = R @ load[:3]
    Fg[turned.dof(n, 3):turned.dof(n, 6)] = R @ load[3:]
    U = turned.solve(Fg, fixed)
    return ref, turned, R, U_ref, U


@pytest.mark.parametrize("etype,sa", [(EB, None), (TS, None), (TS, "channel"), (BS, "channel"), (RD, "channel"),
                                      (LV, "channel")])
@pytest.mark.parametrize("d,twist", [([0.0, 1.0, 0.0], 0.0), ([1.0, 2.0, -2.0], 35.0), ([0.0, 0.0, -1.0], -70.0)])
def test_a_turned_beam_is_the_beam_turned(etype, sa, d, twist):
    sa = channel() if sa else None
    d = np.asarray(d, float) / np.linalg.norm(d)
    ref, turned, R, U_ref, U = _turned_cantilever(etype, sa, d, twist)
    n = len(ref.elements)
    for node in (n // 2, n):
        u_ref = U_ref[ref.dof(node, 0):ref.dof(node, ref.dpn)]
        u = U[turned.dof(node, 0):turned.dof(node, turned.dpn)]
        np.testing.assert_allclose(u[:3], R @ u_ref[:3], atol=1e-10 * np.abs(u_ref).max())
        np.testing.assert_allclose(u[3:6], R @ u_ref[3:6], atol=1e-10 * np.abs(u_ref).max())
        extra_ref = u_ref[6:].copy()
        if etype in (RD, LV):     # a framed mesh holds the slopes as shear strains (framed_element._slope_T)
            extra_ref = extra_ref - np.array([u_ref[5], -u_ref[4]])
        np.testing.assert_allclose(u[6:], extra_ref, atol=1e-10 * np.abs(u_ref).max())   # warping is a scalar
    # mass and geometric stiffness turn the same way
    T = np.kron(np.eye(2 * (n + 1)), R.T) if ref.dpn == 6 else None
    if T is not None:
        np.testing.assert_allclose(turned.M, T.T @ ref.M @ T, atol=1e-10 * np.abs(ref.M).max())
        np.testing.assert_allclose(turned.K_sigma(U), T.T @ ref.K_sigma(U_ref) @ T,
                                   atol=1e-9 * np.abs(ref.K_sigma(U_ref)).max())


def test_buckling_of_a_turned_column():
    sa = channel()
    L, n, P = 2.0, 16, 1.0
    d = np.array([1.0, -1.0, 2.0]) / np.sqrt(6.0)
    out = []
    for coords, R, twist in ((_line([0, 0, 0], [1, 0, 0], L, n), np.eye(3), 0.0),
                             (_line([0, 0, 0], d, L, n), straight_frame([0, 0, 0], d, 25.0), 25.0)):
        m = Mesh(BS, coords, sa=sa, twist=twist)
        F = np.zeros(m.n)
        F[m.dof(n, 0):m.dof(n, 3)] = R @ [-P, 0, 0]
        fixed = [m.dof(0, i) for i in range(7)]
        out.append(m.first_buckling(F, fixed))
    assert out[1] == pytest.approx(out[0], rel=1e-8)


# --------------------------------------------------------- curved cantilever
def _quarter_circle(etype, n, curved, Rc=1.0):
    t = np.linspace(0.0, np.pi / 2, n + 1)
    coords = np.column_stack([Rc * np.sin(t), Rc * (1 - np.cos(t)), 0 * t])   # starts along +x, bends to +y
    b = 0.02
    props = dict(A=b * b, I_y=b**4 / 12, I_z=b**4 / 12, I_x=b**4 / 6, J_t=0.1406 * b**4, kappa=5 / 6)
    m = Mesh(etype, coords, curvature=[0.0, 0.0, 1.0 / Rc] if curved else None, props=props)
    return m, props


@pytest.mark.parametrize("etype", [EB, TS, RD, LV])
@pytest.mark.parametrize("curved,n,tol", [(False, 48, 2e-3), (True, 4, 2e-3)])
def test_quarter_circle_cantilever_out_of_plane(etype, curved, n, tol):
    Rc, P = 1.0, 10.0
    m, p = _quarter_circle(etype, n, curved, Rc)
    F = np.zeros(m.n)
    F[m.dof(n, 2)] = P
    U = m.solve(F, list(range(m.dpn)))
    EI, GJ = E * p["I_y"], G * p["J_t"]
    ref = P * Rc**3 * (np.pi / (4 * EI) + (3 * np.pi / 4 - 2) / GJ)
    assert U[m.dof(n, 2)] == pytest.approx(ref, rel=tol)


@pytest.mark.parametrize("etype", [EB, TS, RD, LV])
def test_quarter_circle_cantilever_in_plane(etype):
    """Radial tip force on a quarter circle: ``delta = P R^3 pi / (4 EI)`` along the force (Roark)."""
    Rc, P, n = 1.0, 10.0, 4
    m, p = _quarter_circle(etype, n, True, Rc)
    F = np.zeros(m.n)
    F[m.dof(n, 0)] = P                      # tip tangent is +y, so x is the radial direction
    U = m.solve(F, list(range(m.dpn)))
    EI = E * p["I_z"]
    # bending only (slender): Castigliano with M = P R cos(theta) gives pi P R^3 / (4 EI)
    assert U[m.dof(n, 0)] == pytest.approx(np.pi * P * Rc**3 / (4 * EI), rel=5e-3)


# --------------------------------------------------------- helical spring
@pytest.mark.parametrize("curved,per_turn", [(False, 48), (True, 8)])
def test_helical_spring(curved, per_turn):
    Rc, alpha, turns, d_w, P = 0.05, np.deg2rad(8.0), 3, 0.005, 20.0
    n = per_turn * turns
    pts, k = _helix(Rc, alpha, turns, n)
    # rigid arms from the coil axis to the wire ends; the load acts on the axis
    axis0 = np.array([0.0, 0.0, pts[0, 2]])
    axis1 = np.array([0.0, 0.0, pts[-1, 2]])
    coords = np.vstack([pts, axis0, axis1])
    a0, a1 = n + 1, n + 2
    conn = [[i, i + 1] for i in range(n)] + [[a0, 0], [a1, n]]
    I = np.pi * d_w**4 / 64
    wire = dict(A=np.pi * d_w**2 / 4, I_y=I, I_z=I, I_x=2 * I, J_t=2 * I, kappa=0.9)
    arm = {kk: v * 1e6 for kk, v in wire.items() if kk != "kappa"}
    props = {kk: np.array([wire[kk]] * n + [arm.get(kk, wire[kk])] * 2) for kk in wire}
    curv = np.vstack([np.tile(k, (n, 1)), np.zeros((2, 3))]) if curved else None
    m = Mesh(TS, coords, conn=conn, curvature=curv, props=props)
    F = np.zeros(m.n)
    F[m.dof(a1, 2)] = P
    fixed = [m.dof(a0, i) for i in range(6)]       # the top end is free: only the axial force acts on it
    U = m.solve(F, fixed)
    EI, GJ = E * I, G * 2 * I
    ref = 2 * np.pi * turns * P * Rc**3 / np.cos(alpha) * (np.cos(alpha) ** 2 / GJ + np.sin(alpha) ** 2 / EI)
    # the closed form leaves out the wire's shear and stretch, about (d / R)^2 / 8 = 0.1 %
    assert U[m.dof(a1, 2)] == pytest.approx(ref, rel=4e-3)


# --------------------------------------------------------- pretwisted strip
@pytest.mark.parametrize("etype", [EB, TS])
def test_pretwisted_strip_untwists_under_tension(etype):
    b, t, L, n, N = 0.05, 0.002, 1.0, 8, 5000.0
    k = np.deg2rad(90.0) / L                  # a quarter turn of pretwist over the length
    A, I_y, I_z = b * t, t**3 * b / 12, b**3 * t / 12
    J = b * t**3 / 3 * (1 - 0.63 * t / b)
    props = dict(A=A, I_y=I_y, I_z=I_z, I_x=I_y + I_z, J_t=J, kappa=5 / 6)
    m = Mesh(etype, _line([0, 0, 0], [1, 0, 0], L, n), curvature=[k, 0.0, 0.0], props=props)
    F = np.zeros(m.n)
    F[m.dof(n, 0)] = N
    U = m.solve(F, list(range(m.dpn)))
    I_p = I_y + I_z
    x, w = np.polynomial.legendre.leggauss(3)            # int r^4 dA over the rectangle
    I_rr = sum(wy * wz * ((0.5 * b * yi) ** 2 + (0.5 * t * zi) ** 2) ** 2 for yi, wy in zip(x, w)
               for zi, wz in zip(x, w)) * A / 4
    K = np.array([[E * A, E * k * I_p], [E * k * I_p, G * J + E * k * k * I_rr]])
    eps, rate = np.linalg.solve(K, [N, 0.0])
    assert rate < 0.0                                       # it untwists
    assert U[m.dof(n, 3)] == pytest.approx(rate * L, rel=1e-3)
    assert U[m.dof(n, 0)] == pytest.approx(eps * L, rel=1e-3)


def test_pretwist_coupling_in_the_section():
    sa = SectionAnalysis.from_isotropic(E=E, G=G, A=1e-4, I_y=2e-8, I_z=8e-9, J_t=1e-9)
    tw = sa.with_pretwist(2.0)
    y, z, wE = sa.axial_points.T
    r2 = y * y + z * z
    assert tw.stiffness[0, 5] == pytest.approx(2.0 * wE @ r2)
    assert tw.stiffness[5, 5] == pytest.approx(sa.stiffness[5, 5] + 4.0 * wE @ (r2 * r2))
    np.testing.assert_allclose(tw.classical_stiffness[3, 0], tw.stiffness[5, 0])
    assert sa.with_pretwist(0.0) is sa


def test_frames_are_refused_where_not_supported():
    with pytest.raises(NotImplementedError, match="linear beams only"):
        Mesh("NonlinearTimoshenkoBeamElement3D", _line([0, 0, 0], [0, 1, 0], 1.0, 2))


def test_rotation_exp_is_a_rotation():
    R = rotation_exp([0.3, -1.2, 0.8])
    np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
    assert np.linalg.det(R) == pytest.approx(1.0)


# --------------------------------------------------------- through the job runner
def _write_quarter_circle_job(job_dir, n, b=0.02, Rc=1.0, P=10.0, twist=None):
    t = np.linspace(0.0, np.pi / 2, n + 1)
    xyz = np.column_stack([Rc * np.sin(t), Rc * (1 - np.cos(t)), 0 * t])
    job_dir.mkdir(parents=True)
    (job_dir / "grid.txt").write_text("[Grid]\n[node_id] [x] [y] [z]\n"
                                      + "".join(f"{i} {x:.15g} {y:.15g} {z:.15g}\n" for i, (x, y, z) in enumerate(xyz)))
    head = ("[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] "
            "[shear_y_order] [shear_z_order] [torsion_order] [load_order]" + (" [twist]" if twist is not None else "") + "\n")
    (job_dir / "element.txt").write_text(head + "".join(
        f"{e} {e} {e + 1} LinearEulerBernoulliBeamElement3D 3 3 3 0 0 3 2" + (f" {twist}" if twist is not None else "")
        + "\n" for e in range(n)))
    (job_dir / "material.txt").write_text("[Material]\n[element_id] [E] [G] [nu] [rho]\n"
                                          + "".join(f"{e} {E} {G} {NU} 7850\n" for e in range(n)))
    I, J = b**4 / 12, 0.1406 * b**4
    (job_dir / "section.txt").write_text("[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
                                         + "".join(f"{e} {b * b} {2 * I} {I} {I} {J} 0.8333 0.0\n" for e in range(n)))
    (job_dir / "precurvature.txt").write_text("[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n"
                                              + "".join(f"{e} 0 0 {1 / Rc}\n" for e in range(n)))
    (job_dir / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(("UX", "UY", "UZ", "RX", "RY", "RZ"))))
    x, y, z = xyz[-1]
    (job_dir / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                                            f"{x:.15g} {y:.15g} {z:.15g} 0 0 {P} 0 0 0\n")
    (job_dir / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n")
    return I, J


def test_curved_cantilever_job(tmp_path):
    import csv

    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    n, Rc, P = 6, 1.0, 10.0
    I, J = _write_quarter_circle_job(tmp_path / "job", n, Rc=Rc, P=P)
    res = setup_job_results_directory("test_curved_cantilever_job")
    process_job(str(tmp_path / "job"), res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="") as fh:
        U = np.array([float(r[1]) for r in list(csv.reader(fh))[1:]]).reshape(n + 1, -1)
    ref = P * Rc**3 * (np.pi / (4 * E * I) + (3 * np.pi / 4 - 2) / (G * J))
    assert U[n, 2] == pytest.approx(ref, rel=2e-3)
    shutil.rmtree(res, ignore_errors=True)


def test_wall_stresses_of_a_turned_beam():
    from pre_processing.section_integration.thin_walled.element_recovery import element_wall_stress

    sa = channel()
    d = np.array([0.3, 0.9, -0.4])
    d /= np.linalg.norm(d)
    ref, turned, R, U_ref, U = _turned_cantilever(BS, sa, d, 20.0)
    wall = sa.section.walls[1]
    for e in (0, 3):
        a = element_wall_stress(ref.elements[e], U_ref[ref.maps[e]], 0.2, sa, wall.wall_id, 0.05)
        b = element_wall_stress(turned.elements[e], U[turned.maps[e]], 0.2, sa, wall.wall_id, 0.05)
        np.testing.assert_allclose(b.resultants, a.resultants, rtol=1e-8, atol=1e-8 * np.abs(a.resultants).max())
