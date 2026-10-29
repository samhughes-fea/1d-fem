"""Linear Euler-Bernoulli and Timoshenko elements with coupled thin-walled section stiffness.

Cantilevers are checked against the exact solution of the same beam theory: the
section forces of a tip load are statically determinate, the strains follow from the
section flexibility, and the displacements from integrating the kinematics.
"""

from __future__ import annotations

import os
import shutil
import tempfile

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from pre_processing.section_integration.section_analysis import SectionAnalysis
from pre_processing.section_integration.thin_walled import (
    Ply,
    PlyMaterial,
    ThinWalledSection,
    ThinWalledSectionAnalysis,
    Wall,
    WallLaminate,
)

E, NU = 70e9, 0.3
CFRP = PlyMaterial(E1=140e9, E2=10e9, G12=5e9, nu12=0.3, G13=5e9, G23=3.5e9, rho=1600.0)

EB = "LinearEulerBernoulliBeamElement3D"
TS = "LinearTimoshenkoBeamElement3D"


def box_section(b=0.3, h=0.15, offset=(0.0, 0.0), top=None, bottom=None, webs=None):
    """Box with separate top/bottom/web laminates (ids 1, 2, 3); walls run counter-clockwise."""
    y0, z0 = offset
    pts = {1: (y0 - b / 2, z0 - h / 2), 2: (y0 + b / 2, z0 - h / 2), 3: (y0 + b / 2, z0 + h / 2), 4: (y0 - b / 2, z0 + h / 2)}
    walls = [Wall(1, 1, 2, 2), Wall(2, 2, 3, 3), Wall(3, 3, 4, 1), Wall(4, 4, 1, 3)]
    lams = {1: top, 2: bottom, 3: webs}
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), lams)


def cas_box(**kw):
    """Circumferentially asymmetric stiffness: +20 on top, -20 on the bottom, gives bending-twist coupling."""
    top = WallLaminate([Ply(CFRP, 20.0, 0.5e-3)] * 4)
    bottom = WallLaminate([Ply(CFRP, -20.0, 0.5e-3)] * 4)
    webs = WallLaminate([Ply(CFRP, a, 0.5e-3) for a in (45, -45, -45, 45)])
    return box_section(top=top, bottom=bottom, webs=webs, **kw)


def cus_box(**kw):
    """Circumferentially uniform stiffness (+20 everywhere along s): extension-twist coupling."""
    lam = WallLaminate([Ply(CFRP, 20.0, 0.5e-3)] * 4)
    return box_section(top=lam, bottom=lam, webs=lam, **kw)


# --------------------------------------------------------------------------- #
# Exact cantilever solution (right-handed), clamped at x = 0, tip load at x = L
# --------------------------------------------------------------------------- #


def exact_tip(flex_rh, L, P, theory):
    """Tip displacement [u, v, w, theta_x, theta_y, theta_z] for tip force/moment P = [Fx, Fy, Fz, Mx, My, Mz].

    Section forces (right-handed) at x: N = Fx, V_y = Fy, V_z = Fz, T = Mx,
    M_y = My + Fz (x - L), M_z = Mz - Fy (x - L).
    """
    Fx, Fy, Fz, Mx, My, Mz = P
    C = np.array(flex_rh, dtype=float)
    if theory == "eb":
        idx = [0, 1, 2, 5]
        C4 = np.linalg.inv(np.linalg.inv(C[np.ix_(idx, idx)]))
        C = np.zeros((6, 6))
        C[np.ix_(idx, idx)] = C4

    def forces(x):
        return np.array([Fx, My + Fz * (x - L), Mz - Fy * (x - L), Fy, Fz, Mx])

    # integrate theta' = kappa, then v' = gamma_xy + theta_z, w' = gamma_xz - theta_y
    xg, wg = leggauss(12)
    def strains(x):
        return C @ forces(x)

    def theta(x):
        s = 0.5 * x * (xg + 1)
        k = np.array([strains(si) for si in s])
        return 0.5 * x * np.array([wg @ k[:, 5], wg @ k[:, 1], wg @ k[:, 2]])  # theta_x, theta_y, theta_z

    s = 0.5 * L * (xg + 1)
    ps = np.array([strains(si) for si in s])
    th = np.array([theta(si) for si in s])
    u = 0.5 * L * wg @ ps[:, 0]
    v = 0.5 * L * wg @ (ps[:, 3] + th[:, 2])
    w = 0.5 * L * wg @ (ps[:, 4] - th[:, 1])
    tx, ty, tz = theta(L)
    return np.array([u, v, w, tx, ty, tz])


# --------------------------------------------------------------------------- #
# Finite element cantilever with the 1d-fem elements
# --------------------------------------------------------------------------- #


def fe_tip(etype, section, L, P, n_el=8, orders=None):
    tmp = tempfile.mkdtemp()
    try:
        jr = os.path.join(tmp, "job")
        for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
            os.makedirs(os.path.join(jr, sub), exist_ok=True)
        x = np.linspace(0.0, L, n_el + 1)
        coords = np.column_stack([x, np.zeros_like(x), np.zeros_like(x)])
        o = orders or {"axial": 2, "bending_y": 2, "bending_z": 2, "shear_y": 2, "shear_z": 2, "torsion": 2, "load": 2}
        ed = {
            "ids": np.arange(n_el),
            "connectivity": np.column_stack([np.arange(n_el), np.arange(1, n_el + 1)]),
            "types": np.array([etype] * n_el),
            "integration_orders": {k: np.full(n_el, v) for k, v in o.items()},
            "section_analysis_per_element": tuple([section] * n_el),
        }
        mat = {k: np.full(n_el, v) for k, v in dict(E=E, G=E / 2.6, nu=NU, rho=0.0).items()}
        sec = {k: np.full(n_el, v) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9, kappa=5 / 6, alpha=0.0).items()}
        from pre_processing.element_library.element_factory import ElementFactory  # noqa: F401  (registry import)
        if etype == EB:
            from pre_processing.element_library.linear.beam.zero_order_shear_deformation_theory.euler_bernoulli.linear_euler_bernoulli_3D import (
                LinearEulerBernoulliBeamElement3D as Cls,
            )
        else:
            from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
                LinearTimoshenkoBeamElement3D as Cls,
            )
        n_dof = 6 * (n_el + 1)
        K = np.zeros((n_dof, n_dof))
        elements = []
        for e in range(n_el):
            el = Cls(
                element_id=e, element_dictionary=ed, grid_dictionary={"coordinates": coords},
                material_dictionary=mat, section_dictionary=sec,
                point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)),
                job_results_dir=jr,
            )
            el.logger_operator = None
            Ke = el.element_stiffness_matrix().K_e
            K[6 * e:6 * e + 12, 6 * e:6 * e + 12] += Ke
            elements.append({"Ke": Ke, "Me": el.element_mass_matrix().M_e if section.mass else None, "D": el._D6()})
        F = np.zeros(n_dof)
        F[-6:] = P
        free = np.arange(6, n_dof)
        U = np.zeros(n_dof)
        U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
        return U[-6:], elements, K
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def rel_err(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


LOADS = {
    "axial": [1e4, 0, 0, 0, 0, 0],
    "torque": [0, 0, 0, 2e2, 0, 0],
    "moment_y": [0, 0, 0, 0, 3e2, 0],
    "moment_z": [0, 0, 0, 0, 0, 3e2],
    "force_z": [0, 0, 2e2, 0, 0, 0],
    "force_y": [0, 1e2, 0, 0, 0, 0],
}


@pytest.mark.parametrize("load", list(LOADS))
@pytest.mark.parametrize("make", [cas_box, cus_box])
def test_euler_bernoulli_cantilever_matches_exact_coupled_solution(load, make):
    tw = make(offset=(0.02, -0.01))
    sa = SectionAnalysis.from_thin_walled(tw)
    L = 2.0
    exact = exact_tip(tw.flexibility, L, np.array(LOADS[load], float), "eb")
    if load.startswith("force"):
        # bending moment varies along the beam; with bending-twist coupling the twist rate
        # varies too, which linear twist interpolation captures at O(h^2)
        e8 = rel_err(fe_tip(EB, sa, L, LOADS[load], n_el=8)[0], exact)
        e32 = rel_err(fe_tip(EB, sa, L, LOADS[load], n_el=32)[0], exact)
        assert e32 < 2e-4
        assert e32 < e8 / 10 or e32 < 1e-9
    else:
        fe, _, _ = fe_tip(EB, sa, L, LOADS[load], n_el=4)
        assert rel_err(fe, exact) < 1e-9


@pytest.mark.parametrize("load", list(LOADS))
@pytest.mark.parametrize("make", [cas_box, cus_box])
def test_timoshenko_cantilever_converges_to_exact_coupled_solution(load, make):
    tw = make(offset=(0.02, -0.01))
    sa = SectionAnalysis.from_thin_walled(tw)
    L = 2.0
    exact = exact_tip(tw.flexibility, L, np.array(LOADS[load], float), "ts")
    e1 = rel_err(fe_tip(TS, sa, L, LOADS[load], n_el=8)[0], exact)
    e2 = rel_err(fe_tip(TS, sa, L, LOADS[load], n_el=32)[0], exact)
    assert e2 < 2e-3
    assert e2 <= e1 + 1e-12


def test_coupling_is_active():
    """The CAS box couples bending and twist: a tip moment about y twists the beam."""
    tw = cas_box()
    sa = SectionAnalysis.from_thin_walled(tw)
    tip, _, _ = fe_tip(EB, sa, 2.0, LOADS["moment_y"], n_el=2)
    assert abs(tip[3]) > 0.05 * abs(tip[4])
    tip, _, _ = fe_tip(EB, SectionAnalysis.from_thin_walled(cus_box()), 2.0, LOADS["axial"], n_el=2)
    assert abs(tip[3]) > 0.0


def test_timoshenko_no_shear_locking_for_slender_coupled_beam():
    tw = cas_box()
    sa = SectionAnalysis.from_thin_walled(tw)
    L = 40.0
    exact = exact_tip(tw.flexibility, L, np.array(LOADS["force_z"], float), "ts")
    assert rel_err(fe_tip(TS, sa, L, LOADS["force_z"], n_el=16)[0], exact) < 5e-3


def test_timoshenko_stiffness_rank_and_symmetry():
    sa = SectionAnalysis.from_thin_walled(cas_box(offset=(0.05, 0.03)))
    _, els, _ = fe_tip(TS, sa, 1.0, LOADS["axial"], n_el=1)
    Ke = els[0]["Ke"]
    np.testing.assert_allclose(Ke, Ke.T, rtol=1e-12, atol=1e-6)
    ev = np.linalg.eigvalsh(Ke)
    assert np.sum(ev < 1e-9 * ev.max()) == 6


@pytest.mark.parametrize("etype", [EB, TS])
def test_section_mass_matrix_rigid_body_kinetic_energy(etype):
    lam = WallLaminate([Ply(CFRP, 0.0, 1e-3), Ply(PlyMaterial.isotropic(E, NU, 0.0), 0.0, 1e-3)])
    tw = box_section(top=lam, bottom=lam, webs=lam, offset=(0.04, -0.03))
    sa = SectionAnalysis.from_thin_walled(tw)
    L = 1.5
    _, els, _ = fe_tip(etype, sa, L, LOADS["axial"], n_el=1)
    Me = els[0]["Me"]
    mp = sa.mass
    m, ym, zm = mp["m"], mp["y_m"], mp["z_m"]
    # rigid translation along z: T = m L / 2
    U = np.zeros(12); U[[2, 8]] = 1.0
    assert U @ Me @ U == pytest.approx(m * L, rel=1e-10)
    # rigid rotation about x: u_y = -z, u_z = y  ->  2T = I_p L
    U = np.zeros(12); U[[3, 9]] = 1.0
    assert U @ Me @ U == pytest.approx(mp["I_p"] * L, rel=1e-10)
    # translation y then rotation x: coupling -m z_m L
    Uy = np.zeros(12); Uy[[1, 7]] = 1.0
    Ux = np.zeros(12); Ux[[3, 9]] = 1.0
    assert Uy @ Me @ Ux == pytest.approx(-m * zm * L, rel=1e-10)
    # rigid rotation about z (theta_z = 1, u_y = x): 2T = int (I_zz + m x^2 - 2 x S...) -> check u_x part
    Uz = np.zeros(12); Uz[[5, 11]] = 1.0; Uz[7] = L
    Uxa = np.zeros(12); Uxa[[0, 6]] = 1.0
    assert Uxa @ Me @ Uz == pytest.approx(-m * ym * L, rel=1e-10)


@pytest.mark.parametrize("etype", [EB, TS])
def test_isotropic_thin_walled_box_matches_section_txt_element(etype):
    """A thin-walled isotropic box gives the same element stiffness as section.txt with its EA, EI, GJ."""
    lam = WallLaminate.isotropic(E, NU, 2e-3)
    tw = box_section(top=lam, bottom=lam, webs=lam)
    sa = SectionAnalysis.from_thin_walled(tw)
    _, els, _ = fe_tip(etype, sa, 1.0, LOADS["axial"], n_el=1)
    D = sa.euler_bernoulli_D() if etype == EB else sa.stiffness
    assert not sa.is_coupled
    np.testing.assert_allclose(els[0]["D"], D)
    # same element stiffness as the isotropic section.txt path with these EA, EI, GA, GJ
    from pre_processing.section_integration.section_analysis import SectionAnalysis as SA
    K = sa.stiffness
    iso = SA(np.diag(np.diag(K)), np.diag(np.diag(sa.classical_stiffness)))
    _, els_iso, _ = fe_tip(etype, iso, 1.0, LOADS["axial"], n_el=1)
    np.testing.assert_allclose(els[0]["Ke"], els_iso[0]["Ke"], rtol=1e-9, atol=1e-9 * np.abs(els[0]["Ke"]).max())


def test_warping_dof_with_thin_walled_section_is_refused():
    sa = SectionAnalysis.from_thin_walled(cas_box())
    from pre_processing.element_library.linear.beam.first_order_shear_deformation_theory.timoshenko.linear_timoshenko_3D import (
        LinearTimoshenkoBeamElement3D,
    )

    tmp = tempfile.mkdtemp()
    try:
        jr = os.path.join(tmp, "job")
        for sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
            os.makedirs(os.path.join(jr, sub), exist_ok=True)
        ed = {
            "ids": np.array([0]), "connectivity": np.array([[0, 1]]), "types": np.array([TS]),
            "warping": np.array([1], dtype=np.int8),
            "integration_orders": {k: np.array([2]) for k in ("axial", "bending_y", "bending_z", "shear_y", "shear_z", "torsion", "load")},
            "section_analysis_per_element": (sa,),
        }
        sec = {k: np.array([v]) for k, v in dict(A=1e-3, I_x=0.0, I_y=1e-8, I_z=1e-8, J_t=1e-9, kappa=5 / 6, alpha=0.0,
                                                y_sc=0.0, z_sc=0.0, Gamma=1e-9).items()}
        with pytest.raises(NotImplementedError, match="warping"):
            LinearTimoshenkoBeamElement3D(
                element_id=0, element_dictionary=ed, grid_dictionary={"coordinates": np.array([[0.0, 0, 0], [1.0, 0, 0]])},
                material_dictionary={k: np.array([v]) for k, v in dict(E=E, G=E / 2.6, nu=NU, rho=0.0).items()},
                section_dictionary=sec, point_load_array=np.empty((0, 9)), distributed_load_array=np.empty((0, 9)),
                job_results_dir=jr,
            )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_process_job_with_thin_walled_section(tmp_path):
    """End to end: element_section.txt drives a coupled cantilever through process_job."""
    import csv

    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    L, n_el = 2.0, 4
    job = tmp_path / "job_thin_walled"
    job.mkdir()
    nodes = "\n".join(f"{i}  {L * i / n_el:.6f}  0.0  0.0" for i in range(n_el + 1))
    (job / "grid.txt").write_text(f"[Grid]\n[node_id] [x] [y] [z]\n{nodes}\n", encoding="utf-8")
    rows = "\n".join(f"{e} {e} {e + 1} {EB} 3 3 3 0 0 3 2" for e in range(n_el))
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
        f"[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order]\n{rows}\n",
        encoding="utf-8",
    )
    (job / "material.txt").write_text(
        "[Material]\n[element_id] [E] [G] [nu] [rho]\n" + "\n".join(f"{e} 7e10 2.7e10 0.3 0" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    (job / "section.txt").write_text(
        "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
        + "\n".join(f"{e} 1e-3 0.0 1e-8 1e-8 1e-9 0.8333 0.0" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "\n".join(f"{i} 0 {d} 0.0 displacement" for i, d in enumerate(["UX", "UY", "UZ", "RX", "RY", "RZ"])) + "\n",
        encoding="utf-8",
    )
    (job / "point_load.txt").write_text(
        "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
        f"{L} 0.0 0.0 0.0 0.0 200.0 0.0 0.0 0.0\n",
        encoding="utf-8",
    )
    (job / "simulation_settings.txt").write_text(
        "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n",
        encoding="utf-8",
    )
    (job / "laminate.txt").write_text(
        "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]\n"
        "1 140e9 10e9 5e9 0.3 5e9 3.5e9 1600\n"
        "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n"
        + "".join(f"1 {k} 1 20.0 0.5e-3\n" for k in range(1, 5))
        + "".join(f"2 {k} 1 -20.0 0.5e-3\n" for k in range(1, 5))
        + "".join(f"3 {k} 1 {a} 0.5e-3\n" for k, a in enumerate((45, -45, -45, 45), start=1)),
        encoding="utf-8",
    )
    (job / "section_geometry.txt").write_text(
        "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
        "1 1 -0.15 -0.075\n1 2 0.15 -0.075\n1 3 0.15 0.075\n1 4 -0.15 0.075\n"
        "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
        "1 1 1 2 2\n1 2 2 3 3\n1 3 3 4 1\n1 4 4 1 3\n",
        encoding="utf-8",
    )
    (job / "element_section.txt").write_text(
        "[Element_Section]\n[element_id] [section_id]\n" + "\n".join(f"{e} 1" for e in range(n_el)) + "\n",
        encoding="utf-8",
    )
    res = setup_job_results_directory("test_thin_walled_section_job")
    process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    U = np.array([float(r[1]) for r in rows[1:]]) if len(rows) > 2 else np.array([float(x) for x in rows[1][1:]])
    tip = U[-6:]
    expected, _, _ = fe_tip(EB, SectionAnalysis.from_thin_walled(cas_box()), L, [0, 0, 200.0, 0, 0, 0], n_el=n_el)
    assert rel_err(tip, expected) < 1e-9
    assert abs(tip[3]) > 0.05 * abs(tip[4])
