"""Nodal and interior point loads, and line loads, on every linear beam element.

A cantilever clamped at x = 0 with an end load ``[F_x, F_y, F_z, M_x, M_y, M_z]``
(global, right-handed, like the nodal DOFs) has the closed-form tip response

    u = F_x L / EA,                 theta_x = M_x L / GJ,
    v = F_y L^3 / 3EI_z + M_z L^2 / 2EI_z,   theta_z = F_y L^2 / 2EI_z + M_z L / EI_z,
    w = F_z L^3 / 3EI_y - M_y L^2 / 2EI_y,   theta_y = -F_z L^2 / 2EI_y + M_y L / EI_y,

plus ``F L / kGA`` shear terms for shear-deformable elements (made negligible here).
Before the fix Euler-Bernoulli elements dropped nodal ``M_y`` and ``M_z`` and the
Timoshenko family applied ``M_y`` with the wrong sign.
"""

from __future__ import annotations

import csv
import io
import contextlib
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from pre_processing.element_library.beam_load_interpolation import (
    beam_field_interpolation,
    distributed_load_vector,
    point_load_vector,
)

E, G, NU = 2.1e11, 8.1e10, 0.3
A, IY, IZ, J = 1e-3, 2e-7, 3e-7, 1e-7
DOFS = ("UX", "UY", "UZ", "RX", "RY", "RZ")

LINEAR = (
    "LinearEulerBernoulliBeamElement3D",
    "LinearTimoshenkoBeamElement3D",
    "LinearReddyBeamElement3D",
    "LinearLevinsonBeamElement3D",
)


def _write_job(job, etype, coords, loads, dist=None, nonlinear=False, iz=IZ):
    job.mkdir(parents=True)
    n = len(coords) - 1
    (job / "grid.txt").write_text("[Grid]\n[node_id] [x] [y] [z]\n"
                                  + "".join(f"{i} {float(x)!r} {float(y)!r} {float(z)!r}\n" for i, (x, y, z) in enumerate(coords)))
    (job / "element.txt").write_text(
        "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] [bending_z_order] "
        "[shear_y_order] [shear_z_order] [torsion_order] [load_order]\n"
        + "".join(f"{e} {e} {e + 1} {etype} 3 3 3 2 2 3 3\n" for e in range(n)))
    (job / "material.txt").write_text("[Material]\n[element_id] [E] [G] [nu] [rho]\n"
                                      + "".join(f"{e} {E} {G} {NU} 0\n" for e in range(n)))
    (job / "section.txt").write_text("[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
                                     + "".join(f"{e} {A} {IY + iz} {IY} {iz} {J} 1e6 0\n" for e in range(n)))
    (job / "precurvature.txt").write_text("[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n"
                                          + "".join(f"{e} 0 0 0\n" for e in range(n)))
    (job / "prescribed_displacement.txt").write_text(
        "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
        + "".join(f"{k} 0 {d} 0.0 displacement\n" for k, d in enumerate(DOFS)))
    (job / "point_load.txt").write_text("[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
                                        + "".join(" ".join(repr(float(v)) for v in row) + "\n" for row in loads))
    if dist is not None:
        (job / "distributed_load.txt").write_text(
            "[Distributed load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
            + "".join(" ".join(repr(float(v)) for v in row) + "\n" for row in dist))
    if nonlinear:
        settings = ("[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-9\nrelative_tolerance = 1e-10\n"
                    "max_iterations = 40\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n\n"
                    "[Nonlinear]\nload_steps = 1\n")
    else:
        settings = "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n"
    (job / "simulation_settings.txt").write_text(settings)


def _run(tmp_path, etype, coords, loads, dist=None, nonlinear=False, iz=IZ):
    from workflow_orchestrator.run_job import process_job, setup_job_results_directory

    job = tmp_path / f"job_{len(list(tmp_path.iterdir()))}"
    _write_job(job, etype, coords, loads, dist, nonlinear, iz)
    res = setup_job_results_directory("test_point_moments")
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            process_job(str(job), res, {}, {}, force_serial=True, max_processes_per_job=1)
        with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="") as fh:
            U = np.array([float(r[1]) for r in list(csv.reader(fh))[1:]])
    finally:
        shutil.rmtree(res, ignore_errors=True)
    return U.reshape(len(coords), -1)[:, :6]


def _exact_tip(P, L):
    Fx, Fy, Fz, Mx, My, Mz = P
    return np.array([
        Fx * L / (E * A),
        Fy * L**3 / (3 * E * IZ) + Mz * L**2 / (2 * E * IZ),
        Fz * L**3 / (3 * E * IY) - My * L**2 / (2 * E * IY),
        Mx * L / (G * J),
        -Fz * L**2 / (2 * E * IY) + My * L / (E * IY),
        Fy * L**2 / (2 * E * IZ) + Mz * L / (E * IZ),
    ])


def _line(n, L=2.0, direction=(1.0, 0.0, 0.0)):
    d = np.asarray(direction, float) / np.linalg.norm(direction)
    return np.outer(np.linspace(0.0, L, n + 1), d)


# ----------------------------------------------------------------- interpolation
@pytest.mark.parametrize("bending", ["hermite", "linear"])
def test_field_interpolation_reproduces_nodal_values(bending):
    rng = np.random.default_rng(3)
    U = rng.normal(size=12)
    np.testing.assert_allclose(beam_field_interpolation(-1.0, 1.7, bending) @ U, U[:6], atol=1e-14)
    np.testing.assert_allclose(beam_field_interpolation(1.0, 1.7, bending) @ U, U[6:], atol=1e-14)


def test_hermite_rotations_are_slopes():
    """theta_z = v' and theta_y = -w' everywhere along a Hermite element."""
    L, h = 1.3, 1e-6
    U = np.random.default_rng(5).normal(size=12)
    for xi in (-0.6, 0.1, 0.8):
        N0, Np, Nm = (beam_field_interpolation(x, L, "hermite") @ U for x in (xi, xi + h, xi - h))
        slope = (Np - Nm) / (2 * h) * (2 / L)
        assert N0[5] == pytest.approx(slope[1], rel=1e-7)
        assert N0[4] == pytest.approx(-slope[2], rel=1e-7)


def test_consistent_line_load_end_moments():
    """A uniform q_y over a Hermite element gives q L / 2 forces and -/+ q L^2 / 12 end moments."""
    L, q = 2.0, 3.0
    xg, wg = np.polynomial.legendre.leggauss(4)
    F = distributed_load_vector(xg, wg, np.tile([0, q, 0, 0, 0, 0], (4, 1)), L, "hermite")
    np.testing.assert_allclose(F[[1, 7]], [q * L / 2] * 2)
    np.testing.assert_allclose(F[[5, 11]], [q * L**2 / 12, -q * L**2 / 12])
    G_ = distributed_load_vector(xg, wg, np.tile([0, 0, q, 0, 0, 0], (4, 1)), L, "hermite")
    np.testing.assert_allclose(G_[[4, 10]], [-q * L**2 / 12, q * L**2 / 12])


def test_point_moment_at_a_node_goes_to_its_rotation():
    for bending in ("hermite", "linear"):
        F = point_load_vector(1.0, 2.0, [0, 0, 0, 1.0, 2.0, 3.0], bending)
        np.testing.assert_allclose(F, np.r_[np.zeros(9), 1.0, 2.0, 3.0], atol=1e-15)


# ----------------------------------------------------------------- through the job runner
LOADS = {"Fx": [1e3, 0, 0, 0, 0, 0], "Fy": [0, 10, 0, 0, 0, 0], "Fz": [0, 0, 10, 0, 0, 0],
         "Mx": [0, 0, 0, 10, 0, 0], "My": [0, 0, 0, 0, 10, 0], "Mz": [0, 0, 0, 0, 0, 10]}


@pytest.mark.parametrize("etype", LINEAR)
def test_end_loads_every_direction(tmp_path, etype):
    L, n = 2.0, 16
    P = np.array(list(LOADS.values()), float)
    loads = [np.r_[L, 0, 0, row] for row in P]
    for row in loads:
        U = _run(tmp_path, etype, _line(n, L), [row])
        ref = _exact_tip(row[3:], L)
        # Timoshenko-type linear elements with reduced shear integration: 1 / (4 n^2) in bending
        np.testing.assert_allclose(U[-1], ref, rtol=0, atol=2e-3 * np.abs(ref).max())


def test_euler_bernoulli_end_moments_are_exact(tmp_path):
    L = 2.0
    for row in ([L, 0, 0, 0, 0, 0, 0, 7.0, 0], [L, 0, 0, 0, 0, 0, 0, 0, -4.0]):
        U = _run(tmp_path, LINEAR[0], _line(3, L), [row])
        np.testing.assert_allclose(U[-1], _exact_tip(row[3:], L), rtol=1e-9, atol=1e-15)


def test_nonlinear_euler_bernoulli_end_moment(tmp_path):
    L = 2.0
    row = [L, 0, 0, 0, 0, 0, 0, 0.5, -0.5]
    U = _run(tmp_path, "NonlinearEulerBernoulliBeamElement3D", _line(4, L), [row], nonlinear=True)
    ref = _exact_tip(row[3:], L)
    np.testing.assert_allclose(U[-1, 1:], ref[1:], rtol=2e-3, atol=1e-3 * np.abs(ref).max())


def test_interior_point_moment_counts_once(tmp_path):
    """A moment at an interior node (shared by two elements) and one inside an element."""
    L = 2.0
    for x in (1.0, 0.75):
        U = _run(tmp_path, LINEAR[0], _line(4, L), [[x, 0, 0, 0, 0, 0, 0, 0, 5.0]])
        # the beam beyond x rotates rigidly: theta_z = M x / EI, v(L) = M x^2 / 2EI + theta (L - x)
        th = 5.0 * x / (E * IZ)
        assert U[-1, 5] == pytest.approx(th, rel=1e-9)
        assert U[-1, 1] == pytest.approx(5.0 * x**2 / (2 * E * IZ) + th * (L - x), rel=1e-9)


def test_euler_bernoulli_uniform_line_load_is_exact_at_nodes(tmp_path):
    """Consistent line loads make Hermite elements exact at the nodes: w(L) = q L^4 / 8EI."""
    L, q = 2.0, 30.0
    dist = [[x, 0, 0, 0, 0, q, 0, 0, 0] for x in np.linspace(0, L, 5)]
    U = _run(tmp_path, LINEAR[0], _line(2, L), [], dist=dist)
    assert U[-1, 2] == pytest.approx(q * L**4 / (8 * E * IY), rel=1e-9)
    assert U[-1, 4] == pytest.approx(-q * L**3 / (6 * E * IY), rel=1e-9)


def test_distributed_moment_on_timoshenko(tmp_path):
    """A uniform distributed m_y: theta_y(L) = m L^2 / 2EI, w(L) = -m L^3 / 3EI."""
    L, m = 2.0, 4.0
    dist = [[x, 0, 0, 0, 0, 0, 0, m, 0] for x in np.linspace(0, L, 5)]
    U = _run(tmp_path, LINEAR[1], _line(32, L), [], dist=dist)
    assert U[-1, 4] == pytest.approx(m * L**2 / (2 * E * IY), rel=2e-3)
    assert U[-1, 2] == pytest.approx(-m * L**3 / (3 * E * IY), rel=2e-3)


@pytest.mark.parametrize("etype", LINEAR[:2])
def test_turned_cantilever_end_moment(tmp_path, etype):
    """An end moment on a round-section beam along d = (1, 2, -2) / 3: theta = M L / EI, w = M L^2 / 2EI."""
    L = 2.0
    d = np.array([1.0, 2.0, -2.0]) / 3.0
    M = np.cross(d, [0.0, 0.0, 1.0])
    M *= 6.0 / np.linalg.norm(M)
    coords = _line(8, L, d)
    U = _run(tmp_path, etype, coords, [np.r_[coords[-1], 0, 0, 0, M]], iz=IY)
    np.testing.assert_allclose(U[-1, 3:], M * L / (E * IY), rtol=1e-5, atol=1e-5 * np.linalg.norm(M) * L / (E * IY))
    # the tip moves along M x d by M L^2 / 2EI
    np.testing.assert_allclose(U[-1, :3], np.cross(M, d) * L**2 / (2 * E * IY), rtol=1e-5,
                               atol=1e-6 * np.linalg.norm(M) * L**2 / (E * IY))
