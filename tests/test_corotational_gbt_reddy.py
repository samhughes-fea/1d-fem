"""Corotational GBT, Reddy and Levinson beams.

References:

- At rest the tangent is the linear core's stiffness, with the third-order beams' nodal
  shear strains turned into the core's slopes (``v' = gamma_xy + theta_z``,
  ``w' = gamma_xz - theta_y``).
- The tangent is the derivative of the internal force (central differences).
- Objectivity: a rigid rotation of a deformed state turns the internal forces with it and
  leaves the nodal shear strains' forces unchanged.
- A cantilever under a tip moment rolls up as the corotational Timoshenko beam does
  (pure bending has no shear, so the third-order terms are idle).
- A GBT box under a small load is the linear GBT box; at large deflection its stiffening
  is that of the corotational Timoshenko and Benscoter boxes.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

from onedfem.pre_processing.element_library.nonlinear.large_rotations.corotational.corotational_section import (
    rotation_exp,
    rotation_log,
)
from tests.test_corotational_frames import RECT, _JOB, _newton, _tip
from tests.test_element_frames import Mesh, _line

CR_RD, CR_LV, CR_TS = ("CorotationalReddyBeamElement3D", "CorotationalLevinsonBeamElement3D",
                       "CorotationalTimoshenkoBeamElement3D")


def _mesh(etype, coords):
    mesh = Mesh(etype, coords, props=RECT)
    for el in mesh.elements:
        el.job_results_dir = _JOB
    return mesh


def _element(etype, L=0.25):
    el = _mesh(etype, _line([0, 0, 0], [1, 0, 0], L, 1)).elements[0]
    el._cr_setup()
    return el


@pytest.mark.parametrize("etype", [CR_RD, CR_LV])
def test_tangent_at_rest_is_the_core_with_shear_strain_dofs(etype):
    el = _element(etype)
    T = np.eye(16)
    T[12, 5] = T[14, 11] = 1.0          # v' = gamma_xy + theta_z
    T[13, 4] = T[15, 10] = -1.0         # w' = gamma_xz - theta_y
    K = np.asarray(el.element_stiffness_matrix().K_e)
    np.testing.assert_allclose(el.tangent_stiffness_matrix(np.zeros(16)), T.T @ K @ T,
                               atol=1e-9 * np.abs(K).max())


@pytest.mark.parametrize("etype", [CR_RD, CR_LV])
def test_tangent_is_the_derivative_of_the_internal_force(etype):
    el = _element(etype)
    U = np.random.default_rng(1).normal(size=16) * 0.2
    h = 1e-6
    Kf = np.column_stack([(el.internal_force_vector(U + h * e) - el.internal_force_vector(U - h * e)) / (2 * h)
                          for e in np.eye(16)])
    np.testing.assert_allclose(el.tangent_stiffness_matrix(U), Kf, atol=1e-6 * np.abs(Kf).max())


@pytest.mark.parametrize("etype", [CR_RD, CR_LV])
def test_rigid_rotation_turns_the_forces(etype):
    el = _element(etype)
    X = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
    U = np.random.default_rng(2).normal(size=16) * 0.05
    Q = rotation_exp(np.array([0.7, -1.1, 0.4]))
    V = U.copy()
    for i, (t, r) in enumerate(((slice(0, 3), slice(3, 6)), (slice(6, 9), slice(9, 12)))):
        V[t] = Q @ (X[i] + U[t]) - X[i]
        V[r] = rotation_log(Q @ rotation_exp(U[r]))
    f, g = el.internal_force_vector(U), el.internal_force_vector(V)
    scale = np.abs(f).max()
    for s in (slice(0, 3), slice(6, 9)):
        np.testing.assert_allclose(g[s], Q @ f[s], atol=1e-8 * scale)
    np.testing.assert_allclose(g[12:], f[12:], atol=1e-8 * scale)
    # moments: work-conjugate to the rotation vector, so compare the spatial moments
    from onedfem.processing.static.follower_loads import tangent_operator

    for r in (slice(3, 6), slice(9, 12)):
        m_u = np.linalg.solve(tangent_operator(U[r]).T, f[r])
        m_v = np.linalg.solve(tangent_operator(V[r]).T, g[r])
        np.testing.assert_allclose(m_v, Q @ m_u, atol=1e-8 * scale)


@pytest.mark.parametrize("etype", [CR_RD, CR_LV])
def test_tip_moment_rolls_up_as_the_timoshenko_beam(etype):
    from tests.test_second_order import E

    n = 16
    M = E * RECT["I_z"] * np.pi / 2           # a quarter turn of the tip
    ref = _mesh(CR_TS, _line([0, 0, 0], [1, 0, 0], 1.0, n))
    third = _mesh(etype, _line([0, 0, 0], [1, 0, 0], 1.0, n))
    U0 = _newton(ref, _tip(ref, n, [0, 0, 0, 0, 0, M]), list(range(ref.dpn)), 4).reshape(-1, ref.dpn)
    U1 = _newton(third, _tip(third, n, [0, 0, 0, 0, 0, M]), list(range(third.dpn)), 4).reshape(-1, third.dpn)
    assert U1[-1, 5] == pytest.approx(np.pi / 2, rel=3e-3)
    np.testing.assert_allclose(U1[:, :6], U0, atol=1e-6)


# --- whole jobs: the corotational GBT box cantilever ----------------------------------

_NL = ("[Simulation]\n[Type]\nstatic_nonlinear\n\n[Newton]\ntolerance = 1e-6\nrelative_tolerance = 1e-10\n"
       "relative_reference = first_residual\nmax_iterations = 40\ntolerance_delta_u = 1e-14\n\n[Solver]\n"
       "type = direct\ntolerance = 1e-12\n\n[Nonlinear]\nnum_increments = {n}\nline_search = false\n")


def _box_job(tmp_path, etype, load, nonlinear, increments=8):
    import csv
    import logging

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples", "thin_walled"))
    import common

    from onedfem.workflow_orchestrator.run_job import process_job

    # stocky walls (b/t = 7.5): at the large-deflection load the Brazier flattening of
    # CorotationalGBT stays small, so the test isolates the large-rotation kinematics
    b, h, t = 0.06, 0.06, 8e-3
    pts = {1: (-b / 2, -h / 2), 2: (b / 2, -h / 2), 3: (b / 2, h / 2), 4: (-b / 2, h / 2)}
    walls = [(1, 1, 2, 1), (2, 2, 3, 1), (3, 3, 4, 1), (4, 4, 1, 1)]
    mat = {1: common.Material(70e9, 70e9, 26.9e9, 0.3, 26.9e9, 26.9e9, 2700.0)}
    sec = common.SectionSpec(pts, walls, mat, {1: [(1, 0.0, t)]}, n_sub=2)
    gbt = "GBT" in etype
    spec = common.CantileverSpec(1.5, 12, etype, sec, tip_load=load, modes=2 if gbt else 0, warping=1 if gbt else 0)
    name = f"{etype}_{load[2]:g}_{int(nonlinear)}"
    jd = str(tmp_path / "jobs" / name)
    common.write_job(jd, spec)
    common.write_supports(jd, spec)
    if nonlinear:
        with open(os.path.join(jd, "simulation_settings.txt"), "w", encoding="utf-8") as f:
            f.write(_NL.format(n=increments))
    res = common.results_dir(str(tmp_path / "res" / name))
    logging.getLogger().setLevel(logging.WARNING)
    process_job(jd, res, {}, {}, force_serial=True, max_processes_per_job=1)
    if nonlinear:
        with open(os.path.join(res, "logs", "newton_history.csv"), newline="", encoding="utf-8") as fh:
            last = list(csv.DictReader(fh))[-1]
        assert float(last["load_factor"]) == pytest.approx(1.0) and last["residual_ok"] == "1", last
    with open(os.path.join(res, "primary_results", "global", "U_global.csv"), newline="", encoding="utf-8") as f:
        U = np.array([float(r[1]) for r in list(csv.reader(f))[1:]])
    return U.reshape(spec.n_el + 1, -1)


def test_corotational_gbt_job_small_load_is_the_linear_gbt(tmp_path):
    load = [0.0, 0.0, 12.0, 0.0, 0.0, 0.0]
    lin = _box_job(tmp_path, "LinearGBTBeamElement3D", load, False)
    cr = _box_job(tmp_path, "CorotationalGBTBeamElement3D", load, True, increments=1)
    # the two modes' amplitudes and slopes (last four columns) also take the Brazier
    # flattening, which is quadratic in the load and absent from the linear element
    beam = slice(0, cr.shape[1] - 4)
    assert np.abs(cr[:, beam] - lin[:, beam]).max() < 1e-3 * np.abs(lin).max()
    assert np.abs(cr[:, -4:]).max() < 1e-2 * np.abs(lin).max()


def test_corotational_gbt_job_large_deflection_follows_the_timoshenko_beam(tmp_path):
    """At a tip deflection of ~14 % of the span the corotational GBT box shortens its
    deflection relative to its linear solution as the corotational Timoshenko beam does."""
    load = [0.0, 0.0, 1.4e4, 0.0, 0.0, 0.0]
    ratio = {}
    for core in ("GBT", "Timoshenko"):
        lin = _box_job(tmp_path, f"Linear{core}BeamElement3D", load, False)
        cr = _box_job(tmp_path, f"Corotational{core}BeamElement3D", load, True)
        ratio[core] = cr[-1, 2] / lin[-1, 2]
        assert cr[-1, 0] < 0.0  # the tip draws in
    assert ratio["GBT"] < 0.995
    assert ratio["GBT"] == pytest.approx(ratio["Timoshenko"], abs=5e-3)
