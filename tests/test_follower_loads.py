"""Follower loads (processing/static/follower_loads.py and follower_load.txt)."""

import importlib.util
import os

import numpy as np
import pytest

from processing.static.follower_loads import FollowerLoads, parse_follower_load, rotation_exp, tangent_operator

ROOT = os.path.join(os.path.dirname(__file__), "..")


def test_tangent_operator_is_the_derivative_of_the_exponential():
    th = np.array([0.4, -1.1, 0.7])
    J = tangent_operator(th)
    h = 1e-6
    R = rotation_exp(th)
    for j, e in enumerate(np.eye(3) * h):
        dR = (rotation_exp(th + e) - rotation_exp(th - e)) / (2 * h)
        W = dR @ R.T
        assert np.allclose([W[2, 1], W[0, 2], W[1, 0]], J[:, j], atol=1e-8)


def test_load_stiffness_is_minus_the_derivative_of_the_load():
    xyz = np.array([[0.0, 0, 0], [1.0, 0, 0]])
    fl = FollowerLoads([[1.0, 0, 0, 3.0, -2.0, 1.0, 0.5, 0.2, -0.7]], xyz, 6, 12)
    rng = np.random.default_rng(1)
    U = rng.normal(size=12) * 0.5
    K = fl.stiffness(U).toarray()
    h = 1e-6
    D = np.column_stack([(fl.forces(U + e) - fl.forces(U - e)) / (2 * h) for e in np.eye(12) * h])
    assert np.abs(K + D).max() < 1e-6 * np.abs(K).max()
    # undeformed: the load is the given one
    assert np.allclose(fl.forces(np.zeros(12))[6:], [3.0, -2.0, 1.0, 0.5, 0.2, -0.7])


def test_parser(tmp_path):
    p = tmp_path / "follower_load.txt"
    p.write_text("# tip\n[Follower load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n1 0 0 0 5 0 0 0 0\n")
    assert np.allclose(parse_follower_load(str(p)), [[1, 0, 0, 0, 5, 0, 0, 0, 0]])


def _example():
    spec = importlib.util.spec_from_file_location("follower_run", os.path.join(ROOT, "examples", "follower_loads",
                                                                                "run.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("follower", [True, False])
def test_cantilever_matches_the_elastica(tmp_path, follower):
    """Tip force P L^2 / EI = 3 staying normal to the tip (follower) or vertical (dead), 40
    corotational Timoshenko elements (the error falls as h^4: 6e-3 L at 20 elements)."""
    import logging

    mod = _example()
    level = logging.getLogger().level
    try:
        U = mod.run_job(mod.write_job(str(tmp_path / "job"), 3.0, follower), str(tmp_path))
    finally:
        logging.getLogger().setLevel(level)
    _, x, y, th = mod.elastica(3.0, follower)
    assert U[-1, 0] == pytest.approx(x[-1] - mod.L, abs=5e-4 * mod.L)
    assert U[-1, 1] == pytest.approx(y[-1], abs=5e-4 * mod.L)
    assert U[-1, 5] == pytest.approx(th[-1], abs=1e-3)


def test_follower_load_needs_a_nonlinear_job(tmp_path, capfd):
    import shutil

    from workflow_orchestrator.run_job import process_job

    job = tmp_path / "job"
    shutil.copytree(os.path.join(ROOT, "jobs", "job_smoke_eigen"), job)
    (job / "simulation_settings.txt").write_text("[Simulation]\n[Type]\nstatic\n")
    (job / "follower_load.txt").write_text("[Follower load]\n2 0 0 0 1 0 0 0 0\n")
    res = tmp_path / "res"
    for sub in ("element_stiffness_matrices", "element_force_vectors", "primary_results", "secondary_results",
                "logs", "maps", "diagnostics"):
        os.makedirs(res / sub)
    process_job(str(job), str(res), {}, {}, force_serial=True, max_processes_per_job=1)
    out = capfd.readouterr()
    assert "needs a nonlinear static job" in out.out + out.err
    assert not (res / "primary_results" / "global" / "U_global.csv").exists()
