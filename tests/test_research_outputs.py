"""Schema checks for run manifest, inner solve history, and join keys to Newton CSV."""

import csv
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from processing.static.operations.solver import SolveCondensedSystem
from workflow_orchestrator.run_manifest import write_run_manifest


def test_write_run_manifest_schema(tmp_path: Path):
    job_dir = tmp_path / "jobs" / "case_a"
    job_dir.mkdir(parents=True)
    settings_txt = job_dir / "simulation_settings.txt"
    settings_txt.write_text("[Simulation]\n[Type]\nstatic_nonlinear\n", encoding="utf-8")

    results = tmp_path / "out"
    results.mkdir()
    (results / "logs").mkdir(parents=True)
    (results / "logs" / "newton_history.csv").write_text("h\n1\n", encoding="utf-8")

    p = write_run_manifest(
        job_results_dir=results,
        job_name="case_a",
        job_dir=job_dir,
        wall_time_sec=12.5,
        simulation_settings={"type": "static_nonlinear"},
    )
    assert p is not None and p.is_file()
    with open(p, encoding="utf-8") as f:
        m = json.load(f)
    assert m["job_name"] == "case_a"
    assert m["wall_time_sec"] == 12.5
    assert m["simulation_settings_txt_sha256"] is not None
    assert len(m["simulation_settings_txt_sha256"]) == 64
    assert m["simulation_settings_resolved"]["type"] == "static_nonlinear"
    assert m["paths"]["newton_history_csv"] is not None


def test_inner_solve_history_join_keys(tmp_path: Path):
    """Condensed solve writes ``inner_solve_history.csv`` with NR join keys."""
    K = sp.eye(4, format="csr", dtype=np.float64)
    F = np.ones(4, dtype=np.float64)
    jr = tmp_path / "primary_results"
    jr.mkdir(parents=True)

    solver = SolveCondensedSystem(
        K_cond=K,
        F_cond=F,
        solver_name="direct",
        job_results_dir=str(jr),
        preconditioner=None,
        tolerance=1e-10,
        max_iterations=10,
    )
    U = solver.solve(
        load_increment_index=1,
        newton_iter=2,
        load_factor=0.5,
    )
    assert U is not None
    log_dir = tmp_path / "logs"
    csv_path = log_dir / "inner_solve_history.csv"
    assert csv_path.is_file()
    with open(csv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert len(rows) == 2
    header, data = rows
    assert "load_increment_index" in header
    assert "newton_iter" in header
    assert "fallback_superlu" in header
    idx_inc = header.index("load_increment_index")
    idx_nr = header.index("newton_iter")
    idx_fb = header.index("fallback_superlu")
    assert data[idx_inc] == "1"
    assert data[idx_nr] == "2"
    assert data[idx_fb] == "0"


@pytest.mark.parametrize("solver", ["cg", "direct"])
def test_inner_solve_history_row_exists(tmp_path: Path, solver: str):
    np.random.seed(0)
    n = 32
    A = sp.random(n, n, density=0.15, format="csr", dtype=np.float64)
    K = (A + A.T) / 2 + sp.eye(n, format="csr") * float(n)
    F = np.random.randn(n)
    jr = tmp_path / "res"
    jr.mkdir()

    sol = SolveCondensedSystem(
        K_cond=K,
        F_cond=F,
        solver_name=solver,
        job_results_dir=str(jr),
        preconditioner="jacobi" if solver == "cg" else None,
        tolerance=1e-6,
        max_iterations=200,
    )
    out = sol.solve()
    assert out is not None
    p = tmp_path / "logs" / "inner_solve_history.csv"
    assert p.is_file()
    with open(p, newline="", encoding="utf-8") as f:
        r = list(csv.DictReader(f))
    assert len(r) == 1
    assert int(r[0]["fallback_superlu"]) in (0, 1)
