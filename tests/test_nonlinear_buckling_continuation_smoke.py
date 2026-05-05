from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from workflow_orchestrator.run_job import process_job, setup_job_results_directory
from processing.common.nonlinear_buckling_continuation import (
    ContinuationConfig,
    ImperfectionConfig,
    predictor_load_factor_for_increment,
    seed_initial_imperfection,
)
from processing.common.nonlinear_equilibrium import solve_arc_length_corrector
from simulation_runner.buckling.nonlinear_buckling_simulation import NonlinearBucklingSimulationRunner


@pytest.mark.integration
def test_nonlinear_buckling_continuation_smoke_outputs_history(tmp_path: Path) -> None:
    src = PROJECT_ROOT / "jobs" / "job_smoke_buckling"
    job_copy = tmp_path / "job_nl_buckling_smoke"
    shutil.copytree(src, job_copy)

    element = job_copy / "element.txt"
    element.write_text(
        element.read_text(encoding="utf-8").replace(
            "LinearEulerBernoulliBeamElement3D",
            "NonlinearEulerBernoulliBeamElement3D",
        ),
        encoding="utf-8",
    )

    sim = job_copy / "simulation_settings.txt"
    sim.write_text(
        sim.read_text(encoding="utf-8").rstrip()
        + "\nnonlinear_buckling = true\nnum_increments = 2\n",
        encoding="utf-8",
    )

    res_dir = setup_job_results_directory("pytest_nl_buckling_continuation_smoke")
    process_job(
        str(job_copy),
        res_dir,
        {},
        {},
        force_serial=True,
        max_processes_per_job=1,
    )

    history = Path(res_dir) / "primary_results" / "nonlinear_buckling_results" / "continuation_history.csv"
    text = history.read_text(encoding="utf-8")
    assert "increment_index" in text
    assert "load_factor" in text


@pytest.mark.integration
def test_nonlinear_buckling_arc_length_mode_writes_predictor_metadata(tmp_path: Path) -> None:
    src = PROJECT_ROOT / "jobs" / "job_smoke_buckling"
    job_copy = tmp_path / "job_nl_buckling_arc"
    shutil.copytree(src, job_copy)

    element = job_copy / "element.txt"
    element.write_text(
        element.read_text(encoding="utf-8").replace(
            "LinearEulerBernoulliBeamElement3D",
            "NonlinearEulerBernoulliBeamElement3D",
        ),
        encoding="utf-8",
    )

    sim = job_copy / "simulation_settings.txt"
    sim.write_text(
        sim.read_text(encoding="utf-8").rstrip()
        + "\nnonlinear_buckling = true\ncontinuation_method = arc_length\narc_length_radius = 0.01\nnum_increments = 2\n",
        encoding="utf-8",
    )

    res_dir = setup_job_results_directory("pytest_nl_buckling_arc_length_smoke")
    process_job(
        str(job_copy),
        res_dir,
        {},
        {},
        force_serial=True,
        max_processes_per_job=1,
    )

    history = Path(res_dir) / "primary_results" / "nonlinear_buckling_results" / "continuation_history.csv"
    text = history.read_text(encoding="utf-8")
    assert "continuation_method" in text
    assert "predictor_load_factor" in text


def test_nonlinear_buckling_imperfection_seed_from_linear_mode_shape(tmp_path: Path) -> None:
    jrd = tmp_path / "job_out"
    modal_dir = jrd / "primary_results" / "modal_results"
    modal_dir.mkdir(parents=True)
    (modal_dir / "seed_job_buckling_mode_shapes.txt").write_text("0.0\n0.2\n0.4\n", encoding="utf-8")

    runner = NonlinearBucklingSimulationRunner(
        {
            "job_results_dir": str(jrd),
            "simulation_settings": {
                "buckling": {
                    "imperfection_source": "linear_buckling",
                    "imperfection_mode_index": 0,
                    "imperfection_scale": 0.5,
                }
            },
        },
        "seed_job",
    )
    seeded, meta = seed_initial_imperfection(
        job_name="seed_job",
        primary_results_dir=str(jrd / "primary_results"),
        job_results_dir=str(jrd),
        U_global=np.zeros(3, dtype=float),
        config=ImperfectionConfig(source="linear_buckling", mode_index=0, scale=0.5),
    )
    assert np.allclose(seeded, np.array([0.0, 0.1, 0.2]))
    assert meta is not None
    assert meta["mode_index"] == 0


def test_predictor_load_factor_helper_supports_arc_length_metadata() -> None:
    cfg = ContinuationConfig(
        continuation_method="arc_length",
        load_factors=np.array([0.5, 1.0], dtype=float),
        arc_length_radius=0.1,
        arc_length_alpha_scale=1.0,
    )
    predictor = predictor_load_factor_for_increment(
        config=cfg,
        increment_index=0,
        current_U=np.zeros(3, dtype=float),
        tip_dof=1,
        reference_load_vector=np.array([0.0, 1.0, 0.0], dtype=float),
    )
    assert predictor >= 0.0


def test_arc_length_corrector_solver_converges_on_trivial_state() -> None:
    def build_system_from_state(U_state: np.ndarray, load_factor_state: float):
        residual = np.zeros_like(U_state)
        return 0.0, 0.0, residual

    def solve_condensed_step_from_state(iteration: int, load_factor_state: float) -> np.ndarray:
        return np.zeros(2, dtype=float)

    def reconstruct_delta(delta_u_cond: np.ndarray) -> np.ndarray:
        return np.asarray(delta_u_cond, dtype=float)

    seen = []

    result = solve_arc_length_corrector(
        U_prev=np.zeros(2, dtype=float),
        load_factor_prev=0.0,
        predictor_displacement=np.array([0.1, 0.0], dtype=float),
        reference_load_vector=np.array([0.0, 1.0], dtype=float),
        arc_length_radius=0.2,
        alpha_scale=1.0,
        newton_tol=1e-8,
        newton_max_iter=3,
        build_system_from_state=build_system_from_state,
        solve_condensed_step_from_state=solve_condensed_step_from_state,
        reconstruct_delta=reconstruct_delta,
        iteration_callback=lambda rec: seen.append(rec),
        load_increment_index=1,
    )
    assert result.converged is True
    assert result.iterations_used >= 1
    assert len(seen) >= 1


@pytest.mark.integration
def test_imperfect_column_benchmark_job_writes_acceptance_artifacts(tmp_path: Path) -> None:
    src = PROJECT_ROOT / "jobs" / "job_benchmark_nl_buckling_imperfect_column"
    job_copy = tmp_path / "job_benchmark_nl_buckling_imperfect_column"
    shutil.copytree(src, job_copy)

    modal_dir = job_copy / "primary_results" / "modal_results"
    modal_dir.mkdir(parents=True)
    (modal_dir / "job_benchmark_nl_buckling_imperfect_column_buckling_mode_shapes.txt").write_text(
        "0.0\n1.0\n0.0\n0.0\n0.0\n0.0\n0.0\n0.5\n0.0\n0.0\n0.0\n0.0\n0.0\n0.25\n0.0\n0.0\n0.0\n0.0\n0.0\n0.125\n0.0\n0.0\n0.0\n0.0\n0.0\n0.0625\n0.0\n0.0\n0.0\n0.0\n",
        encoding="utf-8",
    )

    sim = job_copy / "simulation_settings.txt"
    sim.write_text(
        sim.read_text(encoding="utf-8").rstrip()
        + "\ncontinuation_method = arc_length\n",
        encoding="utf-8",
    )

    res_dir = setup_job_results_directory("pytest_nl_buckling_imperfect_column_benchmark")
    process_job(
        str(job_copy),
        res_dir,
        {},
        {},
        force_serial=True,
        max_processes_per_job=1,
    )

    history = Path(res_dir) / "primary_results" / "nonlinear_buckling_results" / "continuation_history.csv"
    summary = Path(res_dir) / "diagnostics" / "nonlinear_buckling_summary.json"
    assert history.is_file()
    assert summary.is_file()
    htxt = history.read_text(encoding="utf-8")
    stxt = summary.read_text(encoding="utf-8")
    assert ",load_control," in htxt or ",arc_length," in htxt
    assert "history_csv" in stxt
    rows = htxt.strip().splitlines()
    assert len(rows) >= 2
    first = rows[1].split(",")
    assert abs(float(first[1]) - 1.0) <= 1.0e-12
    assert abs(float(first[6]) - 0.0) <= 1.0e-12
