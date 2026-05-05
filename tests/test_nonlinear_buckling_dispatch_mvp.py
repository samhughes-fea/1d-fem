"""Nonlinear buckling MVP: continuation outputs and dispatch wiring."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from simulation_runner.buckling.nonlinear_buckling_simulation import NonlinearBucklingSimulationRunner
from workflow_orchestrator.run_job import process_job, setup_job_results_directory


def test_nonlinear_buckling_runner_requires_full_settings(tmp_path: Path) -> None:
    jrd = tmp_path / "job_out"
    jrd.mkdir()
    settings = {"job_results_dir": str(jrd), "simulation_settings": {}}
    with pytest.raises(KeyError):
        NonlinearBucklingSimulationRunner(settings, "mvp").run()


@pytest.mark.integration
def test_process_job_nonlinear_buckling_flag_dispatches_continuation(tmp_path: Path) -> None:
    src = PROJECT_ROOT / "jobs" / "job_smoke_buckling"
    job_copy = tmp_path / "job_nl_buck_mvp"
    shutil.copytree(src, job_copy)
    sim = job_copy / "simulation_settings.txt"
    content = sim.read_text(encoding="utf-8").rstrip() + "\nnonlinear_buckling = true\n"
    sim.write_text(content, encoding="utf-8")

    res_dir = setup_job_results_directory("pytest_nl_buck_mvp_dispatch")
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
    assert history.is_file(), "nonlinear buckling continuation should write history CSV"
    assert summary.is_file(), "nonlinear buckling continuation should write summary JSON"


def test_deprecated_runner_aliases_are_subclasses() -> None:
    from simulation_runner.buckling.buckling_simulation import (
        BucklingSimulationRunner,
        LinearBucklingSimulationRunner,
    )
    from simulation_runner.transient.dynamic_simulation import (
        DynamicSimulationRunner,
        TransientSimulationRunner,
    )

    assert issubclass(BucklingSimulationRunner, LinearBucklingSimulationRunner)
    assert issubclass(DynamicSimulationRunner, TransientSimulationRunner)
