"""End-to-end smoke: canonical ``jobs/job_smoke_eigen`` and ``jobs/job_smoke_buckling``."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from workflow_orchestrator.run_job import process_job, setup_job_results_directory


@pytest.mark.integration
@pytest.mark.parametrize(
    "job_folder",
    ["job_smoke_eigen", "job_smoke_buckling"],
)
def test_canonical_smoke_job_process_job(job_folder: str) -> None:
    job_dir = PROJECT_ROOT / "jobs" / job_folder
    assert job_dir.is_dir(), f"Missing job directory: {job_dir}"

    case_name = f"pytest_{job_folder}"
    res_dir = setup_job_results_directory(case_name)
    jt: dict = {}
    je: dict = {}
    process_job(
        str(job_dir),
        res_dir,
        jt,
        je,
        force_serial=True,
        max_processes_per_job=1,
    )
    assert (Path(res_dir) / "logs" / "run_manifest.json").is_file()
    assert (Path(res_dir) / "primary_results").is_dir()
