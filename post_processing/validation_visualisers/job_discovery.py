# post_processing/validation_visualisers/job_discovery.py
"""
Discover job names from post_processing/results (timestamped result dirs).
Used by run_all_abaqus_jobs.py --from-results and can be reused by comparison scripts.
"""
from __future__ import annotations

from pathlib import Path

from post_processing.validation_visualisers.abaqus.config import (
    FEM_RESULTS_DIR,
    RESULT_DIR_PATTERN,
)


def discover_job_names_from_results(
    fem_results_dir: Path | None = None,
) -> list[str]:
    """
    Scan fem_results_dir for timestamped result directories, extract canonical
    job names job_XXXX_nN, and return sorted unique list.

    Result dirs match: job_<base_id>_n<n>_<timestamp>_pid<pid>_<hash>
    """
    base = Path(fem_results_dir) if fem_results_dir is not None else FEM_RESULTS_DIR
    if not base.is_dir():
        return []
    seen: set[str] = set()
    for p in base.iterdir():
        if not p.is_dir():
            continue
        m = RESULT_DIR_PATTERN.match(p.name)
        if not m:
            continue
        base_id = int(m.group("base_id"))
        n = int(m.group("n"))
        job_name = f"job_{base_id:04d}_n{n}"
        seen.add(job_name)
    return sorted(seen)
