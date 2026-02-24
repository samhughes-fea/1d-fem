# post_processing/graphical_visualisers/job_discovery_utils.py

"""
Shared job result discovery and mesh path resolution for post-processing visualisers.

Result directories under post_processing/results are named:
  {case_name}_{timestamp}_pid{pid}_{uid}
where case_name is the job input folder name (e.g. job_0000 or job_0000_n10).
This module parses those names and resolves mesh paths to the correct job input dir.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Tuple

# Matches job_0000_ts... and job_0000_n10_ts... (optional _nN suffix after numeric id).
# Captures: id (digits), optional suffix (_n10 etc), timestamp (rest until end).
_JOB_RESULT_DIR_PATTERN = re.compile(
    r"^job_(?P<id>\d+)(?:(?P<suffix>_n\d+))?_(?P<ts>[\d\-_]+_pid\d+_[a-f0-9]+)$"
)


def parse_job_result_dir_name(dir_name: str) -> Optional[Tuple[str, str, str]]:
    """
    Parse a result directory name under post_processing/results.

    Supports both legacy (job_0000_...) and current (job_0000_n10_...) naming.

    Parameters
    ----------
    dir_name : str
        Basename of the result directory, e.g. "job_0000_n10_2026-02-23_12-00-00_pid1234_abc12345".

    Returns
    -------
    Optional[Tuple[str, str, str]]
        If the name matches: (job_id, job_folder_name, timestamp).
        - job_id: numeric id string, e.g. "0000".
        - job_folder_name: full job input folder name for mesh lookup, e.g. "job_0000_n10".
        - timestamp: remainder for filenames, e.g. "2026-02-23_12-00-00_pid1234_abc12345".
        Returns None if the name does not match.
    """
    m = _JOB_RESULT_DIR_PATTERN.match(dir_name.strip())
    if not m:
        return None
    job_id = m.group("id")
    suffix = m.group("suffix") or ""
    job_folder_name = f"job_{job_id}{suffix}"
    timestamp = m.group("ts")
    return (job_id, job_folder_name, timestamp)


def get_mesh_paths_for_result_dir(
    result_dir_name: str, jobs_dir: Path
) -> Optional[Tuple[Path, Path]]:
    """
    Resolve grid.txt and element.txt paths for a result directory name.

    Parameters
    ----------
    result_dir_name : str
        Basename of the result directory (e.g. job_0000_n10_2026-...).
    jobs_dir : Path
        Path to the jobs directory (project root / "jobs").

    Returns
    -------
    Optional[Tuple[Path, Path]]
        (grid_file_path, element_file_path) for the corresponding job input dir,
        or None if the result dir name is unrecognised.
    """
    parsed = parse_job_result_dir_name(result_dir_name)
    if not parsed:
        return None
    _job_id, job_folder_name, _ts = parsed
    job_input_dir = jobs_dir / job_folder_name
    return (job_input_dir / "grid.txt", job_input_dir / "element.txt")
