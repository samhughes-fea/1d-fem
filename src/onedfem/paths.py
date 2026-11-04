"""Where 1d-fem reads jobs and writes results when no directory is given.

In a source checkout this is the repository root (``jobs/`` and ``results/``). For an
installed package it is the working directory, unless ``ONEDFEM_PROJECT_ROOT`` is set.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["project_root", "jobs_root", "results_root"]


def project_root() -> Path:
    env = os.environ.get("ONEDFEM_PROJECT_ROOT")
    if env:
        return Path(env).resolve()
    checkout = Path(__file__).resolve().parents[2]
    if (checkout / "pyproject.toml").is_file() and (checkout / "src" / "onedfem").is_dir():
        return checkout
    return Path.cwd()


def jobs_root() -> Path:
    return project_root() / "jobs"


def results_root() -> Path:
    return project_root() / "results"
