"""Machine-readable manifest of primary result files (Sections 2–5)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

SCHEMA_VERSION = "1.0"


def write_primary_artifact_manifest(
    job_results_root: str | Path,
    *,
    family: str,
    job_name: str,
    artifacts: Mapping[str, str],
) -> Path | None:
    """
    Write ``<job_results_root>/logs/primary_artifacts.json``.

    *artifacts* maps logical keys to paths **relative to** *job_results_root* (POSIX slashes).
    Returns the written path, or None on failure.
    """
    root = Path(job_results_root)
    logs = root / "logs"
    try:
        logs.mkdir(parents=True, exist_ok=True)
        out = logs / "primary_artifacts.json"
        payload = {
            "schema_version": SCHEMA_VERSION,
            "family": str(family),
            "job_name": str(job_name),
            "artifacts": {str(k): str(Path(v).as_posix()) for k, v in artifacts.items()},
        }
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    except OSError:
        return None
    return out
