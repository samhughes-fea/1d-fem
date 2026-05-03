"""Write ``logs/run_manifest.json`` for reproducibility and cross-run comparison."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from pre_processing.parsing.simulation_settings_resolution import finalize_simulation_settings

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore

try:
    import scipy
except ImportError:
    scipy = None


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _git_commit(repo: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return None


def _json_safe(obj: Any) -> Any:
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if np is not None:
        if isinstance(obj, (np.floating, np.integer)):
            return float(obj) if isinstance(obj, np.floating) else int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
    if isinstance(obj, Mapping):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(x) for x in obj]
    return str(obj)


def write_run_manifest(
    *,
    job_results_dir: str | Path,
    job_name: str,
    job_dir: str | Path,
    wall_time_sec: float,
    simulation_settings: Mapping[str, Any] | None = None,
) -> Path | None:
    """
    Write ``<job_results_dir>/logs/run_manifest.json``.

    Returns the written path, or None if the write failed.
    """
    root = Path(job_results_dir)
    logs = root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    out_path = logs / "run_manifest.json"

    repo = _repo_root()
    settings_txt = Path(job_dir) / "simulation_settings.txt"
    sha256 = None
    if settings_txt.is_file():
        try:
            h = hashlib.sha256()
            with open(settings_txt, "rb") as f:
                for chunk in iter(lambda: f.read(65536), b""):
                    h.update(chunk)
            sha256 = h.hexdigest()
        except OSError:
            sha256 = None

    nh = root / "logs" / "newton_history.csv"
    ish = root / "logs" / "inner_solve_history.csv"
    ps = root / "primary_results" / "primary_summary.csv"
    pap = root / "logs" / "primary_artifacts.json"

    manifest: dict[str, Any] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "job_name": job_name,
        "results_directory": str(root.resolve()),
        "job_input_directory": str(Path(job_dir).resolve()),
        "wall_time_sec": float(wall_time_sec),
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "git_commit": _git_commit(repo),
        "numpy_version": getattr(np, "__version__", None) if np else None,
        "scipy_version": getattr(scipy, "__version__", None) if scipy else None,
        "simulation_settings_txt_sha256": sha256,
        "simulation_settings_txt_path": str(settings_txt.resolve())
        if settings_txt.is_file()
        else None,
        "paths": {
            "newton_history_csv": str(nh.resolve()) if nh.is_file() else None,
            "inner_solve_history_csv": str(ish.resolve()) if ish.is_file() else None,
            "primary_summary_csv": str(ps.resolve()) if ps.is_file() else None,
            "primary_artifacts_json": str(pap.resolve()) if pap.is_file() else None,
            "process_job_log": str((root / "logs" / "process_job.log").resolve())
            if (root / "logs" / "process_job.log").is_file()
            else None,
        },
    }
    if simulation_settings is not None:
        resolved = dict(simulation_settings)
        finalize_simulation_settings(resolved, type_line_explicit=("type" in resolved))
        manifest["simulation_settings_resolved"] = _json_safe(resolved)

    try:
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)
    except OSError:
        return None
    return out_path
