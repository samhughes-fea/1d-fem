# post_processing/validation_visualisers/abaqus/run_abaqus_cae.py
"""
Generate Abaqus script for given job(s) and run it with project Python (abqpy).
This is the project's canonical way to generate Abaqus results: the Python Abaqus
package (abqpy) is used to run the generated CAE script; abqpy's saveAs() launches
Abaqus (abaqus cae noGUI=<script>), which builds the model, runs the job, and
exports ODB to CSV. Writes results to validation_visualisers/abaqus_results/job_XXXX_nX/.
Run from project root. Requires abqpy and Abaqus (ABAQUS_BAT_PATH or on PATH).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = SCRIPT_DIR.parent
PROJECT_ROOT = VALIDATION_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from post_processing.validation_visualisers.abaqus.config import (
    JOBS_DIR,
    ABAQUS_GENERATED_DIR,
    ABAQUS_RESULTS_DIR,
    ABAQUS_CAE_CMD,
    ABAQUS_LAUNCHER_PATH,
)
from post_processing.validation_visualisers.abaqus.job_to_abaqus_script import _parse_job, _generate_script_content


def _abaqus_available() -> bool:
    p = Path(ABAQUS_CAE_CMD)
    if p.is_absolute() and p.is_file():
        return True
    return shutil.which(ABAQUS_CAE_CMD) is not None


def _abqpy_available() -> bool:
    try:
        import driverUtils  # noqa: F401
        return True
    except ImportError:
        return False


def generate_script_only(job_name: str, regenerate: bool = True) -> bool:
    """
    Generate Abaqus CAE script for job_name only (no Abaqus run).
    Return True if script was generated or already exists.
    """
    job_dir = JOBS_DIR / job_name
    if not job_dir.is_dir():
        print(f"Job directory not found: {job_dir}", file=sys.stderr)
        return False

    script_path = ABAQUS_GENERATED_DIR / f"run_{job_name}.py"
    if regenerate or not script_path.exists():
        try:
            data = _parse_job(job_dir)
            out_csv_dir = str(ABAQUS_RESULTS_DIR / job_name)
            content = _generate_script_content(data, out_csv_dir)
            ABAQUS_GENERATED_DIR.mkdir(parents=True, exist_ok=True)
            script_path.write_text(content, encoding="utf-8")
        except Exception as e:
            print(f"Failed to generate script: {e}", file=sys.stderr)
            return False
    return True


def run_job(job_name: str, regenerate: bool = True) -> bool:
    """
    Generate script for job_name (e.g. job_0000_n8), run Abaqus CAE, return True if successful.
    """
    job_dir = JOBS_DIR / job_name
    if not job_dir.is_dir():
        print(f"Job directory not found: {job_dir}", file=sys.stderr)
        return False

    script_path = ABAQUS_GENERATED_DIR / f"run_{job_name}.py"
    if regenerate or not script_path.exists():
        try:
            data = _parse_job(job_dir)
            out_csv_dir = str(ABAQUS_RESULTS_DIR / job_name)
            content = _generate_script_content(data, out_csv_dir)
            ABAQUS_GENERATED_DIR.mkdir(parents=True, exist_ok=True)
            script_path.write_text(content, encoding="utf-8")
        except Exception as e:
            print(f"Failed to generate script: {e}", file=sys.stderr)
            return False

    if not _abqpy_available():
        print(
            "abqpy not found. Install with: pip install abqpy==2021.*",
            file=sys.stderr,
        )
        return False
    if not _abaqus_available():
        print(
            f"Abaqus not found (launcher: {ABAQUS_CAE_CMD}). "
            "Script was generated; set ABAQUS_CAE_ROOT or run manually with project Python.",
            file=sys.stderr,
        )
        return False

    env = os.environ.copy()
    if ABAQUS_LAUNCHER_PATH:
        env["ABAQUS_BAT_PATH"] = ABAQUS_LAUNCHER_PATH
    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(PROJECT_ROOT),
        env=env,
    )
    return result.returncode == 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate and run Abaqus CAE for given job(s)")
    parser.add_argument("--job", type=str, action="append", dest="jobs", help="Job name (e.g. job_0000_n8); repeat for multiple")
    parser.add_argument("--no-regenerate", action="store_true", help="Skip script generation if script already exists")
    args = parser.parse_args()

    if not args.jobs:
        # Default: run job_0000_n8
        args.jobs = ["job_0000_n8"]

    for job_name in args.jobs:
        print(f"Running Abaqus for {job_name} ...")
        ok = run_job(job_name, regenerate=not args.no_regenerate)
        print(f"  {'OK' if ok else 'FAILED'}")
    print("Done.")


if __name__ == "__main__":
    main()
