#!/usr/bin/env python
"""
Batch validation: run Abaqus for a fixed set of jobs (when Abaqus is available),
then run comparison scripts and check that output files exist.

Use for regression checks on a runner with Abaqus license. Run from project root:
  python post_processing/validation_visualisers/run_batch_validation.py
  python post_processing/validation_visualisers/run_batch_validation.py --jobs job_0000_n8 job_0005_n16
  python post_processing/validation_visualisers/run_batch_validation.py --compare-only
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Default jobs for batch validation (point-load and distributed-load coverage)
DEFAULT_JOBS = ["job_0000_n8", "job_0005_n16"]


def _abaqus_can_run() -> bool:
    """True if abqpy and Abaqus launcher are available."""
    try:
        from post_processing.validation_visualisers.abaqus.run_abaqus_cae import (
            _abqpy_available,
            _abaqus_available,
        )
        return _abqpy_available() and _abaqus_available()
    except Exception:
        return False


def _run_abaqus_jobs(job_names: list[str], regenerate: bool = True) -> bool:
    """Run run_abaqus_cae for each job. Return True if all succeeded."""
    from post_processing.validation_visualisers.abaqus.run_abaqus_cae import run_job

    all_ok = True
    for job_name in job_names:
        print(f"Running Abaqus for {job_name} ...")
        ok = run_job(job_name, regenerate=regenerate)
        print(f"  {'OK' if ok else 'FAILED'}")
        if not ok:
            all_ok = False
    return all_ok


def _run_comparisons() -> bool:
    """Run all validation visualisers. Return True if all exited 0."""
    run_all = SCRIPT_DIR / "run_all_validation_visualisers.py"
    if not run_all.is_file():
        print(f"Not found: {run_all}", file=sys.stderr)
        return False
    result = subprocess.run(
        [sys.executable, str(run_all)],
        cwd=str(PROJECT_ROOT),
        capture_output=False,
        timeout=300,
    )
    return result.returncode == 0


def _check_outputs(job_names: list[str]) -> tuple[bool, list[str]]:
    """
    Check that expected output files exist.
    Returns (all_found, list of missing paths).
    """
    abaqus_results = SCRIPT_DIR / "abaqus_results"
    output_dir = SCRIPT_DIR / "output"
    missing = []

    for job_name in job_names:
        job_dir = abaqus_results / job_name
        u_csv = job_dir / "U_global.csv"
        sf_csv = job_dir / "section_forces.csv"
        if not u_csv.is_file():
            missing.append(str(u_csv))
        if not sf_csv.is_file():
            missing.append(str(sf_csv))

    # At least one comparison output should exist if comparisons were run
    if output_dir.is_dir():
        outputs = list(output_dir.glob("*.png")) + list(output_dir.glob("*.csv"))
        if not outputs and job_names:
            missing.append(f"{output_dir} (no .png or .csv)")

    return (len(missing) == 0, missing)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Batch validation: run Abaqus jobs, then comparisons, then check outputs."
    )
    parser.add_argument(
        "--jobs",
        type=str,
        nargs="*",
        default=None,
        help=f"Job names to run (default: {DEFAULT_JOBS}). Use empty to skip Abaqus run.",
    )
    parser.add_argument(
        "--compare-only",
        action="store_true",
        help="Only run comparison scripts and check outputs; do not run Abaqus.",
    )
    parser.add_argument(
        "--no-regenerate",
        action="store_true",
        help="Skip Abaqus script regeneration if script already exists.",
    )
    args = parser.parse_args()

    jobs = args.jobs if args.jobs is not None else DEFAULT_JOBS
    run_abaqus = not args.compare_only and len(jobs) > 0

    if run_abaqus:
        if not _abaqus_can_run():
            print(
                "Abaqus or abqpy not available; skipping Abaqus run. "
                "Use --compare-only to only run comparisons.",
                file=sys.stderr,
            )
        else:
            if not _run_abaqus_jobs(jobs, regenerate=not args.no_regenerate):
                print("One or more Abaqus jobs failed.", file=sys.stderr)
                # Continue to run comparisons and checks

    print("\n--- Running comparison scripts ---")
    if not _run_comparisons():
        print("Comparison scripts reported failures.", file=sys.stderr)

    print("\n--- Checking output files ---")
    all_found, missing = _check_outputs(jobs)
    if all_found:
        print("All expected output files present.")
    else:
        for path in missing:
            print(f"  Missing: {path}", file=sys.stderr)
        print("Some expected output files are missing.", file=sys.stderr)

    return 0 if all_found else 1


if __name__ == "__main__":
    sys.exit(main())
