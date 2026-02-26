#!/usr/bin/env python
"""
Run Abaqus for all validation jobs (every job_XXXX_nN under jobs/).

Use to re-run all Abaqus results in one go. Run from project root:
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py --dry-run
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py --jobs job_0000_n8 job_0005_n16
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from post_processing.validation_visualisers.abaqus.config import JOBS_DIR
from post_processing.validation_visualisers.abaqus.run_abaqus_cae import run_job

# job_XXXX_nN (e.g. job_0000_n8, job_0007_n128)
JOB_DIR_PATTERN = re.compile(r"^job_\d{4}_n\d+$")


def discover_all_jobs() -> list[str]:
    """Return sorted list of job names (job_XXXX_nN) that exist under JOBS_DIR."""
    if not JOBS_DIR.is_dir():
        return []
    names = []
    for p in JOBS_DIR.iterdir():
        if p.is_dir() and JOB_DIR_PATTERN.match(p.name):
            names.append(p.name)
    return sorted(names)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run Abaqus for all (or selected) validation jobs."
    )
    parser.add_argument(
        "--jobs",
        type=str,
        nargs="*",
        default=None,
        help="Job names to run (default: all discovered job_XXXX_nN under jobs/).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print job list and exit; do not run Abaqus.",
    )
    parser.add_argument(
        "--no-regenerate",
        action="store_true",
        help="Skip script regeneration; use existing generated scripts.",
    )
    args = parser.parse_args()

    if args.jobs is not None:
        job_list = sorted(set(args.jobs))
        for j in job_list:
            if not JOB_DIR_PATTERN.match(j):
                print(f"Warning: '{j}' does not match job_XXXX_nN.", file=sys.stderr)
            job_dir = JOBS_DIR / j
            if not job_dir.is_dir():
                print(f"Error: job directory not found: {job_dir}", file=sys.stderr)
                return 1
    else:
        job_list = discover_all_jobs()
        if not job_list:
            print(f"No job directories found under {JOBS_DIR}.", file=sys.stderr)
            return 1

    print(f"Jobs to run ({len(job_list)}): {', '.join(job_list)}")
    if args.dry_run:
        return 0

    regenerate = not args.no_regenerate
    all_ok = True
    for i, job_name in enumerate(job_list, start=1):
        print(f"[{i}/{len(job_list)}] Running Abaqus for {job_name} ...")
        ok = run_job(job_name, regenerate=regenerate)
        print(f"  {'OK' if ok else 'FAILED'}")
        if not ok:
            all_ok = False

    print("Done.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
