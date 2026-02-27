#!/usr/bin/env python
"""
Run Abaqus for all validation jobs (every job_XXXX_nN under jobs/).
Abaqus results are generated using the Python Abaqus package (abqpy): each job's
CAE script is run with project Python; abqpy launches Abaqus, which runs the
script and writes results to abaqus_results/job_XXXX_nN/.
With --from-results, discover jobs from post_processing/results instead.

Use to re-run all Abaqus results in one go. Run from project root:
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py --dry-run
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py --from-results
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py --from-results --script-only
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py --jobs job_0000_n8 job_0005_n16
  python post_processing/validation_visualisers/run_all_abaqus_jobs.py --n500-only   # Abaqus converged reference batch
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

from post_processing.validation_visualisers.abaqus.config import (
    ABAQUS_REFERENCE_N,
    FEM_RESULTS_DIR,
    JOBS_DIR,
)
from post_processing.validation_visualisers.abaqus.run_abaqus_cae import (
    generate_script_only,
    run_job,
)
from post_processing.validation_visualisers.job_discovery import (
    discover_job_names_from_results,
)

# job_XXXX_nN (e.g. job_0000_n8, job_0007_n128)
JOB_DIR_PATTERN = re.compile(r"^job_\d{4}_n\d+$")

# Validation base jobs 0-11 (used by deformation/section_forces comparisons)
VALIDATION_BASE_IDS = 12


def discover_n500_jobs() -> list[str]:
    """Return job_0000_n500 ... job_0011_n500 that exist under JOBS_DIR (Abaqus converged reference batch)."""
    job_list = [
        f"job_{i:04d}_n{ABAQUS_REFERENCE_N}"
        for i in range(VALIDATION_BASE_IDS)
    ]
    return [j for j in job_list if (JOBS_DIR / j).is_dir()]


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
    parser.add_argument(
        "--from-results",
        action="store_true",
        help="Discover job names from post_processing/results (timestamped dirs); only run for jobs that also exist under jobs/.",
    )
    parser.add_argument(
        "--script-only",
        action="store_true",
        help="Only generate Abaqus CAE scripts; do not run Abaqus (no license required).",
    )
    parser.add_argument(
        "--n500-only",
        action="store_true",
        help="Run only the Abaqus n500 reference batch (job_0000_n500 ... job_0011_n500). Use for validation converged reference.",
    )
    args = parser.parse_args()

    if args.n500_only:
        job_list = discover_n500_jobs()
        if not job_list:
            print(
                f"No job_XXXX_n{ABAQUS_REFERENCE_N} directories found under {JOBS_DIR}. "
                "Run the mesh variant scripts (create_point_load_mesh_variants.py, create_distributed_mesh_variants.py, create_timoshenko_mesh_variants.py) first.",
                file=sys.stderr,
            )
            return 1
        missing = [f"job_{i:04d}_n{ABAQUS_REFERENCE_N}" for i in range(VALIDATION_BASE_IDS) if f"job_{i:04d}_n{ABAQUS_REFERENCE_N}" not in job_list]
        if missing:
            print(f"Note: missing n500 job dirs: {', '.join(missing)}", file=sys.stderr)
    elif args.jobs is not None:
        job_list = sorted(set(args.jobs))
        for j in job_list:
            if not JOB_DIR_PATTERN.match(j):
                print(f"Warning: '{j}' does not match job_XXXX_nN.", file=sys.stderr)
            job_dir = JOBS_DIR / j
            if not job_dir.is_dir():
                print(f"Error: job directory not found: {job_dir}", file=sys.stderr)
                return 1
    elif args.from_results:
        from_results = discover_job_names_from_results(FEM_RESULTS_DIR)
        if not from_results:
            print(
                f"No result directories found under {FEM_RESULTS_DIR}.",
                file=sys.stderr,
            )
            return 1
        job_list = []
        skipped = []
        for j in from_results:
            if (JOBS_DIR / j).is_dir():
                job_list.append(j)
            else:
                skipped.append(j)
        if skipped:
            print(
                f"Skipped (no jobs/ input): {', '.join(skipped)}",
                file=sys.stderr,
            )
        if not job_list:
            print(
                "No jobs from results have a matching directory under jobs/.",
                file=sys.stderr,
            )
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
        if args.script_only:
            print(f"[{i}/{len(job_list)}] Generating script for {job_name} ...")
            ok = generate_script_only(job_name, regenerate=regenerate)
        else:
            print(f"[{i}/{len(job_list)}] Running Abaqus for {job_name} ...")
            ok = run_job(job_name, regenerate=regenerate)
        print(f"  {'OK' if ok else 'FAILED'}")
        if not ok:
            all_ok = False

    print("Done.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
