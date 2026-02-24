#!/usr/bin/env python
"""
Run all verification visualisers to produce .png (and CSV) outputs.
Uses Agg backend so scripts save and exit without blocking.
Run from project root: python post_processing/verification_visualisers/run_all_verification_visualisers.py
"""
from __future__ import annotations

import glob
import subprocess
import sys
from pathlib import Path

# Project root (parent of post_processing)
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _check_results_available() -> int:
    """Return number of U_global.csv files under post_processing/results."""
    results_dir = PROJECT_ROOT / "post_processing" / "results"
    pattern = str(results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
    files = glob.glob(pattern)
    return len(files)


def main() -> None:
    # Pre-check: scripts need primary_results/global/U_global.csv from completed runs
    n_results = _check_results_available()
    results_dir = PROJECT_ROOT / "post_processing" / "results"
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Results dir:  {results_dir}")
    print(f"U_global.csv files found: {n_results}")
    if n_results == 0:
        print(
            "\nNo primary results found. Verification scripts need completed job runs that write\n"
            "  post_processing/results/job_*/primary_results/global/U_global.csv\n"
            "Run your simulation pipeline (e.g. run jobs for job_0000, job_0001, ...) so that\n"
            "primary results exist, then run this script again. No .png files will be generated\n"
            "until then.\n"
        )
        return

    scripts = [
        "deflection_tables/deformation_convergence.py",
        "deflection_tables/distributed_load_convergence.py",
        "deflection_tables/gci_richardson_roark_report.py",
        "shear_deformable_verification.py",
        "roarks_formulas/roark_verification.py",
        "roarks_formulas/roark_section_forces_verification.py",
    ]
    for rel_path in scripts:
        path = SCRIPT_DIR / rel_path
        if not path.is_file():
            print(f"Skip (not found): {rel_path}")
            continue
        print(f"\n--- {rel_path} ---")
        try:
            result = subprocess.run(
                [sys.executable, str(path)],
                cwd=str(PROJECT_ROOT),
                capture_output=False,
                timeout=120,
            )
            if result.returncode != 0:
                print(f"Exit code: {result.returncode}")
        except subprocess.TimeoutExpired:
            print("Timed out after 120s")
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
    print("\nDone.")


if __name__ == "__main__":
    main()
