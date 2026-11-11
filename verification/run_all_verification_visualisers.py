#!/usr/bin/env python
"""
Run all verification visualisers to produce .png (and CSV) outputs.
Uses Agg backend so scripts save and exit without blocking.
Run from project root: python verification/run_all_verification_visualisers.py
"""
from __future__ import annotations

import glob
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _check_results_available() -> int:
    """Return number of U_global.csv files under results."""
    results_dir = PROJECT_ROOT / "results"
    pattern = str(results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
    return len(glob.glob(pattern))


def _run_scripts(scripts: list[str]) -> None:
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


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Run all verification visualisers")
    parser.add_argument(
        "--with-validation",
        action="store_true",
        help="Also run validation visualisers (FEM vs Abaqus comparison)",
    )
    args = parser.parse_args()

    n_results = _check_results_available()
    results_dir = PROJECT_ROOT / "results"
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Results dir:  {results_dir}")
    print(f"U_global.csv files found: {n_results}")
    if n_results == 0:
        print(
            "\nNo primary results found. Verification scripts need completed job runs that write\n"
            "  results/job_*/primary_results/global/U_global.csv\n"
            "Run your simulation pipeline so that primary results exist, then run this script again.\n"
        )
        if not args.with_validation:
            return
        print("Proceeding with validation-only scripts.\n")

    if n_results > 0:
        _run_scripts(
            [
                "roark/deformation_convergence.py",
                "roark/distributed_load_convergence.py",
                "roark/roark_verification.py",
                "roark/roark_section_forces_verification.py",
                "grid_convergence_index/gci_richardson_roark_report.py",
                "roark/shear_deformable_verification.py",
                "elastica/cantilever_elastica.py",
                "nonlinear/nl_static_verification.py",
                "nonlinear/nl_static_convergence.py",
                "nonlinear/gci_richardson_nl_static_report.py",
                "buckling/buckling_verification.py",
                "buckling/buckling_convergence.py",
                "buckling/gci_richardson_buckling_report.py",
                "eigen/eigen_verification.py",
                "eigen/eigen_convergence.py",
                "eigen/gci_richardson_eigen_report.py",
                "harmonic/harmonic_verification.py",
                "harmonic/harmonic_convergence.py",
                "harmonic/gci_richardson_harmonic_report.py",
            ]
        )

    if args.with_validation:
        validation_dir = PROJECT_ROOT / "validation"
        validation_scripts = [
            "deformation/deformation_comparison.py",
            "section_forces/section_forces_comparison.py",
            "grid_convergence_study/gci_richardson_abaqus_report.py",
            "grid_convergence_study/u_global_largest_mesh_review.py",
            "grid_convergence_study/csv_to_latex_table.py",
            "static/run_nonlinear_static_validation.py",
            "transient/run_transient_validation.py",
            "nonlinear_buckling/run_nonlinear_buckling_validation.py",
        ]
        for rel_path in validation_scripts:
            path = validation_dir / rel_path
            if not path.is_file():
                print(f"Skip (not found): validation/{rel_path}")
                continue
            print(f"\n--- validation/{rel_path} ---")
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
