#!/usr/bin/env python
"""
Run all validation visualisers (FEM vs Abaqus comparison).
Uses Agg backend. Run from project root:
  python validation/run_all_validation_visualisers.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))


def main() -> None:
    scripts = [
        "deformation/deformation_comparison.py",
        "section_forces/section_forces_comparison.py",
        "static/run_nonlinear_static_validation.py",
        "transient/run_transient_validation.py",
        "nonlinear_buckling/run_nonlinear_buckling_validation.py",
        "grid_convergence_study/gci_richardson_abaqus_report.py",
        "grid_convergence_study/u_global_largest_mesh_review.py",
        "grid_convergence_study/csv_to_latex_table.py",
        "static/gci_richardson_nl_static_validation.py",
        "eigen/run_eigen_validation.py",
        "eigen/gci_richardson_eigen_validation.py",
        "buckling/run_buckling_validation.py",
        "buckling/gci_richardson_buckling_validation.py",
        "harmonic/run_harmonic_validation.py",
        "harmonic/gci_richardson_harmonic_validation.py",
        "transient/gci_richardson_transient_validation.py",
        "nonlinear_buckling/gci_richardson_nl_buckling_validation.py",
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
    print("\nDone.")


if __name__ == "__main__":
    main()
