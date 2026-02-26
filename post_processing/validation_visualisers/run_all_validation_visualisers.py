#!/usr/bin/env python
"""
Run all validation visualisers (FEM vs Abaqus comparison).
Uses Agg backend. Run from project root:
  python post_processing/validation_visualisers/run_all_validation_visualisers.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def main() -> None:
    scripts = [
        "deflection_tables/deformation_comparison.py",
        "deflection_tables/gci_richardson_abaqus_report.py",
        "deflection_tables/u_global_largest_mesh_review.py",
        "section_forces/section_forces_comparison.py",
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
