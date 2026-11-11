"""
Buckling load-factor verification: FEM vs Euler analytical vs Engesser (Timoshenko) analytical.

Reads primary_results/buckling/buckling_load_factors.csv from result directories under
results/. Loads properties from the pinned benchmark job at
jobs/job_benchmark_linear_buckling_column/.

Outputs:
  - buckling_plots/buckling_load_factor_verification.png
  - buckling_plots/buckling_load_factor_verification.csv
"""
from __future__ import annotations

import csv
import glob
import sys
from pathlib import Path
from typing import Final

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_BENCHMARK_JOB_DIR: Final[Path] = PROJECT_ROOT / "jobs" / "job_benchmark_linear_buckling_column"
_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"


def _read_load_factors(csv_path: Path) -> list[tuple[int, float]]:
    rows: list[tuple[int, float]] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            rows.append((int(row["mode_index"]), float(row["load_factor"])))
    return rows


def run_buckling_verification() -> None:
    from verification.buckling.analytical_buckling import (
        euler_buckling_load_factor,
        load_euler_buckling_benchmark,
    )
    from verification.buckling.timoshenko_analytical_buckling import (
        engesser_buckling_load_factor,
        load_timoshenko_column_benchmark,
    )

    out_dir = SCRIPT_DIR / "buckling_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(_RESULTS_DIR / "*" / "primary_results" / "buckling" / "buckling_load_factors.csv")
    found = sorted(glob.glob(pattern))
    if not found:
        print("No FEM buckling results found. Run the benchmark job first.")
        return

    euler_bm = load_euler_buckling_benchmark(_BENCHMARK_JOB_DIR)
    engesser_bm = load_timoshenko_column_benchmark(_BENCHMARK_JOB_DIR)
    lf_euler = euler_buckling_load_factor(euler_bm)
    lf_engesser = engesser_buckling_load_factor(engesser_bm)

    csv_rows: list[list] = []

    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        result_name = csv_path.parents[2].name
        fem_rows = _read_load_factors(csv_path)
        if not fem_rows:
            print(f"  {result_name}: empty buckling_load_factors.csv, skipping.")
            continue
        n_modes = len(fem_rows)

        mode_labels = [f"Mode {r[0]}" for r in fem_rows]
        fem_vals = np.array([abs(r[1]) for r in fem_rows])

        # Mode-1 relative errors against Euler and Engesser
        lf_mode1 = fem_vals[0]
        err_euler = 100.0 * (lf_mode1 - lf_euler) / lf_euler
        err_engesser = 100.0 * (lf_mode1 - lf_engesser) / lf_engesser

        fig, (ax_lf, ax_err) = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(
            f"Buckling load-factor verification — {result_name}\n"
            f"Benchmark: L={euler_bm.length_m:.1f} m, K={euler_bm.effective_length_factor:.1f}, "
            f"P_ref={euler_bm.reference_axial_load_n:.0f} N",
            fontsize=11,
        )

        x = np.arange(n_modes)
        ax_lf.bar(x, fem_vals, 0.5, label="FEM |λ|", color="#4F81BD")
        ax_lf.axhline(lf_euler, color="#C0504D", linewidth=1.5, ls="--", label=f"Euler λ={lf_euler:.2f}")
        ax_lf.axhline(lf_engesser, color="#9BBB59", linewidth=1.5, ls="-.", label=f"Engesser λ={lf_engesser:.2f}")
        ax_lf.set_xlabel("Mode")
        ax_lf.set_ylabel("Load factor |λ|")
        ax_lf.set_title("Load factor comparison")
        ax_lf.set_xticks(x)
        ax_lf.set_xticklabels(mode_labels)
        ax_lf.legend()
        ax_lf.grid(axis="y", ls="--", alpha=0.6)

        ax_err.bar([-0.2, 0.2], [err_euler, err_engesser], 0.35,
                   color=["#C0504D", "#9BBB59"],
                   label=["FEM vs Euler", "FEM vs Engesser"])
        ax_err.axhline(0, color="k", linewidth=0.8)
        ax_err.axhline(2, color="orange", linewidth=0.8, ls="--", label="±2% tolerance")
        ax_err.axhline(-2, color="orange", linewidth=0.8, ls="--")
        ax_err.set_xlabel("Reference")
        ax_err.set_ylabel("Mode-1 error [%]")
        ax_err.set_title("Mode-1 relative error vs analytical")
        ax_err.set_xticks([-0.2, 0.2])
        ax_err.set_xticklabels(["Euler", "Engesser"])
        ax_err.legend()
        ax_err.grid(axis="y", ls="--", alpha=0.6)

        fig.tight_layout()
        out_png = out_dir / f"buckling_load_factor_verification_{result_name}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"  Saved: {out_png}")

        for i in range(n_modes):
            csv_rows.append([result_name, fem_rows[i][0], fem_vals[i]])
        csv_rows.append([result_name, "analytical_euler", lf_euler])
        csv_rows.append([result_name, "analytical_engesser", lf_engesser])

        print(
            f"  {result_name}: mode-1 FEM={lf_mode1:.4f}, "
            f"Euler={lf_euler:.4f} (err={err_euler:.2f}%), "
            f"Engesser={lf_engesser:.4f} (err={err_engesser:.2f}%)"
        )

    if csv_rows:
        csv_out = out_dir / "buckling_load_factor_verification.csv"
        with open(csv_out, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["result", "mode_index", "load_factor"])
            writer.writerows(csv_rows)
        print(f"  Saved: {csv_out}")

    print("Buckling verification done.")


if __name__ == "__main__":
    run_buckling_verification()
