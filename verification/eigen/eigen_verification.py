"""
Eigen frequency verification: FEM vs Euler-Bernoulli analytical vs Timoshenko analytical.

Reads primary_results/eigen/eigen_frequencies.csv from result directories under
results/. Loads properties from the pinned benchmark job at
jobs/job_benchmark_eigen_cantilever/.

Outputs:
  - eigen_plots/eigen_frequency_verification.png  (grouped bar chart per mode)
  - eigen_plots/eigen_frequency_verification.csv
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

_BENCHMARK_JOB_DIR: Final[Path] = PROJECT_ROOT / "jobs" / "job_benchmark_eigen_cantilever"
_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"


def _read_eigen_frequencies(csv_path: Path) -> list[tuple[int, float]]:
    rows: list[tuple[int, float]] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            rows.append((int(row["mode_index"]), float(row["frequency_hz"])))
    return rows


def run_eigen_verification() -> None:
    from verification.eigen.analytical_eigen import (
        cantilever_eigen_frequencies_hz,
        load_cantilever_eigen_benchmark,
    )
    from verification.eigen.timoshenko_analytical_eigen import (
        load_timoshenko_cantilever_benchmark,
        timoshenko_eigen_frequencies_hz,
    )

    out_dir = SCRIPT_DIR / "eigen_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(_RESULTS_DIR / "*" / "primary_results" / "eigen" / "eigen_frequencies.csv")
    found = sorted(glob.glob(pattern))
    if not found:
        print("No FEM eigen results found. Run the benchmark job first.")
        return

    eb_benchmark = load_cantilever_eigen_benchmark(_BENCHMARK_JOB_DIR)
    t_benchmark = load_timoshenko_cantilever_benchmark(_BENCHMARK_JOB_DIR)

    csv_rows: list[list] = []

    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        result_name = csv_path.parents[2].name
        fem_rows = _read_eigen_frequencies(csv_path)
        if not fem_rows:
            print(f"  {result_name}: empty eigen_frequencies.csv, skipping.")
            continue
        n_modes = len(fem_rows)
        eb_freqs = cantilever_eigen_frequencies_hz(eb_benchmark, num_modes=n_modes)
        t_freqs = timoshenko_eigen_frequencies_hz(t_benchmark, num_modes=n_modes)

        mode_labels = [f"Mode {r[0]}" for r in fem_rows]
        fem_vals = np.array([r[1] for r in fem_rows])
        eb_vals = np.array(eb_freqs)
        t_vals = np.array(t_freqs)

        eb_errors = 100.0 * (fem_vals - eb_vals) / eb_vals
        t_errors = 100.0 * (fem_vals - t_vals) / t_vals

        x = np.arange(n_modes)
        width = 0.25

        fig, (ax_freq, ax_err) = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(
            f"Eigen frequency verification — {result_name}\n"
            f"Benchmark: L={eb_benchmark.length_m:.1f} m, "
            f"E={eb_benchmark.youngs_modulus_pa:.2e} Pa, "
            f"I_y={eb_benchmark.second_moment_m4:.3e} m⁴",
            fontsize=11,
        )

        bars_fem = ax_freq.bar(x - width, fem_vals, width, label="FEM", color="#4F81BD")
        bars_eb = ax_freq.bar(x, eb_vals, width, label="EB analytical", color="#C0504D")
        bars_t = ax_freq.bar(x + width, t_vals, width, label="Timoshenko", color="#9BBB59")
        ax_freq.set_xlabel("Mode")
        ax_freq.set_ylabel("Frequency [Hz]")
        ax_freq.set_title("Frequency comparison")
        ax_freq.set_xticks(x)
        ax_freq.set_xticklabels(mode_labels)
        ax_freq.legend()
        ax_freq.grid(axis="y", ls="--", alpha=0.6)

        ax_err.bar(x - width / 2, eb_errors, width, label="FEM vs EB", color="#C0504D", alpha=0.8)
        ax_err.bar(x + width / 2, t_errors, width, label="FEM vs Timoshenko", color="#9BBB59", alpha=0.8)
        ax_err.axhline(0, color="k", linewidth=0.8)
        ax_err.axhline(2, color="orange", linewidth=0.8, ls="--", label="±2% tolerance")
        ax_err.axhline(-2, color="orange", linewidth=0.8, ls="--")
        ax_err.set_xlabel("Mode")
        ax_err.set_ylabel("Error [%]")
        ax_err.set_title("Relative error vs analytical")
        ax_err.set_xticks(x)
        ax_err.set_xticklabels(mode_labels)
        ax_err.legend()
        ax_err.grid(axis="y", ls="--", alpha=0.6)

        fig.tight_layout()
        out_png = out_dir / f"eigen_frequency_verification_{result_name}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"  Saved: {out_png}")

        for i in range(n_modes):
            csv_rows.append([
                result_name,
                fem_rows[i][0],
                fem_vals[i],
                eb_vals[i],
                t_vals[i],
                eb_errors[i],
                t_errors[i],
            ])
        print(
            f"  {result_name}: modes={n_modes}, "
            f"max |err vs EB|={np.max(np.abs(eb_errors)):.2f}%, "
            f"max |err vs Timoshenko|={np.max(np.abs(t_errors)):.2f}%"
        )

    if csv_rows:
        csv_out = out_dir / "eigen_frequency_verification.csv"
        header = "result,mode_index,freq_fem_hz,freq_eb_hz,freq_timoshenko_hz,error_vs_eb_pct,error_vs_timoshenko_pct"
        with open(csv_out, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(header.split(","))
            writer.writerows(csv_rows)
        print(f"  Saved: {csv_out}")

    print("Eigen verification done.")


if __name__ == "__main__":
    run_eigen_verification()
