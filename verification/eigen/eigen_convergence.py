"""
Eigen frequency mesh convergence: mode-1 relative error vs n_elements (log-log).

Reads primary_results/eigen/eigen_frequencies.csv from result directories under
results/ that correspond to job_benchmark_eigen_cantilever_n* jobs.
Computes relative error in mode-1 frequency against Euler-Bernoulli analytical formula.

Outputs:
  - eigen_plots/eigen_convergence.png
  - eigen_plots/eigen_convergence.csv
"""
from __future__ import annotations

import csv
import glob
import re
import sys
from pathlib import Path
from typing import Final

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_BENCHMARK_JOB_DIR: Final[Path] = PROJECT_ROOT / "jobs" / "job_benchmark_eigen_cantilever"
_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
_N_RE: Final[re.Pattern] = re.compile(r"benchmark_eigen_cantilever_n(\d+)")


def _read_f1(csv_path: Path) -> float | None:
    try:
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                return float(row["frequency_hz"])
    except Exception:
        return None
    return None


def run_eigen_convergence() -> None:
    from verification.eigen.analytical_eigen import (
        cantilever_eigen_frequencies_hz,
        load_cantilever_eigen_benchmark,
    )

    out_dir = SCRIPT_DIR / "eigen_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(
        _RESULTS_DIR / "*benchmark_eigen_cantilever_n*" / "primary_results" / "eigen" / "eigen_frequencies.csv"
    )
    found = sorted(glob.glob(pattern))
    if not found:
        print("No FEM eigen convergence results found — run benchmark n-variant jobs first.")
        return

    benchmark = load_cantilever_eigen_benchmark(_BENCHMARK_JOB_DIR)
    f1_ref = cantilever_eigen_frequencies_hz(benchmark, num_modes=1)[0]

    # Collect one (n, f1) per mesh size — keep the last (most recent) found per n
    n_to_f1: dict[int, float] = {}
    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        m = _N_RE.search(csv_path.parents[2].name)
        if not m:
            continue
        n = int(m.group(1))
        f1 = _read_f1(csv_path)
        if f1 is not None:
            n_to_f1[n] = f1

    if not n_to_f1:
        print("No valid f₁ values found in eigen convergence results.")
        return

    ns = sorted(n_to_f1.keys())
    f1s = [n_to_f1[n] for n in ns]
    errs = [100.0 * abs(f - f1_ref) / f1_ref for f in f1s]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(ns, errs, "o-", color="#4F81BD", linewidth=1.5, markersize=6, label="FEM vs EB analytical")
    ax.set_xlabel("Number of elements $n$")
    ax.set_ylabel(r"$|f_1^{\rm FEM} - f_1^{\rm ref}| / f_1^{\rm ref}$ [%]")
    ax.set_title(
        f"Eigen mode-1 mesh convergence\n"
        f"EB ref: $f_1$ = {f1_ref:.3f} Hz  "
        f"(L = {benchmark.length_m:.1f} m, "
        f"$I_y$ = {benchmark.second_moment_m4:.3e} m$^4$)"
    )
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    out_png = out_dir / "eigen_convergence.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_png}")

    csv_out = out_dir / "eigen_convergence.csv"
    with open(csv_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["n_elements", "f1_fem_hz", "f1_ref_hz", "abs_rel_error_pct"])
        for n, f1 in zip(ns, f1s):
            err = 100.0 * abs(f1 - f1_ref) / f1_ref
            writer.writerow([n, f"{f1:.10e}", f"{f1_ref:.10e}", f"{err:.6f}"])
    print(f"Saved: {csv_out}")

    print(f"\n--- Eigen f₁ convergence  (EB analytical f₁ = {f1_ref:.4f} Hz) ---")
    for n, f1 in zip(ns, f1s):
        print(f"  n={n:4d}:  f₁ = {f1:.6f} Hz   |error| = {100.0 * abs(f1 - f1_ref) / f1_ref:.4f}%")


if __name__ == "__main__":
    run_eigen_convergence()
