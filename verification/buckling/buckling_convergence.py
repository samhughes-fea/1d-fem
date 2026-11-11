"""
Buckling load-factor mesh convergence: mode-1 relative error vs n_elements (log-log).

Reads primary_results/buckling/buckling_load_factors.csv from result directories under
results/ that correspond to job_benchmark_linear_buckling_column_n* jobs.
Computes relative error in mode-1 load factor against the Euler analytical value.

Outputs:
  - buckling_plots/buckling_convergence.png
  - buckling_plots/buckling_convergence.csv
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

_BENCHMARK_JOB_DIR: Final[Path] = PROJECT_ROOT / "jobs" / "job_benchmark_linear_buckling_column"
_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
_N_RE: Final[re.Pattern] = re.compile(r"benchmark_linear_buckling_column_n(\d+)")


def _read_lambda1(csv_path: Path) -> float | None:
    try:
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                return abs(float(row["load_factor"]))
    except Exception:
        return None
    return None


def run_buckling_convergence() -> None:
    from verification.buckling.analytical_buckling import (
        euler_buckling_load_factor,
        load_euler_buckling_benchmark,
    )

    out_dir = SCRIPT_DIR / "buckling_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(
        _RESULTS_DIR
        / "*benchmark_linear_buckling_column_n*"
        / "primary_results"
        / "buckling"
        / "buckling_load_factors.csv"
    )
    found = sorted(glob.glob(pattern))
    if not found:
        print("No FEM buckling convergence results found — run benchmark n-variant jobs first.")
        return

    benchmark = load_euler_buckling_benchmark(_BENCHMARK_JOB_DIR)
    lf_ref = euler_buckling_load_factor(benchmark)

    n_to_lf: dict[int, float] = {}
    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        m = _N_RE.search(csv_path.parents[2].name)
        if not m:
            continue
        n = int(m.group(1))
        lf = _read_lambda1(csv_path)
        if lf is not None:
            n_to_lf[n] = lf

    if not n_to_lf:
        print("No valid λ₁ values found in buckling convergence results.")
        return

    ns = sorted(n_to_lf.keys())
    lfs = [n_to_lf[n] for n in ns]
    errs = [100.0 * abs(lf - lf_ref) / lf_ref for lf in lfs]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(ns, errs, "s-", color="#C0504D", linewidth=1.5, markersize=6, label="FEM vs Euler analytical")
    ax.set_xlabel("Number of elements $n$")
    ax.set_ylabel(r"$|\lambda_1^{\rm FEM} - \lambda_1^{\rm Euler}| / \lambda_1^{\rm Euler}$ [%]")
    ax.set_title(
        f"Buckling mode-1 mesh convergence\n"
        f"Euler ref: $\\lambda_1$ = {lf_ref:.3f}  "
        f"(L = {benchmark.length_m:.1f} m, "
        f"$I$ = {benchmark.second_moment_m4:.3e} m$^4$)"
    )
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    out_png = out_dir / "buckling_convergence.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_png}")

    csv_out = out_dir / "buckling_convergence.csv"
    with open(csv_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["n_elements", "lambda1_fem", "lambda1_euler", "abs_rel_error_pct"])
        for n, lf in zip(ns, lfs):
            err = 100.0 * abs(lf - lf_ref) / lf_ref
            writer.writerow([n, f"{lf:.10e}", f"{lf_ref:.10e}", f"{err:.6f}"])
    print(f"Saved: {csv_out}")

    print(f"\n--- Buckling λ₁ convergence  (Euler analytical λ₁ = {lf_ref:.4f}) ---")
    for n, lf in zip(ns, lfs):
        print(f"  n={n:4d}:  λ₁ = {lf:.6f}   |error| = {100.0 * abs(lf - lf_ref) / lf_ref:.4f}%")


if __name__ == "__main__":
    run_buckling_convergence()
