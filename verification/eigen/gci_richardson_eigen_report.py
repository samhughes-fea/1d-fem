"""
GCI and Richardson extrapolation report for eigen mode-1 frequency.

Uses three mesh levels n=32, 64, 128 (refinement ratio R=2, safety factor F_s=1.25).
Compares Richardson extrapolate to Euler-Bernoulli analytical formula.

Outputs:
  - eigen_plots/gci_tables/gci_richardson_eigen.csv
  Console summary with p_obs, GCI_12, GCI_23, phi_ext, Roark delta.
"""
from __future__ import annotations

import csv
import glob
import math
import re
import sys
from pathlib import Path
from typing import Final

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_BENCHMARK_JOB_DIR: Final[Path] = PROJECT_ROOT / "jobs" / "job_benchmark_eigen_cantilever"
_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
_N_RE: Final[re.Pattern] = re.compile(r"benchmark_eigen_cantilever_n(\d+)")

N_FINE, N_MED, N_COARSE = 128, 64, 32
R: Final[float] = 2.0
F_S: Final[float] = 1.25


def _read_f1(csv_path: Path) -> float | None:
    try:
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                return float(row["frequency_hz"])
    except Exception:
        return None
    return None


def _get_f1_for_n(n: int, found: list[str]) -> float | None:
    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        m = _N_RE.search(csv_path.parents[2].name)
        if m and int(m.group(1)) == n:
            return _read_f1(csv_path)
    return None


def _richardson_gci(
    phi_fine: float, phi_med: float, phi_coarse: float
) -> tuple[float, float, float, float]:
    e21 = phi_med - phi_fine
    e32 = phi_coarse - phi_med
    if abs(e21) < 1e-30 and abs(e32) < 1e-30:
        return float("nan"), phi_fine, 0.0, 0.0
    ratio = e32 / e21 if abs(e21) >= 1e-30 else float("nan")
    if ratio <= 0 or not math.isfinite(ratio):
        return float("nan"), phi_fine, float("nan"), float("nan")
    p_obs = math.log(abs(ratio)) / math.log(R)
    denom = (R ** p_obs) - 1.0
    if denom <= 0 or not math.isfinite(denom):
        return float("nan"), phi_fine, float("nan"), float("nan")
    phi_ext = phi_fine + (phi_fine - phi_med) / denom
    gci_12 = 100.0 * F_S * abs(e21 / phi_fine) / denom if abs(phi_fine) >= 1e-30 else float("nan")
    gci_23 = 100.0 * F_S * abs(e32 / phi_med) / denom if abs(phi_med) >= 1e-30 else float("nan")
    return p_obs, phi_ext, gci_12, gci_23


def _run_report() -> None:
    from verification.eigen.analytical_eigen import (
        cantilever_eigen_frequencies_hz,
        load_cantilever_eigen_benchmark,
    )

    out_dir = SCRIPT_DIR / "eigen_plots" / "gci_tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    pattern = str(
        _RESULTS_DIR / "*benchmark_eigen_cantilever_n*" / "primary_results" / "eigen" / "eigen_frequencies.csv"
    )
    found = sorted(glob.glob(pattern))

    missing = [n for n in (N_FINE, N_MED, N_COARSE) if _get_f1_for_n(n, found) is None]
    if missing:
        print(
            f"GCI/Richardson skipped: missing results for n={missing}.\n"
            f"Run job_benchmark_eigen_cantilever_n{missing[0]} (and others) first."
        )
        return

    f1_fine = _get_f1_for_n(N_FINE, found)
    f1_med = _get_f1_for_n(N_MED, found)
    f1_coarse = _get_f1_for_n(N_COARSE, found)

    p_obs, phi_ext, gci_12, gci_23 = _richardson_gci(f1_fine, f1_med, f1_coarse)

    benchmark = load_cantilever_eigen_benchmark(_BENCHMARK_JOB_DIR)
    f1_ref = cantilever_eigen_frequencies_hz(benchmark, num_modes=1)[0]

    error_fine_vs_ref = 100.0 * (f1_fine - f1_ref) / f1_ref if f1_ref else float("nan")
    error_ext_vs_ref = 100.0 * (phi_ext - f1_ref) / f1_ref if f1_ref and math.isfinite(phi_ext) else float("nan")
    r_p = (R ** p_obs) if math.isfinite(p_obs) else float("nan")
    gci_ratio = gci_23 / gci_12 if math.isfinite(gci_12) and abs(gci_12) > 1e-30 and math.isfinite(gci_23) else float("nan")

    csv_path = out_dir / "gci_richardson_eigen.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "QoI", "n_fine", "n_med", "n_coarse",
            "phi_fine", "phi_med", "phi_coarse",
            "p_obs", "r_p", "phi_ext", "f1_ref",
            "error_fine_vs_ref_pct", "error_ext_vs_ref_pct",
            "GCI_12_pct", "GCI_23_pct", "GCI_23_12_ratio",
        ])
        writer.writerow([
            "f1_hz", N_FINE, N_MED, N_COARSE,
            f"{f1_fine:.12e}", f"{f1_med:.12e}", f"{f1_coarse:.12e}",
            f"{p_obs:.6f}", f"{r_p:.6f}", f"{phi_ext:.12e}", f"{f1_ref:.12e}",
            f"{error_fine_vs_ref_pct:.6f}" if math.isfinite(error_fine_vs_ref) else "nan",
            f"{error_ext_vs_ref:.6f}" if math.isfinite(error_ext_vs_ref) else "nan",
            f"{gci_12:.6f}" if math.isfinite(gci_12) else "nan",
            f"{gci_23:.6f}" if math.isfinite(gci_23) else "nan",
            f"{gci_ratio:.6f}" if math.isfinite(gci_ratio) else "nan",
        ])
    print(f"Saved: {csv_path}")

    def _fmt(v: float, d: int = 4) -> str:
        return f"{v:.{d}f}" if math.isfinite(v) else "N/A"

    print(f"\n--- GCI / Richardson: eigen mode-1 frequency ---")
    print(f"  Grids: n={N_COARSE}, {N_MED}, {N_FINE}   R={R}   F_s={F_S}")
    print(f"  phi_coarse (n={N_COARSE}) = {f1_coarse:.6f} Hz")
    print(f"  phi_med    (n={N_MED})  = {f1_med:.6f} Hz")
    print(f"  phi_fine   (n={N_FINE}) = {f1_fine:.6f} Hz")
    print(f"  phi_ext    (Richardson) = {phi_ext:.6f} Hz")
    print(f"  f1_ref     (EB analyt.) = {f1_ref:.6f} Hz")
    print(f"  p_obs   = {_fmt(p_obs)}  (theoretical: 4 for cubic EB element)")
    print(f"  r^p     = {_fmt(r_p)}")
    print(f"  GCI_12  = {_fmt(gci_12)}%  (fine-grid GCI)")
    print(f"  GCI_23  = {_fmt(gci_23)}%  (medium-grid GCI)")
    print(f"  GCI_23/GCI_12 = {_fmt(gci_ratio)}  (expect ≈ r^p in asymptotic range)")
    print(f"  Δ(phi_fine, ref) = {_fmt(error_fine_vs_ref)}%")
    print(f"  Δ(phi_ext,  ref) = {_fmt(error_ext_vs_ref)}%")
    print("")


if __name__ == "__main__":
    _run_report()
    sys.exit(0)
