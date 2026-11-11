"""
GCI and Richardson extrapolation validation report for buckling mode-1 load factor.

Uses FEM n=32/64/128 vs Abaqus n=500 reference.
QoI: mode-1 load factor (absolute value).

Gracefully skips if FEM or Abaqus reference results are absent.

Output:
  buckling/plots/gci_tables/gci_richardson_buckling_validation.csv
  Console summary.
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
VALIDATION_DIR: Final[Path] = SCRIPT_DIR.parent
PROJECT_ROOT: Final[Path] = VALIDATION_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from validation.abaqus.config import (
    ABAQUS_RESULTS_DIR,
    FEM_RESULTS_DIR,
    ABAQUS_REFERENCE_N,
)

_REF_JOB: Final[str] = f"job_benchmark_linear_buckling_column_n{ABAQUS_REFERENCE_N}"
_OUT_DIR: Final[Path] = SCRIPT_DIR / "plots" / "gci_tables"

N_FINE, N_MED, N_COARSE = 128, 64, 32
R: Final[float] = 2.0
F_S: Final[float] = 1.25


def _read_lambda1(csv_path: Path) -> float | None:
    try:
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                return abs(float(row["load_factor"]))
    except Exception:
        return None
    return None


def _get_lf_for_n(n: int) -> float | None:
    pattern = str(
        FEM_RESULTS_DIR
        / f"job_benchmark_linear_buckling_column_n{n}_*"
        / "primary_results"
        / "buckling"
        / "buckling_load_factors.csv"
    )
    found = sorted(glob.glob(pattern), key=lambda p: Path(p).stat().st_mtime, reverse=True)
    return _read_lambda1(Path(found[0])) if found else None


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
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    phi_fine = _get_lf_for_n(N_FINE)
    phi_med = _get_lf_for_n(N_MED)
    phi_coarse = _get_lf_for_n(N_COARSE)

    missing = [n for n, phi in ((N_FINE, phi_fine), (N_MED, phi_med), (N_COARSE, phi_coarse)) if phi is None]
    if missing:
        print(
            f"Buckling GCI/Richardson validation skipped: missing FEM results for n={missing}.\n"
            f"Run job_benchmark_linear_buckling_column_n{missing[0]} (and others) first."
        )
        return

    ref_csv = ABAQUS_RESULTS_DIR / _REF_JOB / "buckling_load_factors.csv"
    if not ref_csv.is_file():
        print(
            f"Buckling GCI/Richardson validation skipped: no Abaqus reference at {ref_csv}.\n"
            f"Stage ODB-extracted buckling_load_factors.csv into abaqus_results/{_REF_JOB}/ first."
        )
        return
    abaqus_ref = _read_lambda1(ref_csv)
    if abaqus_ref is None:
        print(f"Buckling GCI/Richardson validation skipped: could not read mode-1 load factor from {ref_csv}.")
        return

    p_obs, phi_ext, gci_12, gci_23 = _richardson_gci(phi_fine, phi_med, phi_coarse)
    r_p = (R ** p_obs) if math.isfinite(p_obs) else float("nan")
    error_fine_pct = 100.0 * (phi_fine - abaqus_ref) / abaqus_ref if abaqus_ref else float("nan")
    error_ext_pct = (
        100.0 * (phi_ext - abaqus_ref) / abaqus_ref
        if abaqus_ref and math.isfinite(phi_ext)
        else float("nan")
    )
    gci_ratio = (
        gci_23 / gci_12
        if math.isfinite(gci_12) and abs(gci_12) > 1e-30 and math.isfinite(gci_23)
        else float("nan")
    )

    def _fmt(v: float, d: int = 6) -> str:
        return f"{v:.{d}f}" if math.isfinite(v) else "nan"

    csv_path = _OUT_DIR / "gci_richardson_buckling_validation.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "QoI", "n_fine", "n_med", "n_coarse",
            "phi_fine", "phi_med", "phi_coarse",
            "p_obs", "r_p", "phi_ext", "abaqus_ref",
            "error_fine_pct", "error_ext_pct",
            "GCI_12_pct", "GCI_23_pct", "GCI_23_12_ratio",
        ])
        writer.writerow([
            "mode1_load_factor", N_FINE, N_MED, N_COARSE,
            f"{phi_fine:.12e}", f"{phi_med:.12e}", f"{phi_coarse:.12e}",
            _fmt(p_obs), _fmt(r_p),
            f"{phi_ext:.12e}" if math.isfinite(phi_ext) else "nan",
            f"{abaqus_ref:.12e}",
            _fmt(error_fine_pct), _fmt(error_ext_pct),
            _fmt(gci_12), _fmt(gci_23), _fmt(gci_ratio),
        ])
    print(f"Saved: {csv_path}")

    print(f"\n--- GCI / Richardson: buckling mode-1 load factor (vs Abaqus n={ABAQUS_REFERENCE_N}) ---")
    print(f"  Grids: n={N_COARSE}, {N_MED}, {N_FINE}   R={R}   F_s={F_S}")
    print(f"  phi_coarse (n={N_COARSE}) = {phi_coarse:.6f}")
    print(f"  phi_med    (n={N_MED})  = {phi_med:.6f}")
    print(f"  phi_fine   (n={N_FINE}) = {phi_fine:.6f}")
    print(f"  phi_ext    (Richardson) = {_fmt(phi_ext, 6)}")
    print(f"  abaqus_ref (n={ABAQUS_REFERENCE_N})  = {abaqus_ref:.6f}")
    print(f"  p_obs   = {_fmt(p_obs)}")
    print(f"  r^p     = {_fmt(r_p)}")
    print(f"  GCI_12  = {_fmt(gci_12)}%")
    print(f"  GCI_23  = {_fmt(gci_23)}%")
    print(f"  GCI_23/GCI_12 = {_fmt(gci_ratio)}")
    print(f"  Δ(phi_fine, Abaqus) = {_fmt(error_fine_pct)}%")
    print(f"  Δ(phi_ext,  Abaqus) = {_fmt(error_ext_pct)}%")
    print("")


if __name__ == "__main__":
    _run_report()
    sys.exit(0)
