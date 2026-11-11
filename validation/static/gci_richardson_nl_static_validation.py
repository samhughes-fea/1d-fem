"""
GCI and Richardson extrapolation validation report for NL static tip displacement.

For each of the 6 NL static load cases, uses FEM n=32/64/128 vs Abaqus n=500 reference.
QoI: tip displacement at max load factor (final row of tip_load_history.csv).

Gracefully skips per-case if FEM or Abaqus reference results are absent.

Output:
  static/gci_tables/gci_richardson_nl_static_validation.csv
  Console summary per load case.
"""
from __future__ import annotations

import csv
import glob
import math
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

_LOAD_CASES: Final[list[dict]] = [
    {"case": "cantilever_tip",  "label": "Tip load"},
    {"case": "midspan_point",   "label": "Mid-span load"},
    {"case": "quarter_point",   "label": "Quarter load"},
    {"case": "udl",             "label": "UDL"},
    {"case": "triangular",      "label": "Triangular"},
    {"case": "parabolic",       "label": "Parabolic"},
]

_OUT_DIR: Final[Path] = SCRIPT_DIR / "gci_tables"

N_FINE, N_MED, N_COARSE = 128, 64, 32
R: Final[float] = 2.0
F_S: Final[float] = 1.25


def _read_final_tip(csv_path: Path) -> float | None:
    try:
        rows = []
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                rows.append((float(row["load_factor"]), float(row["tip_displacement"])))
        return abs(rows[-1][1]) if rows else None
    except Exception:
        return None


def _get_fem_tip_for_n(case: str, n: int) -> float | None:
    pattern = str(
        FEM_RESULTS_DIR
        / f"job_benchmark_nl_static_{case}_n{n}_*"
        / "primary_results"
        / "nonlinear_static_validation"
        / "*_tip_load_history.csv"
    )
    found = sorted(glob.glob(pattern), key=lambda p: Path(p).stat().st_mtime, reverse=True)
    return _read_final_tip(Path(found[0])) if found else None


def _get_abaqus_tip(case: str) -> float | None:
    ref_csv = (
        ABAQUS_RESULTS_DIR
        / f"job_benchmark_nl_static_{case}_n{ABAQUS_REFERENCE_N}"
        / "tip_load_history.csv"
    )
    if not ref_csv.is_file():
        return None
    return _read_final_tip(ref_csv)


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

    fieldnames = [
        "job_name", "n_fine", "n_med", "n_coarse",
        "phi_fine", "phi_med", "phi_coarse",
        "p_obs", "r_p", "phi_ext", "abaqus_ref",
        "error_fine_pct", "error_ext_pct",
        "GCI_12_pct", "GCI_23_pct", "GCI_23_12_ratio",
    ]

    rows_out: list[dict] = []
    any_done = False

    def _fmt(v: float, d: int = 6) -> str:
        return f"{v:.{d}f}" if math.isfinite(v) else "nan"

    for cfg in _LOAD_CASES:
        case = cfg["case"]
        base_job = f"job_benchmark_nl_static_{case}"

        phi_fine = _get_fem_tip_for_n(case, N_FINE)
        phi_med = _get_fem_tip_for_n(case, N_MED)
        phi_coarse = _get_fem_tip_for_n(case, N_COARSE)

        missing_fem = [n for n, phi in ((N_FINE, phi_fine), (N_MED, phi_med), (N_COARSE, phi_coarse)) if phi is None]
        if missing_fem:
            print(f"Skip {case}: missing FEM results for n={missing_fem}")
            continue

        abaqus_ref = _get_abaqus_tip(case)
        if abaqus_ref is None:
            print(
                f"Skip {case}: no Abaqus reference at "
                f"abaqus_results/{base_job}_n{ABAQUS_REFERENCE_N}/tip_load_history.csv"
            )
            continue

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

        rows_out.append({
            "job_name": base_job,
            "n_fine": N_FINE,
            "n_med": N_MED,
            "n_coarse": N_COARSE,
            "phi_fine": f"{phi_fine:.12e}",
            "phi_med": f"{phi_med:.12e}",
            "phi_coarse": f"{phi_coarse:.12e}",
            "p_obs": _fmt(p_obs),
            "r_p": _fmt(r_p),
            "phi_ext": f"{phi_ext:.12e}" if math.isfinite(phi_ext) else "nan",
            "abaqus_ref": f"{abaqus_ref:.12e}",
            "error_fine_pct": _fmt(error_fine_pct),
            "error_ext_pct": _fmt(error_ext_pct),
            "GCI_12_pct": _fmt(gci_12),
            "GCI_23_pct": _fmt(gci_23),
            "GCI_23_12_ratio": _fmt(gci_ratio),
        })

        any_done = True
        print(f"\n--- GCI / Richardson: {cfg['label']} ({case}) vs Abaqus n={ABAQUS_REFERENCE_N} ---")
        print(f"  Grids: n={N_COARSE}, {N_MED}, {N_FINE}   R={R}   F_s={F_S}")
        print(f"  phi_coarse (n={N_COARSE}) = {phi_coarse:.8e} m")
        print(f"  phi_med    (n={N_MED})  = {phi_med:.8e} m")
        print(f"  phi_fine   (n={N_FINE}) = {phi_fine:.8e} m")
        if math.isfinite(phi_ext):
            print(f"  phi_ext    (Richardson) = {phi_ext:.8e} m")
        else:
            print("  phi_ext    (Richardson) = N/A")
        print(f"  abaqus_ref (n={ABAQUS_REFERENCE_N})  = {abaqus_ref:.8e} m")
        print(f"  p_obs   = {_fmt(p_obs)}")
        print(f"  r^p     = {_fmt(r_p)}")
        print(f"  GCI_12  = {_fmt(gci_12)}%")
        print(f"  GCI_23  = {_fmt(gci_23)}%")
        print(f"  GCI_23/GCI_12 = {_fmt(gci_ratio)}")
        print(f"  err(phi_fine, Abaqus) = {_fmt(error_fine_pct)}%")
        print(f"  err(phi_ext,  Abaqus) = {_fmt(error_ext_pct)}%")

    if not any_done:
        print(
            "\nNL static GCI/Richardson validation skipped for all cases. "
            f"Need FEM at n=32/64/128 and Abaqus n={ABAQUS_REFERENCE_N} results."
        )
        return

    csv_path = _OUT_DIR / "gci_richardson_nl_static_validation.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"\nSaved: {csv_path}")
    print("NL static GCI/Richardson validation done.")


if __name__ == "__main__":
    _run_report()
    sys.exit(0)
