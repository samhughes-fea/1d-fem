"""
GCI and Richardson extrapolation report for NL static tip displacement.

For each of the 6 NL static load cases, uses n=32, 64, 128 (R=2, F_s=1.25).
Loads target δ/L ≈ 20 % (large-deflection regime). No analytical reference is
used here — Richardson extrapolation is applied purely to the FEM results.

Output:
  nonlinear_plots/gci_tables/gci_richardson_nl_static.csv
  Console summary per load case.
"""
from __future__ import annotations

import csv
import glob
import math
import re
import sys
from pathlib import Path
from typing import Final

import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from verification.roark.roark_utilities.roarks_formulas_euler_bernoulli_point import (
    RoarksFormulaePointLoad,
)
from verification.roark.roark_utilities.roarks_formulas_euler_bernoulli_distributed import (
    RoarksFormulaeDistributedLoad,
)

_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
_JOBS_DIR: Final[Path] = PROJECT_ROOT / "jobs"
_OUT_DIR: Final[Path] = SCRIPT_DIR / "nonlinear_plots" / "gci_tables"

_LOAD_CASES: Final[list[dict]] = [
    {"case": "cantilever_tip",  "label": "Tip load",      "roark_type": "point", "roark_kwarg": {"load_type": "end"}},
    {"case": "midspan_point",   "label": "Mid-span load", "roark_type": "point", "roark_kwarg": {"load_type": "mid"}},
    {"case": "quarter_point",   "label": "Quarter load",  "roark_type": "point", "roark_kwarg": {"load_type": "quarter"}},
    {"case": "udl",             "label": "UDL",           "roark_type": "dist",  "roark_kwarg": {"load_type": "udl"}},
    {"case": "triangular",      "label": "Triangular",    "roark_type": "dist",  "roark_kwarg": {"load_type": "triangular"}},
    {"case": "parabolic",       "label": "Parabolic",     "roark_type": "dist",  "roark_kwarg": {"load_type": "parabolic"}},
]

N_FINE, N_MED, N_COARSE = 128, 64, 32
R: Final[float] = 2.0
F_S: Final[float] = 1.25


def _beam_params(base_job: str) -> dict | None:
    job_dir = _JOBS_DIR / base_job
    if not job_dir.is_dir():
        return None
    try:
        grid_lines = (job_dir / "grid.txt").read_text(encoding="utf-8").splitlines()
        xs = []
        for line in grid_lines:
            if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    xs.append(float(parts[1]))
                except ValueError:
                    pass
        L = max(xs) if xs else None

        mat_lines = (job_dir / "material.txt").read_text(encoding="utf-8").splitlines()
        E = None
        for line in mat_lines:
            if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                E = float(parts[1])
                break

        sec_lines = (job_dir / "section.txt").read_text(encoding="utf-8").splitlines()
        I_y = None
        for line in sec_lines:
            if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                I_y = float(parts[3])  # 0=eid, 1=A, 2=I_x, 3=I_y
                break

        if None in (L, E, I_y):
            return None

        load_mag = None
        pt_path = job_dir / "point_load.txt"
        dist_path = job_dir / "distributed_load.txt"
        if pt_path.is_file():
            for line in pt_path.read_text(encoding="utf-8").splitlines():
                if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    load_mag = abs(float(parts[4]))
                    break
        elif dist_path.is_file():
            max_fy = 0.0
            for line in dist_path.read_text(encoding="utf-8").splitlines():
                if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    max_fy = max(max_fy, abs(float(parts[4])))
            load_mag = max_fy if max_fy > 0 else None

        return {"L": L, "E": E, "I": I_y, "P": load_mag} if load_mag else None
    except Exception:
        return None


def _roark_tip(params: dict, cfg: dict) -> float | None:
    try:
        L, E, I, P = params["L"], params["E"], params["I"], params["P"]
        x_eval = np.linspace(0.0, L, 1000)
        if cfg["roark_type"] == "point":
            r = RoarksFormulaePointLoad(L, E, I, P, **cfg["roark_kwarg"])
        else:
            r = RoarksFormulaeDistributedLoad(L, E, I, P, **cfg["roark_kwarg"])
        return abs(float(r.deflection(x_eval)[-1]))
    except Exception:
        return None


def _get_tip_for_n(case: str, n: int) -> float | None:
    pattern = str(
        _RESULTS_DIR
        / f"job_benchmark_nl_static_{case}_n{n}_*"
        / "primary_results"
        / "nonlinear_static_validation"
        / "*_tip_load_history.csv"
    )
    found = sorted(glob.glob(pattern), key=lambda p: Path(p).stat().st_mtime, reverse=True)
    if not found:
        return None
    try:
        rows = []
        with open(found[0], "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                rows.append((float(row["load_factor"]), float(row["tip_displacement"])))
        return abs(rows[-1][1]) if rows else None
    except Exception:
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
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "job_name", "n_fine", "n_med", "n_coarse",
        "phi_fine", "phi_med", "phi_coarse",
        "p_obs", "r_p", "phi_ext", "roark_ref",
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

        phi_fine = _get_tip_for_n(case, N_FINE)
        phi_med = _get_tip_for_n(case, N_MED)
        phi_coarse = _get_tip_for_n(case, N_COARSE)

        missing = [n for n, phi in ((N_FINE, phi_fine), (N_MED, phi_med), (N_COARSE, phi_coarse)) if phi is None]
        if missing:
            print(f"Skip {case}: missing results for n={missing}")
            continue

        params = _beam_params(base_job)
        if params is None:
            print(f"Skip {case}: could not load beam parameters")
            continue
        roark_ref = _roark_tip(params, cfg)
        if roark_ref is None:
            print(f"Skip {case}: Roark computation failed")
            continue

        p_obs, phi_ext, gci_12, gci_23 = _richardson_gci(phi_fine, phi_med, phi_coarse)

        error_fine_pct = 100.0 * (phi_fine - roark_ref) / roark_ref if roark_ref else float("nan")
        error_ext_pct = (
            100.0 * (phi_ext - roark_ref) / roark_ref
            if roark_ref and math.isfinite(phi_ext)
            else float("nan")
        )
        r_p = (R ** p_obs) if math.isfinite(p_obs) else float("nan")
        gci_ratio = (
            gci_23 / gci_12
            if (math.isfinite(gci_12) and abs(gci_12) > 1e-30 and math.isfinite(gci_23))
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
            "roark_ref": f"{roark_ref:.12e}",
            "error_fine_pct": _fmt(error_fine_pct),
            "error_ext_pct": _fmt(error_ext_pct),
            "GCI_12_pct": _fmt(gci_12),
            "GCI_23_pct": _fmt(gci_23),
            "GCI_23_12_ratio": _fmt(gci_ratio),
        })

        any_done = True
        print(f"\n--- GCI / Richardson: {cfg['label']} ({case}) ---")
        print(f"  Grids: n={N_COARSE}, {N_MED}, {N_FINE}   R={R}   F_s={F_S}")
        print(f"  phi_coarse (n={N_COARSE}) = {phi_coarse:.8e} m")
        print(f"  phi_med    (n={N_MED})  = {phi_med:.8e} m")
        print(f"  phi_fine   (n={N_FINE}) = {phi_fine:.8e} m")
        if math.isfinite(phi_ext):
            print(f"  phi_ext    (Richardson) = {phi_ext:.8e} m")
        else:
            print("  phi_ext    (Richardson) = N/A")
        print(f"  roark_ref               = {roark_ref:.8e} m")
        print(f"  p_obs   = {_fmt(p_obs)}  (theoretical: 4 for cubic EB element)")
        print(f"  r^p     = {_fmt(r_p)}")
        print(f"  GCI_12  = {_fmt(gci_12)}%  (fine-grid GCI)")
        print(f"  GCI_23  = {_fmt(gci_23)}%  (medium-grid GCI)")
        print(f"  GCI_23/GCI_12 = {_fmt(gci_ratio)}  (expect ≈ r^p in asymptotic range)")
        print(f"  Δ(phi_fine, Roark) = {_fmt(error_fine_pct)}%")
        print(f"  Δ(phi_ext,  Roark) = {_fmt(error_ext_pct)}%")

    if not any_done:
        print(
            "\nNo NL static GCI results computed. "
            "Run FEM solver for nl_static n=32, 64, and 128 jobs first."
        )
        return

    csv_path = _OUT_DIR / "gci_richardson_nl_static.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"\nSaved: {csv_path}")
    print("NL static GCI/Richardson report done.")


if __name__ == "__main__":
    _run_report()
    sys.exit(0)
