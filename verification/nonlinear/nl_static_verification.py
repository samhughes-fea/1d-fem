"""
Analytical verification overlay for nonlinear static benchmarks.

Compares FEM tip displacement (at max load_factor) against elastica (large-deflection)
analytical values for each of the 6 NL static load cases. Loads are chosen so that
δ/L ≈ 20 % (linearised), exercising the NL solver's geometric nonlinearity.

Output:
  nonlinear_plots/nl_static_verification.png
  nonlinear_plots/nl_static_verification.csv
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
import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from verification.elastica.cantilever_elastica import (
    elastica_tip_deflection,
)

_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
_JOBS_DIR: Final[Path] = PROJECT_ROOT / "jobs"
_OUT_DIR: Final[Path] = SCRIPT_DIR / "nonlinear_plots"

# All 6 NL static load cases
_LOAD_CASES: Final[list[dict]] = [
    {"case": "cantilever_tip",  "label": "Tip load",      "elastica_type": "end"},
    {"case": "midspan_point",   "label": "Mid-span load", "elastica_type": "mid"},
    {"case": "quarter_point",   "label": "Quarter load",  "elastica_type": "quarter"},
    {"case": "udl",             "label": "UDL",           "elastica_type": "udl"},
    {"case": "triangular",      "label": "Triangular",    "elastica_type": "triangular"},
    {"case": "parabolic",       "label": "Parabolic",     "elastica_type": "parabolic"},
]

_N_RE: Final[re.Pattern] = re.compile(r"nl_static_\w+_n(\d+)")


def _beam_params(base_job: str) -> dict | None:
    """Parse L, E, I_y, load_magnitude from the base job directory."""
    job_dir = _JOBS_DIR / base_job
    if not job_dir.is_dir():
        return None
    try:
        # Grid: L = max node x
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

        # Material: E from first data row, col 1
        mat_lines = (job_dir / "material.txt").read_text(encoding="utf-8").splitlines()
        E = None
        for line in mat_lines:
            if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                E = float(parts[1])
                break

        # Section: I_z = ∫y²dA (bending about Z from F_y loading)
        # Format: eid, A, I_x, I_y, I_z, J_t → I_z is at index 4
        sec_lines = (job_dir / "section.txt").read_text(encoding="utf-8").splitlines()
        I_z = None
        for line in sec_lines:
            if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 5:
                I_z = float(parts[4])  # 0=eid, 1=A, 2=I_x, 3=I_y, 4=I_z
                break

        if None in (L, E, I_z):
            return None

        # Load magnitude
        load_mag = None
        pt_path = job_dir / "point_load.txt"
        dist_path = job_dir / "distributed_load.txt"
        if pt_path.is_file():
            for line in pt_path.read_text(encoding="utf-8").splitlines():
                if line.lstrip().startswith("[") or line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    load_mag = abs(float(parts[4]))  # F_y column (5th token: x y z Fx Fy)
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

        return {"L": L, "E": E, "I": I_z, "P": load_mag} if load_mag else None
    except Exception:
        return None


def _elastica_tip(params: dict, cfg: dict) -> float | None:
    try:
        return elastica_tip_deflection(
            params["L"], params["E"], params["I"],
            params["P"], cfg["elastica_type"],
        )
    except Exception:
        return None


def _discover_fem_tip(case: str) -> tuple[float | None, int | None]:
    """Return (tip_displacement_at_max_load_factor, n_elements) for the highest-n result."""
    pattern = str(_RESULTS_DIR / f"job_benchmark_nl_static_{case}_n*" / "primary_results" / "nonlinear_static_validation")
    dirs = sorted(glob.glob(pattern), key=lambda p: Path(p).stat().st_mtime, reverse=True)
    best_n = -1
    best_tip = None
    for d in dirs:
        m = _N_RE.search(str(d))
        n = int(m.group(1)) if m else 0
        if n < best_n:
            continue
        # find csv inside this dir
        csv_candidates = list(Path(d).glob(f"*_tip_load_history.csv"))
        if not csv_candidates:
            continue
        try:
            rows = []
            with open(csv_candidates[0], "r", encoding="utf-8", newline="") as f:
                for row in csv.DictReader(f):
                    rows.append((float(row["load_factor"]), float(row["tip_displacement"])))
            if rows and rows[-1][0] >= 0.99:
                tip = abs(rows[-1][1])  # displacement at full load
                best_tip = tip
                best_n = n
        except Exception:
            pass
    return best_tip, best_n if best_n > 0 else None


def run_nl_static_verification() -> None:
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows_out: list[dict] = []
    labels, fem_vals, elastica_vals = [], [], []

    for cfg in _LOAD_CASES:
        case = cfg["case"]
        base_job = f"job_benchmark_nl_static_{case}"
        params = _beam_params(base_job)
        if params is None:
            print(f"Skip {case}: could not load beam parameters from {base_job}")
            continue
        elastica_tip = _elastica_tip(params, cfg)
        if elastica_tip is None:
            print(f"Skip {case}: elastica computation failed")
            continue
        fem_tip, n = _discover_fem_tip(case)
        if fem_tip is None:
            print(f"Skip {case}: no FEM results found")
            continue
        error_m = fem_tip - elastica_tip
        error_pct = 100.0 * error_m / elastica_tip if elastica_tip > 0 else float("nan")
        rows_out.append({
            "job_name": base_job,
            "n_elements": n,
            "fem_tip_m": fem_tip,
            "elastica_tip_m": elastica_tip,
            "error_m": error_m,
            "error_pct": error_pct,
        })
        labels.append(cfg["label"])
        fem_vals.append(fem_tip * 1000.0)   # m → mm
        elastica_vals.append(elastica_tip * 1000.0)

    if not rows_out:
        print("No NL static FEM results found. Run FEM solver for nl_static jobs first.")
        return

    # --- Plot ---
    x = np.arange(len(labels))
    w = 0.35
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("NL static verification: FEM vs elastica tip deflection", fontsize=12, fontweight="bold")

    axes[0].bar(x - w / 2, fem_vals, w, label="FEM (NL solver)", color="steelblue")
    axes[0].bar(x + w / 2, elastica_vals, w, label="Elastica (large-deflection)", color="tomato")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=15, ha="right")
    axes[0].set_ylabel("Tip deflection |u_y| (mm)")
    axes[0].set_title("Tip deflection")
    axes[0].legend()
    axes[0].grid(True, axis="y", alpha=0.3)

    err_pcts = [r["error_pct"] for r in rows_out]
    colors = ["green" if abs(e) < 2.0 else "orange" for e in err_pcts]
    axes[1].bar(x, err_pcts, color=colors)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=15, ha="right")
    axes[1].set_ylabel("Error (%)")
    axes[1].set_title("% error (FEM − elastica) / elastica")
    axes[1].grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    out_png = _OUT_DIR / "nl_static_verification.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_png}")

    # --- CSV ---
    out_csv = _OUT_DIR / "nl_static_verification.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["job_name", "n_elements", "fem_tip_m", "elastica_tip_m", "error_m", "error_pct"])
        writer.writeheader()
        for row in rows_out:
            writer.writerow({
                "job_name": row["job_name"],
                "n_elements": row["n_elements"],
                "fem_tip_m": f"{row['fem_tip_m']:.10e}",
                "elastica_tip_m": f"{row['elastica_tip_m']:.10e}",
                "error_m": f"{row['error_m']:.10e}",
                "error_pct": f"{row['error_pct']:.6f}" if isinstance(row["error_pct"], float) and row["error_pct"] == row["error_pct"] else "nan",
            })
    print(f"Wrote {out_csv}")
    print("NL static verification done.")


if __name__ == "__main__":
    run_nl_static_verification()
    sys.exit(0)
