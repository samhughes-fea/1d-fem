"""
Nonlinear static tip displacement mesh convergence: |error %| vs n_elements (log-log).

For each of the 6 NL static load cases, discovers all n-variant FEM results and computes
the % error in final tip displacement against the elastica (large-deflection) analytical value.
Loads target δ/L ≈ 20 % to exercise the NL solver's geometric nonlinearity.

Output:
  nonlinear_plots/nl_static_convergence.png
  nonlinear_plots/nl_static_convergence.csv
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

        # I_z = ∫y²dA (bending about Z from F_y loading); format: eid A I_x I_y I_z J_t
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


def _read_final_tip(csv_path: Path) -> float | None:
    try:
        rows = []
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                rows.append((float(row["load_factor"]), float(row["tip_displacement"])))
        return abs(rows[-1][1]) if rows and rows[-1][0] >= 0.99 else None
    except Exception:
        return None


def run_nl_static_convergence() -> None:
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows_out: list[dict] = []
    any_found = False

    fig, ax = plt.subplots(figsize=(9, 6))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for idx, cfg in enumerate(_LOAD_CASES):
        case = cfg["case"]
        base_job = f"job_benchmark_nl_static_{case}"
        params = _beam_params(base_job)
        if params is None:
            print(f"Skip {case}: could not load beam parameters")
            continue
        elastica_ref = _elastica_tip(params, cfg)
        if elastica_ref is None:
            print(f"Skip {case}: elastica computation failed")
            continue

        pattern = str(
            _RESULTS_DIR
            / f"job_benchmark_nl_static_{case}_n*"
            / "primary_results"
            / "nonlinear_static_validation"
            / "*_tip_load_history.csv"
        )
        found = sorted(glob.glob(pattern))
        n_to_tip: dict[int, float] = {}
        for csv_path_str in found:
            csv_path = Path(csv_path_str)
            m = _N_RE.search(csv_path.parents[2].name)
            if not m:
                continue
            n = int(m.group(1))
            tip = _read_final_tip(csv_path)
            if tip is not None:
                if n not in n_to_tip or tip > 0:
                    n_to_tip[n] = tip

        if not n_to_tip:
            print(f"Skip {case}: no FEM results found")
            continue

        any_found = True
        ns = sorted(n_to_tip.keys())
        errs = [100.0 * abs(n_to_tip[n] - elastica_ref) / elastica_ref for n in ns]
        color = colors[idx % len(colors)]
        ax.loglog(ns, errs, "o-", color=color, linewidth=1.5, markersize=5, label=cfg["label"])

        for n, err in zip(ns, errs):
            rows_out.append({
                "job_name": f"{base_job}_n{n}",
                "n_elements": n,
                "fem_tip_m": f"{n_to_tip[n]:.10e}",
                "elastica_tip_m": f"{elastica_ref:.10e}",
                "error_pct": f"{err:.6f}",
            })

    if not any_found:
        print("No NL static FEM results found. Run FEM solver for nl_static jobs first.")
        plt.close(fig)
        return

    ax.set_xlabel("Number of elements $n$")
    ax.set_ylabel(r"$|u_y^{\rm FEM} - u_y^{\rm elastica}| / u_y^{\rm elastica}$ [%]")
    ax.set_title("NL static tip displacement mesh convergence (vs elastica large-deflection)")
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend(fontsize=9)
    fig.tight_layout()

    out_png = _OUT_DIR / "nl_static_convergence.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")

    out_csv = _OUT_DIR / "nl_static_convergence.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["job_name", "n_elements", "fem_tip_m", "elastica_tip_m", "error_pct"])
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"Wrote {out_csv}")
    print("NL static convergence done.")


if __name__ == "__main__":
    run_nl_static_convergence()
    sys.exit(0)
