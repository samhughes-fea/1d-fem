#!/usr/bin/env python3
"""
GCI (Grid Convergence Index) and Richardson extrapolation report for Abaqus vs FEM
tip deflection and tip rotation.

Same structure as verification gci_richardson_roark_deflection_rotation.csv but uses
Abaqus results at the fine mesh (n=128) as reference instead of Roark analytical.
Uses three FEM mesh levels (n=32, 64, 128) for jobs 0,1,2 (point load) and 5,6,7 (distributed).
Requires: FEM result dirs job_XXXX_n32, n64, n128; Abaqus result dirs job_XXXX_n128.
Output: validation_visualisers/output/gci_richardson_abaqus_deflection_rotation.csv

Run from project root:
  python post_processing/validation_visualisers/deflection_tables/gci_richardson_abaqus_report.py
"""
from __future__ import annotations

import glob
import math
import re
import sys
from pathlib import Path
from typing import Final, Optional

import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
VALIDATION_DIR: Final[Path] = SCRIPT_DIR.parent
PROJECT_ROOT: Final[Path] = VALIDATION_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.parsing.grid_parser import GridParser  # type: ignore

FEM_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
JOBS_DIR: Final[Path] = PROJECT_ROOT / "jobs"
ABAQUS_RESULTS_DIR: Final[Path] = VALIDATION_DIR / "abaqus_results"
OUT_DIR: Final[Path] = VALIDATION_DIR / "output"

# Same grid levels as Roark report
N_FINE, N_MED, N_COARSE = 128, 64, 32
R: Final[float] = 2.0
F_S: Final[float] = 1.25

POINT_JOBS: Final[list[tuple[int, str, str]]] = [
    (0, "End", "end"),
    (1, "Mid-span", "mid"),
    (2, "Quarter-span", "quarter"),
]
DIST_JOBS: Final[list[tuple[int, str, str]]] = [
    (5, "UDL", "udl"),
    (6, "Triangular", "triangular"),
    (7, "Parabolic", "parabolic"),
]


def _get_node_coordinates(grid_obj: object) -> np.ndarray:
    if isinstance(grid_obj, dict) and "grid_dictionary" in grid_obj:
        inner = grid_obj["grid_dictionary"]
        if isinstance(inner, dict) and "coordinates" in inner:
            return inner["coordinates"]
    if isinstance(grid_obj, dict) and "node_coordinates" in grid_obj:
        return grid_obj["node_coordinates"]
    if hasattr(grid_obj, "node_coordinates"):
        return getattr(grid_obj, "node_coordinates")
    raise KeyError("grid data does not contain coordinates")


def _read_U_global(file: Path) -> Optional[np.ndarray]:
    """Read U_global.csv; return (n_nodes, 6) array (UX,UY,UZ,RX,RY,RZ)."""
    try:
        data = np.genfromtxt(file, delimiter=",", skip_header=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 2:
            return None
        vals = data[:, 1].astype(float)
        if vals.size % 6 != 0:
            return None
        return vals.reshape(-1, 6)
    except Exception:
        return None


def _get_tip_values(csv_file: Path, jobs_dir: Path, job_input_name: str) -> Optional[tuple[float, float, float]]:
    """Return (L, tip_u_y_mm, tip_theta_z_deg) from U_global CSV and grid. Works for FEM or Abaqus CSV.

    Expects U_global.csv in SI: displacement (UX,UY,UZ) in metres, rotation (RX,RY,RZ) in radians.
    Converts to tip_u_y_mm (mm) and tip_theta_z_deg (degrees) for the report.
    """
    grid_file = jobs_dir / job_input_name / "grid.txt"
    if not grid_file.is_file():
        return None
    U = _read_U_global(csv_file)
    if U is None:
        return None
    grid = GridParser(str(grid_file), str(jobs_dir / job_input_name)).parse()
    try:
        node_coords = _get_node_coordinates(grid)
    except Exception:
        return None
    x = node_coords[:, 0]
    if x.shape[0] != U.shape[0]:
        return None
    L = float(np.max(x))
    tip_idx = int(np.argmax(x))
    tip_uy_mm = float(U[tip_idx, 1] * 1000)
    tip_theta_deg = float(np.degrees(U[tip_idx, 5]))
    return (L, tip_uy_mm, tip_theta_deg)


def _richardson_gci(
    phi_fine: float, phi_med: float, phi_coarse: float
) -> tuple[float, float, float, float]:
    """
    Three-grid Richardson and GCI. Returns (p_obs, phi_ext, gci_12_pct, gci_23_pct).
    """
    e21 = phi_med - phi_fine
    e32 = phi_coarse - phi_med
    if abs(e21) < 1e-30 and abs(e32) < 1e-30:
        return (float("nan"), phi_fine, 0.0, 0.0)
    ratio = e32 / e21 if abs(e21) >= 1e-30 else float("nan")
    if ratio <= 0 or not math.isfinite(ratio):
        p_obs = float("nan")
        phi_ext = phi_fine
        gci_12_pct = float("nan")
        gci_23_pct = float("nan")
    else:
        p_obs = math.log(abs(ratio)) / math.log(R)
        denom = (R ** p_obs) - 1.0
        if denom <= 0 or not math.isfinite(denom):
            phi_ext = phi_fine
            gci_12_pct = float("nan")
            gci_23_pct = float("nan")
        else:
            phi_ext = phi_fine + (phi_fine - phi_med) / denom
            gci_12_pct = 100.0 * F_S * abs(e21 / phi_fine) / denom if abs(phi_fine) >= 1e-30 else float("nan")
            gci_23_pct = 100.0 * F_S * abs(e32 / phi_med) / denom if abs(phi_med) >= 1e-30 else float("nan")
    return (p_obs, phi_ext, gci_12_pct, gci_23_pct)


def _run_report() -> None:
    out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover FEM results: post_processing/results/job_XXXX_nN_.../primary_results/global/U_global.csv
    pattern = str(FEM_RESULTS_DIR / "job_*" / "primary_results" / "global" / "U_global.csv")
    csv_files = sorted(glob.glob(pattern))
    by_base: dict[int, list[tuple[int, Path, str]]] = {}
    for csv_path in csv_files:
        csv_file = Path(csv_path)
        job_dir = csv_file.parent.parent.parent
        name = job_dir.name
        m = re.match(r"job_(?P<id>\d+)_n(?P<n>\d+)_[\d\-_]+_pid\d+_[a-f0-9]+", name)
        if not m:
            continue
        base_id = int(m.group("id"))
        n = int(m.group("n"))
        if base_id not in (0, 1, 2, 5, 6, 7):
            continue
        job_input_name = f"job_{base_id:04d}_n{n}"
        by_base.setdefault(base_id, []).append((n, csv_file, job_input_name))

    n_map_per_base: dict[int, dict[int, tuple[Path, str]]] = {}
    for base_id in by_base:
        n_map_per_base[base_id] = {}
        for n, cf, jn in sorted(by_base[base_id], key=lambda t: t[0]):
            if n not in n_map_per_base[base_id]:
                n_map_per_base[base_id][n] = (cf, jn)

    rows: list[dict] = []
    all_jobs = [(j, lbl, _) for j, lbl, _ in POINT_JOBS] + [(j, lbl, _) for j, lbl, _ in DIST_JOBS]

    for job_id, load_label, _ in all_jobs:
        if job_id not in n_map_per_base:
            continue
        n_map = n_map_per_base[job_id]
        if N_FINE not in n_map or N_MED not in n_map or N_COARSE not in n_map:
            for n_need in (N_FINE, N_MED, N_COARSE):
                if n_need not in n_map:
                    print(f"WARNING: job_{job_id:04d} missing FEM n={n_need}, skipping.")
                    break
            continue

        # Abaqus reference: fine mesh (n=100) only
        abaqus_job_name = f"job_{job_id:04d}_n{N_FINE}"
        abaqus_csv = ABAQUS_RESULTS_DIR / abaqus_job_name / "U_global.csv"
        if not abaqus_csv.is_file():
            print(f"WARNING: Abaqus result missing for {abaqus_job_name}, skipping job_{job_id:04d}.")
            continue

        abaqus_tip = _get_tip_values(abaqus_csv, JOBS_DIR, abaqus_job_name)
        if abaqus_tip is None:
            print(f"WARNING: failed to get Abaqus tip for {abaqus_job_name}, skipping.")
            continue
        _L_abaqus, abaqus_uy_mm, abaqus_theta_deg = abaqus_tip

        U_abaqus = _read_U_global(abaqus_csv)

        # Treat all-zero Abaqus deflection as missing reference; warn when rotation missing (UR not in ODB)
        abaqus_uy_mm_eff = abaqus_uy_mm
        abaqus_theta_deg_eff = abaqus_theta_deg
        if U_abaqus is not None:
            if np.max(np.abs(U_abaqus[:, 1])) < 1e-12:
                abaqus_uy_mm_eff = float("nan")
                print(f"WARNING: {abaqus_job_name} Abaqus U_global has all-zero UY (no valid deflection reference); error vs Abaqus set to N/A.")
            if np.max(np.abs(U_abaqus[:, 5])) < 1e-12 and job_id == 0:
                print("WARNING: Abaqus U_global has all-zero rotation (UR). Request U and UR in the Abaqus step and re-extract ODB for rotation comparison.")

        phi_uy: dict[int, float] = {}
        phi_theta: dict[int, float] = {}
        for n in (N_FINE, N_MED, N_COARSE):
            cf, jn = n_map[n]
            out = _get_tip_values(cf, JOBS_DIR, jn)
            if out is None:
                print(f"WARNING: job_{job_id:04d}_n{n} failed to get FEM tip values, skipping.")
                break
            _L, uy_mm, theta_deg = out
            phi_uy[n] = uy_mm
            phi_theta[n] = theta_deg
        else:
            for qoi_name, phi_dict, abaqus_val in [
                ("tip_deflection_mm", phi_uy, abaqus_uy_mm_eff),
                ("tip_rotation_deg", phi_theta, abaqus_theta_deg_eff),
            ]:
                p_obs, phi_ext, gci_12_pct, gci_23_pct = _richardson_gci(
                    phi_dict[N_FINE], phi_dict[N_MED], phi_dict[N_COARSE]
                )
                phi_fine = phi_dict[N_FINE]
                if math.isfinite(abaqus_val) and abs(abaqus_val) >= 1e-30:
                    err_ext_pct = 100.0 * (phi_ext - abaqus_val) / abaqus_val
                    error_phi1_pct = 100.0 * (phi_fine - abaqus_val) / abaqus_val
                else:
                    err_ext_pct = float("nan")
                    error_phi1_pct = float("nan")

                if abs(phi_ext) >= 1e-30 and math.isfinite(phi_ext):
                    fine_vs_rich_pct = 100.0 * abs(phi_fine - phi_ext) / abs(phi_ext)
                else:
                    fine_vs_rich_pct = float("nan")
                if (
                    math.isfinite(gci_12_pct)
                    and abs(gci_12_pct) >= 1e-30
                    and math.isfinite(gci_23_pct)
                ):
                    gci_23_12_ratio = gci_23_pct / gci_12_pct
                else:
                    gci_23_12_ratio = float("nan")
                r_p = (R ** p_obs) if math.isfinite(p_obs) else float("nan")

                rows.append({
                    "job_id": job_id,
                    "load_type": load_label,
                    "QoI": qoi_name,
                    "n_fine": N_FINE,
                    "n_med": N_MED,
                    "n_coarse": N_COARSE,
                    "phi_fine": phi_fine,
                    "phi_med": phi_dict[N_MED],
                    "phi_coarse": phi_dict[N_COARSE],
                    "p_obs": p_obs,
                    "phi_ext": phi_ext,
                    "Abaqus": abaqus_val,
                    "error_ext_Abaqus_pct": err_ext_pct,
                    "error_phi1_Abaqus_pct": error_phi1_pct,
                    "GCI_12_pct": gci_12_pct,
                    "GCI_23_pct": gci_23_pct,
                    "fine_vs_rich_pct": fine_vs_rich_pct,
                    "GCI_23_12_ratio": gci_23_12_ratio,
                    "r_p": r_p,
                })

    if not rows:
        print(
            "No GCI/Richardson vs Abaqus data. Need FEM results for job_0,1,2,5,6,7 at n=32,64,128 "
            "and Abaqus results at n=128 (run run_abaqus_cae.py for job_0000_n128, job_0001_n128, etc.)."
        )
        return

    csv_path = out_dir / "gci_richardson_abaqus_deflection_rotation.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(
            "job_id,load_type,QoI,n_fine,n_med,n_coarse,phi_fine,phi_med,phi_coarse,p_obs,phi_ext,Abaqus,"
            "error_ext_Abaqus_pct,error_phi1_Abaqus_pct,GCI_12_pct,GCI_23_pct,fine_vs_rich_pct,GCI_23_12_ratio,r_p\n"
        )
        for r in rows:
            f.write(
                f"{r['job_id']},{r['load_type']},{r['QoI']},{r['n_fine']},{r['n_med']},{r['n_coarse']},"
                f"{r['phi_fine']:.12e},{r['phi_med']:.12e},{r['phi_coarse']:.12e},"
                f"{r['p_obs']:.6f},{r['phi_ext']:.12e},{r['Abaqus']:.12e},"
                f"{r['error_ext_Abaqus_pct']:.6f},{r['error_phi1_Abaqus_pct']:.6f},"
                f"{r['GCI_12_pct']:.6f},{r['GCI_23_pct']:.6f},"
                f"{r['fine_vs_rich_pct']:.6f},{r['GCI_23_12_ratio']:.6f},{r['r_p']:.6f}\n"
            )
    print(f"Saved: {csv_path}")

    print("\n--- GCI and Richardson extrapolation (Abaqus vs FEM tip deflection/rotation) ---")
    for r in rows:
        p_str = f"{r['p_obs']:.3f}" if math.isfinite(r["p_obs"]) else "N/A"
        err_str = f"{r['error_ext_Abaqus_pct']:.4f}%" if math.isfinite(r["error_ext_Abaqus_pct"]) else str(r["error_ext_Abaqus_pct"])
        gci12_str = f"{r['GCI_12_pct']:.4f}%" if math.isfinite(r["GCI_12_pct"]) else "N/A"
        gci23_str = f"{r['GCI_23_pct']:.4f}%" if math.isfinite(r["GCI_23_pct"]) else "N/A"
        print(
            f"  job_{r['job_id']:04d} {r['load_type']} {r['QoI']}: p_obs={p_str} phi_ext={r['phi_ext']:.6e} "
            f"Abaqus={r['Abaqus']:.6e} error_ext={err_str} GCI_12={gci12_str} GCI_23={gci23_str}"
        )
    print()


if __name__ == "__main__":
    _run_report()
    sys.exit(0)
