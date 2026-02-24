# post_processing/verification_visualisers/deflection_tables/distributed_load_convergence.py
"""
Distributed-load deformation convergence: compute error vs mesh density (n_elements).

Discovers result dirs matching job_XXXX_nX (base 5, 6, 7 = UDL, triangular, parabolic),
computes max |u_y − u_roark| and RMS at node positions, and writes CSV.
Requires job input dirs job_0005_n4, job_0005_n8, job_0005_n16, job_0005_n100 (and
similarly for 0006, 0007). Run jobs first, then this script from project root.
"""

import glob
import re
import sys
from pathlib import Path
from typing import Final, Optional

import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = next(
    (p for p in SCRIPT_DIR.parents if (p / "pre_processing").is_dir()),
    SCRIPT_DIR.parents[4],
)
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR.parent / "roarks_formulas"))

from pre_processing.parsing.grid_parser import GridParser  # type: ignore

# Import Roark distributed (script lives in roarks_formulas one level up in verification_visualisers)
_roarks = SCRIPT_DIR.parent / "roarks_formulas"
sys.path.insert(0, str(_roarks))
from roarks_formulas_distributed import roark_distributed_load_response  # type: ignore

E: Final[float] = 2.0e11   # Pa
I_z: Final[float] = 2.08769e-06  # m^4
w_dist: Final[float] = 500.0  # N/m

DIST_LOAD_CASES: Final[list[tuple[int, str, str]]] = [
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
    raise KeyError("grid data does not contain 'grid_dictionary' → 'coordinates'")


def _read_U_global(file: Path) -> Optional[np.ndarray]:
    try:
        vals = np.genfromtxt(file, delimiter=",", skip_header=1, usecols=1)
        if vals.size % 6:
            return None
        return vals.reshape(-1, 6)
    except Exception:
        return None


def run_distributed_convergence() -> None:
    results_dir: Path = PROJECT_ROOT / "post_processing" / "results"
    jobs_dir: Path = PROJECT_ROOT / "jobs"
    out_dir: Path = SCRIPT_DIR / "deformation_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
    csv_files = sorted(glob.glob(pattern))
    # Collect (base_id, n, csv_file, result_dir_name) for base_id in (5,6,7)
    runs: list[tuple[int, int, Path, str]] = []
    for csv_path in csv_files:
        csv_file = Path(csv_path)
        job_dir = csv_file.parent.parent.parent
        name = job_dir.name
        m = re.match(r"job_(?P<id>\d+)_n(?P<n>\d+)_[\d\-_]+_pid\d+_[a-f0-9]+", name)
        if not m:
            continue
        base_id = int(m.group("id"))
        n = int(m.group("n"))
        if base_id not in (5, 6, 7):
            continue
        job_input_name = f"job_{base_id:04d}_n{n}"
        runs.append((base_id, n, csv_file, job_input_name))

    if not runs:
        print("No distributed-load result dirs (job_0005_n*, job_0006_n*, job_0007_n*) found. Run jobs first.")
        return

    # For each (base_id, n) compute max |u_y - u_roark| and RMS
    data: list[tuple[int, str, int, float, float]] = []  # base_id, title, n, max_err_mm, rms_err_mm
    for base_id, title, roark_type in DIST_LOAD_CASES:
        case_runs = [(n, csv_file, job_input_name) for (bid, n, csv_file, job_input_name) in runs if bid == base_id]
        case_runs.sort(key=lambda t: t[0])
        for n, csv_file, job_input_name in case_runs:
            grid_file = jobs_dir / job_input_name / "grid.txt"
            if not grid_file.is_file():
                continue
            U = _read_U_global(csv_file)
            if U is None:
                continue
            grid = GridParser(str(grid_file), str(jobs_dir / job_input_name)).parse()
            try:
                node_coords = _get_node_coordinates(grid)
            except Exception:
                continue
            x = node_coords[:, 0]
            if x.shape[0] != U.shape[0]:
                continue
            L = float(np.max(x))
            order = np.argsort(x)
            x_sorted = x[order]
            roark_at_nodes = roark_distributed_load_response(x_sorted, L, E, I_z, w_dist, roark_type)
            uy_roark = np.interp(x, x_sorted, roark_at_nodes["deflection"] * 1000)  # mm
            uy_fem = U[:, 1] * 1000
            err = uy_fem - uy_roark
            max_err = float(np.max(np.abs(err)))
            rms_err = float(np.sqrt(np.mean(err ** 2)))
            data.append((base_id, title, n, max_err, rms_err))
            print(f"  job_{base_id:04d}_n{n} ({title}): max|u_y err|={max_err:.6f} mm, RMS={rms_err:.6f} mm")

    if not data:
        print("No valid runs to plot.")
        return

    csv_path = out_dir / "distributed_load_convergence_data.csv"
    header = "job_base_id,load_type,n_elements,max_error_uy_mm,rms_error_uy_mm"
    rows = [f"{bid},{title},{n},{max_err:.12e},{rms_err:.12e}" for (bid, title, n, max_err, rms_err) in data]
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        f.write("\n".join(rows))
    print(f"Saved: {csv_path}")
    print("Distributed-load convergence done.")


if __name__ == "__main__":
    run_distributed_convergence()
