# post_processing/validation_visualisers/deflection_tables/deformation_comparison.py
"""
Compare deformation (U) between FEM and Abaqus for matching jobs.
Discovers FEM results by job_XXXX_nN (full name; latest result dir) and Abaqus from
validation_visualisers/abaqus_results/job_XXXX_nN/. Overlays u_y and θ_z (and optionally all 6 DOFs).
Output: validation_visualisers/output/
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path
from typing import Final, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
VALIDATION_DIR: Final[Path] = SCRIPT_DIR.parent
PROJECT_ROOT: Final[Path] = VALIDATION_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.parsing.grid_parser import GridParser

FEM_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
JOBS_DIR: Final[Path] = PROJECT_ROOT / "jobs"
ABAQUS_RESULTS_DIR: Final[Path] = VALIDATION_DIR / "abaqus_results"
OUT_DIR: Final[Path] = VALIDATION_DIR / "output"

# Match timestamped FEM result dir: job_0000_n8_2026-02-22_..._pid123_abc
FEM_DIR_PATTERN = re.compile(r"job_(?P<base_id>\d+)_n(?P<n>\d+)_[\d\-_]+_pid\d+_[a-f0-9]+")
DOF_LABELS = ["UX", "UY", "UZ", "RX", "RY", "RZ"]


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


def _read_U_global_csv(file: Path) -> Optional[np.ndarray]:
    """Read U_global.csv; return (n_nodes, 6) array (UX,UY,UZ,RX,RY,RZ)."""
    try:
        data = np.genfromtxt(file, delimiter=",", skip_header=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 2:
            return None
        vals = data[:, 1].astype(float)  # Value column
        if vals.size % 6 != 0:
            return None
        return vals.reshape(-1, 6)
    except Exception:
        return None


def _discover_pairs() -> list[tuple[str, Path, Path]]:
    """
    Return list of (job_input_name, fem_csv_path, abaqus_csv_path) for jobs that have both.
    job_input_name = e.g. job_0000_n8.
    """
    pairs = []
    if not ABAQUS_RESULTS_DIR.is_dir():
        return pairs
    # Scan Abaqus results (one dir per job)
    for abaqus_job_dir in sorted(ABAQUS_RESULTS_DIR.iterdir()):
        if not abaqus_job_dir.is_dir():
            continue
        name = abaqus_job_dir.name
        m = re.match(r"job_(\d+)_n(\d+)$", name)
        if not m:
            continue
        base_id, n = int(m.group(1)), int(m.group(2))
        job_input_name = f"job_{base_id:04d}_n{n}"
        abaqus_csv = abaqus_job_dir / "U_global.csv"
        if not abaqus_csv.is_file():
            continue
        # Find one FEM result dir matching this job name (latest by mtime)
        pattern = str(FEM_RESULTS_DIR / f"job_{base_id:04d}_n{n}_*" / "primary_results" / "global" / "U_global.csv")
        fem_files = glob.glob(pattern)
        if not fem_files:
            continue
        fem_files = [Path(p) for p in fem_files]
        fem_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        pairs.append((job_input_name, fem_files[0], abaqus_csv))
    return pairs


def run_deformation_comparison(scale: float = 1.0) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pairs = _discover_pairs()
    if not pairs:
        print("No FEM+Abaqus result pairs found.")
        if not ABAQUS_RESULTS_DIR.is_dir() or not any(ABAQUS_RESULTS_DIR.iterdir()):
            print("  No Abaqus results in abaqus_results/. To get Abaqus on the plots, run:")
            print("  python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n8")
            print("  (and repeat --job for other jobs). Results go to abaqus_results/job_XXXX_nN/U_global.csv")
        else:
            print("  Ensure FEM results exist for the same job names as in abaqus_results/ (e.g. job_0000_n8).")
        return

    error_rows = []
    for job_input_name, fem_csv, abaqus_csv in pairs:
        grid_file = JOBS_DIR / job_input_name / "grid.txt"
        if not grid_file.is_file():
            print(f"Skip {job_input_name}: no grid.txt")
            continue
        grid = GridParser(str(grid_file), str(JOBS_DIR / job_input_name)).parse()
        try:
            coords = _get_node_coordinates(grid)
        except KeyError:
            print(f"Skip {job_input_name}: no coordinates in grid")
            continue
        U_fem = _read_U_global_csv(fem_csv)
        U_abaqus = _read_U_global_csv(abaqus_csv)
        if U_fem is None or U_abaqus is None:
            print(f"Skip {job_input_name}: failed to read U")
            continue
        n_fem, n_abaqus = U_fem.shape[0], U_abaqus.shape[0]
        if n_fem != coords.shape[0] or n_abaqus != coords.shape[0]:
            print(f"Skip {job_input_name}: node count mismatch (grid={coords.shape[0]}, fem={n_fem}, abaqus={n_abaqus})")
            continue

        x = coords[:, 0]
        L = float(np.max(x))

        # Plot u_y (index 1) and θ_z (index 5)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        fig.suptitle(f"Deformation: FEM vs Abaqus ({job_input_name})", fontsize=12, fontweight="bold")

        axes[0].plot(x, U_fem[:, 1] * scale, "b-o", label="FEM", markersize=4)
        axes[0].plot(x, U_abaqus[:, 1] * scale, "r--s", label="Abaqus", markersize=4)
        axes[0].set_xlabel("x (m)")
        axes[0].set_ylabel("u_y (m)")
        axes[0].set_title("Transverse displacement")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(x, U_fem[:, 5] * scale, "b-o", label="FEM", markersize=4)
        axes[1].plot(x, U_abaqus[:, 5] * scale, "r--s", label="Abaqus", markersize=4)
        axes[1].set_xlabel("x (m)")
        axes[1].set_ylabel("θ_z (rad)")
        axes[1].set_title("Rotation")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = OUT_DIR / f"deformation_fem_vs_abaqus_{job_input_name}.png"
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Wrote {out_png}")

        # Collect error row for CSV
        err_uy = np.abs(U_fem[:, 1] - U_abaqus[:, 1])
        err_rz = np.abs(U_fem[:, 5] - U_abaqus[:, 5])
        error_rows.append([
            job_input_name,
            float(np.max(err_uy)), float(np.mean(err_uy)),
            float(np.max(err_rz)), float(np.mean(err_rz)),
        ])

    if error_rows:
        err_csv = OUT_DIR / "deformation_comparison_errors.csv"
        with open(err_csv, "w") as f:
            f.write("job,max_err_uy,mean_err_uy,max_err_rz,mean_err_rz\n")
            for row in error_rows:
                f.write(",".join(str(v) for v in row) + "\n")
        print(f"Wrote {err_csv}")

    print("Deformation comparison done.")


if __name__ == "__main__":
    run_deformation_comparison()
