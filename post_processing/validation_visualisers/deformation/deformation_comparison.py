# post_processing/validation_visualisers/deformation/deformation_comparison.py
"""
Compare deformation (U) between FEM (n128) and Abaqus reference (n500).
FEM results at job_XXXX_n128; Abaqus converged reference at job_XXXX_n500.
Abaqus U is interpolated to FEM node positions for overlay. Overlays u_y and θ_z.
Output: validation_visualisers/deformation/deformation_plots/
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
from pre_processing.parsing.element_parser import ElementParser

from post_processing.validation_visualisers.abaqus.config import ELEMENT_TYPE_MAP, ABAQUS_REFERENCE_N

FEM_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
JOBS_DIR: Final[Path] = PROJECT_ROOT / "jobs"
ABAQUS_RESULTS_DIR: Final[Path] = VALIDATION_DIR / "abaqus_results"
OUT_DIR: Final[Path] = VALIDATION_DIR / "deformation" / "deformation_plots"

# FEM mesh for key comparison; Abaqus reference at ABAQUS_REFERENCE_N (500)
N_ELEMENTS_FILTER: Final[int] = 128

# Job base_id -> load case description for plot titles (aligned with verification/validation job mapping)
LOAD_CASE_BY_JOB_ID: Final[dict[int, str]] = {
    0: "Point (end)",
    1: "Point (mid-span)",
    2: "Point (quarter-span)",
    3: "Point (end)",
    4: "Point (mid-span)",
    5: "UDL",
    6: "Triangular",
    7: "Parabolic",
    8: "Point (end)",
    9: "UDL",
    10: "Triangular",
    11: "Parabolic",
}

# Match timestamped FEM result dir: job_0000_n8_2026-02-22_..._pid123_abc
FEM_DIR_PATTERN = re.compile(r"job_(?P<base_id>\d+)_n(?P<n>\d+)_[\d\-_]+_pid\d+_[a-f0-9]+")
DOF_LABELS = ["UX", "UY", "UZ", "RX", "RY", "RZ"]


def _get_load_case_label(job_input_name: str) -> str:
    """Return load case description for job (e.g. 'Point (end)', 'UDL')."""
    m = re.match(r"job_(\d+)_n\d+$", job_input_name)
    if not m:
        return "—"
    base_id = int(m.group(1))
    return LOAD_CASE_BY_JOB_ID.get(base_id, f"Load case {base_id}")


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


def _get_element_labels(job_input_name: str) -> tuple[str, str]:
    """Return (FEM element name, Abaqus element type) for plot labels, e.g. (EulerBernoulliBeamElement3D, B33)."""
    elem_file = JOBS_DIR / job_input_name / "element.txt"
    if not elem_file.is_file():
        return ("FEM", "Abaqus")
    try:
        ed = ElementParser(str(elem_file), str(JOBS_DIR / job_input_name)).parse()["element_dictionary"]
        types = ed.get("types")
        if types is None or len(types) == 0:
            return ("FEM", "Abaqus")
        elem_type = str(types[0])
        abaqus_elem = ELEMENT_TYPE_MAP.get(elem_type)
        if abaqus_elem is None:
            return ("FEM", "Abaqus")
        return (elem_type, abaqus_elem)
    except Exception:
        return ("FEM", "Abaqus")


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


def _discover_pairs() -> list[tuple[str, Path, Path, str]]:
    """
    Return list of (job_input_name, fem_csv_path, abaqus_csv_path, abaqus_job_name) for jobs that have both.
    FEM at n128; Abaqus reference at n500. job_input_name = e.g. job_0000_n128; abaqus_job_name = job_0000_n500.
    """
    pairs = []
    if not ABAQUS_RESULTS_DIR.is_dir():
        return pairs
    # For each base_id: require FEM job_XXXX_n128 and Abaqus job_XXXX_n500 (reference)
    for base_id in range(12):
        job_input_name = f"job_{base_id:04d}_n{N_ELEMENTS_FILTER}"
        abaqus_job_name = f"job_{base_id:04d}_n{ABAQUS_REFERENCE_N}"
        abaqus_csv = ABAQUS_RESULTS_DIR / abaqus_job_name / "U_global.csv"
        if not abaqus_csv.is_file():
            continue
        pattern = str(FEM_RESULTS_DIR / f"job_{base_id:04d}_n{N_ELEMENTS_FILTER}_*" / "primary_results" / "global" / "U_global.csv")
        fem_files = glob.glob(pattern)
        if not fem_files:
            continue
        fem_files = [Path(p) for p in fem_files]
        fem_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        pairs.append((job_input_name, fem_files[0], abaqus_csv, abaqus_job_name))
    return pairs


def run_deformation_comparison(scale: float = 1.0) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pairs = _discover_pairs()
    if not pairs:
        print(f"No FEM (n128) + Abaqus (n{ABAQUS_REFERENCE_N}) result pairs found.")
        if not ABAQUS_RESULTS_DIR.is_dir() or not any(ABAQUS_RESULTS_DIR.iterdir()):
            print("  No Abaqus results. Run Abaqus for reference jobs, e.g.:")
            print(f"  python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n{ABAQUS_REFERENCE_N}")
        else:
            print(f"  Ensure FEM job_XXXX_n{N_ELEMENTS_FILTER} and Abaqus job_XXXX_n{ABAQUS_REFERENCE_N} exist.")
        return

    error_rows = []
    for job_input_name, fem_csv, abaqus_csv, abaqus_job_name in pairs:
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
        U_abaqus_raw = _read_U_global_csv(abaqus_csv)
        if U_fem is None or U_abaqus_raw is None:
            print(f"Skip {job_input_name}: failed to read U")
            continue
        n_fem = U_fem.shape[0]
        if n_fem != coords.shape[0]:
            print(f"Skip {job_input_name}: FEM node count mismatch")
            continue
        # Abaqus reference (n500) has different mesh; interpolate to FEM x
        grid_abaqus_file = JOBS_DIR / abaqus_job_name / "grid.txt"
        if not grid_abaqus_file.is_file():
            print(f"Skip {job_input_name}: no grid for {abaqus_job_name}")
            continue
        grid_abaqus = GridParser(str(grid_abaqus_file), str(JOBS_DIR / abaqus_job_name)).parse()
        try:
            coords_abaqus = _get_node_coordinates(grid_abaqus)
        except KeyError:
            print(f"Skip {job_input_name}: no coordinates in Abaqus grid")
            continue
        x_fem = coords[:, 0]
        x_abaqus = coords_abaqus[:, 0]
        U_abaqus = np.column_stack([
            np.interp(x_fem, x_abaqus, U_abaqus_raw[:, j]) for j in range(6)
        ])

        x = x_fem
        L = float(np.max(x))

        # Plot: FEM n128 vs Abaqus n500 (reference)
        fem_label, abaqus_label = _get_element_labels(job_input_name)
        load_case = _get_load_case_label(job_input_name)
        abaqus_legend = f"{abaqus_label} (n={ABAQUS_REFERENCE_N} ref)"
        u_y_mm = 1000.0  # m -> mm
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        fig.suptitle(f"Deformation: {fem_label} (n128) vs {abaqus_label} (n={ABAQUS_REFERENCE_N} ref) — {load_case} ({job_input_name})", fontsize=12, fontweight="bold")

        axes[0].plot(x, U_fem[:, 1] * u_y_mm * scale, "b-o", label=fem_label, markersize=4)
        axes[0].plot(x, U_abaqus[:, 1] * u_y_mm * scale, "r--s", label=abaqus_legend, markersize=4)
        axes[0].set_xlabel("x (m)")
        axes[0].set_ylabel("u_y (mm)")
        axes[0].set_title("Transverse displacement")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(x, np.degrees(U_fem[:, 5]) * scale, "b-o", label=fem_label, markersize=4)
        axes[1].plot(x, np.degrees(U_abaqus[:, 5]) * scale, "r--s", label=abaqus_legend, markersize=4)
        axes[1].set_xlabel("x (m)")
        axes[1].set_ylabel("θ_z (deg)")
        axes[1].set_title("Rotation")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = OUT_DIR / f"deformation_fem_vs_abaqus_{job_input_name}.png"
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Wrote {out_png}")

        # Collect error row for CSV (FEM vs interpolated Abaqus reference)
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
