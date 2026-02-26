# post_processing/validation_visualisers/deflection_tables/u_global_largest_mesh_review.py
"""
U_global review: compare Abaqus to FEM at the largest mesh per base job.

Discovers all (job_XXXX_nN) pairs with both FEM and Abaqus U_global; for each base_id
keeps only the pair with maximum n. Computes full-field max/mean error per DOF and
tip u_y / theta_z comparison. Writes validation_visualisers/output/u_global_largest_mesh_review.csv.

Run from project root:
  python post_processing/validation_visualisers/deflection_tables/u_global_largest_mesh_review.py
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path
from typing import Final, Optional

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

DOF_NAMES = ["UX", "UY", "UZ", "RX", "RY", "RZ"]


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
        vals = data[:, 1].astype(float)
        if vals.size % 6 != 0:
            return None
        return vals.reshape(-1, 6)
    except Exception:
        return None


def _discover_pairs() -> list[tuple[str, Path, Path]]:
    """Return (job_input_name, fem_csv_path, abaqus_csv_path) for jobs that have both."""
    pairs = []
    if not ABAQUS_RESULTS_DIR.is_dir():
        return pairs
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
        pattern = str(
            FEM_RESULTS_DIR / f"job_{base_id:04d}_n{n}_*" / "primary_results" / "global" / "U_global.csv"
        )
        fem_files = glob.glob(pattern)
        if not fem_files:
            continue
        fem_files = [Path(p) for p in fem_files]
        fem_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        pairs.append((job_input_name, fem_files[0], abaqus_csv))
    return pairs


def _select_largest_mesh_per_base(
    pairs: list[tuple[str, Path, Path]]
) -> list[tuple[str, Path, Path]]:
    """Group by base_id, keep only the pair with maximum n per base."""
    by_base: dict[int, list[tuple[str, Path, Path]]] = {}
    for job_name, fem_csv, abaqus_csv in pairs:
        m = re.match(r"job_(\d+)_n(\d+)$", job_name)
        if not m:
            continue
        base_id, n = int(m.group(1)), int(m.group(2))
        if base_id not in by_base:
            by_base[base_id] = []
        by_base[base_id].append((job_name, fem_csv, abaqus_csv))
    result = []
    for base_id in sorted(by_base.keys()):
        group = by_base[base_id]
        # Sort by n descending, take first
        group_sorted = sorted(
            group,
            key=lambda t: int(re.match(r"job_\d+_n(\d+)$", t[0]).group(1)),
            reverse=True,
        )
        result.append(group_sorted[0])
    return result


def run_u_global_largest_mesh_review() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pairs = _discover_pairs()
    if not pairs:
        print("No FEM+Abaqus result pairs found. Run jobs and Abaqus to populate results.")
        return

    largest = _select_largest_mesh_per_base(pairs)
    rows = []

    for job_input_name, fem_csv, abaqus_csv in largest:
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
        n_nodes = coords.shape[0]
        if n_fem != n_nodes or n_abaqus != n_nodes:
            print(
                f"Skip {job_input_name}: node count mismatch (grid={n_nodes}, fem={n_fem}, abaqus={n_abaqus})"
            )
            continue

        m = re.match(r"job_(\d+)_n(\d+)$", job_input_name)
        base_id = int(m.group(1)) if m else -1
        n_elements = int(m.group(2)) if m else 0

        # Full-field errors per DOF (indices 0..5 = UX,UY,UZ,RX,RY,RZ)
        diff = np.abs(U_fem - U_abaqus)
        max_err_per_dof = [float(np.max(diff[:, j])) for j in range(6)]
        mean_err_per_dof = [float(np.mean(diff[:, j])) for j in range(6)]

        # Tip: node at max x
        x = coords[:, 0]
        tip_idx = int(np.argmax(x))
        tip_uy_fem = float(U_fem[tip_idx, 1])
        tip_uy_abaqus = float(U_abaqus[tip_idx, 1])
        tip_rz_fem = float(U_fem[tip_idx, 5])
        tip_rz_abaqus = float(U_abaqus[tip_idx, 5])
        tip_uy_abs = abs(tip_uy_fem - tip_uy_abaqus)
        tip_rz_abs = abs(tip_rz_fem - tip_rz_abaqus)
        tip_uy_rel_pct = (tip_uy_abs / abs(tip_uy_fem) * 100.0) if abs(tip_uy_fem) > 1e-30 else 0.0
        tip_rz_rel_pct = (tip_rz_abs / abs(tip_rz_fem) * 100.0) if abs(tip_rz_fem) > 1e-30 else 0.0

        row = [
            base_id,
            job_input_name,
            n_elements,
            n_nodes,
            *max_err_per_dof,
            *mean_err_per_dof,
            tip_uy_fem,
            tip_uy_abaqus,
            tip_uy_abs,
            tip_uy_rel_pct,
            tip_rz_fem,
            tip_rz_abaqus,
            tip_rz_abs,
            tip_rz_rel_pct,
        ]
        rows.append(row)

    if not rows:
        print("No pairs passed validation (grid/node mismatch).")
        return

    # CSV header
    max_cols = [f"max_err_{d}" for d in DOF_NAMES]
    mean_cols = [f"mean_err_{d}" for d in DOF_NAMES]
    header = (
        ["base_id", "job_name", "n_elements", "n_nodes"]
        + max_cols
        + mean_cols
        + [
            "tip_u_y_fem_m",
            "tip_u_y_abaqus_m",
            "tip_u_y_abs_err",
            "tip_u_y_rel_pct",
            "tip_theta_z_fem_rad",
            "tip_theta_z_abaqus_rad",
            "tip_theta_z_abs_err",
            "tip_theta_z_rel_pct",
        ]
    )
    out_csv = OUT_DIR / "u_global_largest_mesh_review.csv"
    with open(out_csv, "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(str(v) for v in row) + "\n")
    print(f"Wrote {out_csv} ({len(rows)} base jobs at largest mesh).")
    print("U_global largest-mesh review done.")


if __name__ == "__main__":
    run_u_global_largest_mesh_review()
