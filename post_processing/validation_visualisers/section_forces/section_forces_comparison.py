# post_processing/validation_visualisers/section_forces/section_forces_comparison.py
"""
Compare section forces (SFD/BMD: Vy, Mz) between FEM and Abaqus.
Loads FEM from tertiary_results (nodal or gaussian); Abaqus from
abaqus_results/job_XXXX_nX/nodal_section_forces.csv if present (same format as FEM),
else section_forces.csv (element x, N, Vy, Vz, T, My, Mz) with interpolation to nodes.
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
from pre_processing.parsing.element_parser import ElementParser

FEM_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
JOBS_DIR: Final[Path] = PROJECT_ROOT / "jobs"
ABAQUS_RESULTS_DIR: Final[Path] = VALIDATION_DIR / "abaqus_results"
OUT_DIR: Final[Path] = VALIDATION_DIR / "output"

IDX_VY, IDX_MZ = 1, 5  # Column order: N, Vy, Vz, T, My, Mz


def _read_fem_nodal_section_forces(csv_path: Path) -> Optional[np.ndarray]:
    """(n_nodes, 6) or None."""
    try:
        with open(csv_path, encoding="utf-8") as f:
            first = f.readline()
        skip = 2 if "column_order" in first else 1
        data = np.genfromtxt(csv_path, delimiter=",", skip_header=skip)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] != 6:
            return None
        return data
    except Exception:
        return None


def _gather_fem_section_forces_from_gaussian(
    job_dir: Path,
    element_dictionary: dict,
    grid_dictionary: dict,
    n_nodes: int,
) -> Optional[tuple[np.ndarray, np.ndarray]]:
    """(x, forces) forces shape (n_nodes, 6)."""
    section_forces_dir = job_dir / "tertiary_results" / "gaussian" / "section_forces"
    if not section_forces_dir.is_dir():
        return None
    files = sorted(
        section_forces_dir.glob("section_forces_elem_*.csv"),
        key=lambda p: int(p.stem.split("_")[-1]),
    )
    if not files:
        return None
    coords = grid_dictionary.get("coordinates")
    if coords is None:
        return None
    x_all = coords[:, 0]
    connectivity = element_dictionary.get("connectivity")
    forces_nodal = np.zeros((n_nodes, 6))
    weight = np.zeros(n_nodes)
    for i, csv_path in enumerate(files):
        try:
            with open(csv_path, encoding="utf-8") as f:
                first = f.readline()
                skip = 2 if "column_order" in first else 1
            data = np.genfromtxt(csv_path, delimiter=",", skip_header=skip)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] != 6:
                continue
        except Exception:
            continue
        node_ids = connectivity[i]
        elem_mean = np.mean(data, axis=0)
        for nid in node_ids:
            if nid < n_nodes:
                forces_nodal[nid] += elem_mean
                weight[nid] += 1.0
    if not np.any(weight > 0):
        return None
    nonzero = weight > 0
    forces_nodal[nonzero] /= weight[nonzero, np.newaxis]
    return (x_all, forces_nodal)


def _read_abaqus_section_forces(csv_path: Path) -> Optional[tuple[np.ndarray, np.ndarray]]:
    """
    Read Abaqus section forces CSV. Expected format: x or node_id, N, Vy, Vz, T, My, Mz
    or one row per node. Returns (x_or_node, (n, 6)) or None.
    """
    try:
        data = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 7:
            return None
        x = data[:, 0]
        forces = data[:, 1:7]
        return (x, forces)
    except Exception:
        return None


def _discover_pairs() -> list[tuple[str, Path, Path, Optional[Path]]]:
    """(job_input_name, fem_result_dir, grid_file, abaqus_section_csv or None)."""
    by_key: dict[tuple[int, int], Path] = {}
    pattern = str(FEM_RESULTS_DIR / "job_*" / "primary_results" / "global" / "U_global.csv")
    for csv_path in glob.glob(pattern):
        fem_dir = Path(csv_path).parent.parent.parent
        name = fem_dir.name
        m = re.match(r"job_(\d+)_n(\d+)_[\d\-_]+_pid\d+_[a-f0-9]+", name)
        if not m:
            continue
        base_id, n = int(m.group(1)), int(m.group(2))
        key = (base_id, n)
        if key not in by_key or fem_dir.stat().st_mtime > by_key[key].stat().st_mtime:
            by_key[key] = fem_dir
    pairs = []
    for (base_id, n), fem_dir in sorted(by_key.items()):
        job_input_name = f"job_{base_id:04d}_n{n}"
        grid_file = JOBS_DIR / job_input_name / "grid.txt"
        if not grid_file.is_file():
            continue
        abaqus_csv = ABAQUS_RESULTS_DIR / job_input_name / "section_forces.csv"
        abaqus_csv = abaqus_csv if abaqus_csv.is_file() else None
        pairs.append((job_input_name, fem_dir, grid_file, abaqus_csv))
    return pairs


def run_section_forces_comparison() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pairs = _discover_pairs()
    if not pairs:
        print("No FEM result dirs found. Run jobs first.")
        return

    n_with_abaqus = sum(1 for (_, _, _, abaqus_csv) in pairs if abaqus_csv is not None)
    if n_with_abaqus == 0:
        print("No Abaqus section_forces.csv found in abaqus_results/. Plots show FEM only.")
        print("To get Abaqus curves, run: python post_processing/validation_visualisers/abaqus/run_abaqus_cae.py --job job_0000_n8 ...")

    for job_input_name, fem_dir, grid_file, abaqus_csv in pairs:
        grid = GridParser(str(grid_file), str(JOBS_DIR / job_input_name)).parse()
        gd = grid["grid_dictionary"]
        coords = gd["coordinates"]
        n_nodes = coords.shape[0]
        x = coords[:, 0]

        # FEM section forces: nodal first, else gaussian
        elem_file = JOBS_DIR / job_input_name / "element.txt"
        element_dictionary = ElementParser(str(elem_file), str(JOBS_DIR / job_input_name)).parse()["element_dictionary"]
        nodal_csv = fem_dir / "tertiary_results" / "nodal" / "nodal_section_forces.csv"
        fem_sf = _read_fem_nodal_section_forces(nodal_csv)
        if fem_sf is None:
            out = _gather_fem_section_forces_from_gaussian(fem_dir, element_dictionary, gd, n_nodes)
            if out is not None:
                x_fem, fem_sf = out
            else:
                print(f"Skip {job_input_name}: no FEM section forces")
                continue
        else:
            x_fem = x
        if fem_sf.shape[0] != n_nodes:
            print(f"Skip {job_input_name}: FEM section forces length mismatch")
            continue

        abaqus_nodal_csv = ABAQUS_RESULTS_DIR / job_input_name / "nodal_section_forces.csv"
        if abaqus_csv is None or not abaqus_csv.is_file():
            abaqus_nodal_csv = None  # no Abaqus SF at all
        if abaqus_nodal_csv is not None and not abaqus_nodal_csv.is_file():
            abaqus_nodal_csv = None

        if abaqus_csv is None or not abaqus_csv.is_file():
            # Plot FEM only
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            fig.suptitle(f"Section forces ({job_input_name}) — FEM only (no Abaqus SF)", fontsize=12)
            axes[0].plot(x_fem, fem_sf[:, IDX_VY], "b-o", label="FEM Vy", markersize=4)
            axes[0].set_xlabel("x (m)")
            axes[0].set_ylabel("Vy (N)")
            axes[0].legend()
            axes[0].grid(True, alpha=0.3)
            axes[1].plot(x_fem, fem_sf[:, IDX_MZ], "b-o", label="FEM Mz", markersize=4)
            axes[1].set_xlabel("x (m)")
            axes[1].set_ylabel("Mz (N·m)")
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)
        else:
            # Prefer Abaqus nodal_section_forces.csv (same format as FEM) for direct nodal comparison
            if abaqus_nodal_csv is not None:
                abaqus_sf = _read_fem_nodal_section_forces(abaqus_nodal_csv)
                if abaqus_sf is not None and abaqus_sf.shape[0] == n_nodes:
                    x_abaqus = x_fem
                else:
                    abaqus_sf = None
            else:
                abaqus_sf = None
            if abaqus_sf is None and abaqus_csv is not None:
                abaqus_out = _read_abaqus_section_forces(abaqus_csv)
                if abaqus_out is None:
                    print(f"Skip {job_input_name}: could not read Abaqus section forces")
                    continue
                x_abaqus, abaqus_sf = abaqus_out
                if abaqus_sf.shape[0] != n_nodes:
                    if x_abaqus.shape[0] < 2:
                        print(f"Skip {job_input_name}: Abaqus SF too few points")
                        continue
                    abaqus_sf = np.column_stack([
                        np.interp(x_fem, x_abaqus, abaqus_sf[:, j]) for j in range(6)
                    ])
            if abaqus_sf is None:
                print(f"Skip {job_input_name}: could not read Abaqus section forces")
                continue
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            fig.suptitle(f"Section forces: FEM vs Abaqus ({job_input_name})", fontsize=12, fontweight="bold")
            axes[0].plot(x_fem, fem_sf[:, IDX_VY], "b-o", label="FEM Vy", markersize=4)
            axes[0].plot(x_fem, abaqus_sf[:, IDX_VY], "r--s", label="Abaqus Vy", markersize=4)
            axes[0].set_xlabel("x (m)")
            axes[0].set_ylabel("Vy (N)")
            axes[0].legend()
            axes[0].grid(True, alpha=0.3)
            axes[1].plot(x_fem, fem_sf[:, IDX_MZ], "b-o", label="FEM Mz", markersize=4)
            axes[1].plot(x_fem, abaqus_sf[:, IDX_MZ], "r--s", label="Abaqus Mz", markersize=4)
            axes[1].set_xlabel("x (m)")
            axes[1].set_ylabel("Mz (N·m)")
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = OUT_DIR / f"section_forces_fem_vs_abaqus_{job_input_name}.png"
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Wrote {out_png}")

    print("Section forces comparison done.")


if __name__ == "__main__":
    run_section_forces_comparison()
