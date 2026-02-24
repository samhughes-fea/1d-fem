# post_processing/verification_visualisers/shear_deformable_verification.py
"""
Timoshenko and Levinson verification: compare FEM tip u_y and θ_z to analytical.

Discovers job_0003_n* (Timoshenko) and job_0004_n* (Levinson) results,
reads tip displacement/rotation from U_global, and compares to analytical
tip deflection and rotation (same formulas as tests/analytical_*_benchmark.py).
Output: overlay plot (u_y, θ_z along span with analytical tip reference),
tip comparison table, and optional CSV.
"""

import glob
import re
import sys
from pathlib import Path
from typing import Final, Optional

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend so script saves and exits without blocking
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = next(
    (p for p in SCRIPT_DIR.parents if (p / "pre_processing").is_dir()),
    SCRIPT_DIR.parents[4],
)
sys.path.insert(0, str(PROJECT_ROOT))

from pre_processing.parsing.grid_parser import GridParser  # type: ignore

from analytical_timoshenko_benchmark import analytical_timoshenko_tip_deflection  # type: ignore
from analytical_levinson_benchmark import analytical_levinson_tip_deflection  # type: ignore

# Match job parameters (E, G, A, I_z, L, P, kappa) from jobs
P: Final[float] = -500.0  # N
E: Final[float] = 2.1e11  # Pa
G: Final[float] = 8.1e10  # Pa
A: Final[float] = 0.00131  # m^2
I_z: Final[float] = 2.08769e-06  # m^4
L: Final[float] = 2.0  # m
KAPPA: Final[float] = 5.0 / 6.0

# Analytical tip rotation (same for both: θ = PL^2/(2EI))
def _analytical_tip_rotation() -> float:
    return (P * L**2) / (2 * E * I_z)


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


def run_shear_deformable_verification() -> None:
    results_dir: Path = PROJECT_ROOT / "post_processing" / "results"
    jobs_dir: Path = PROJECT_ROOT / "jobs"
    out_dir: Path = SCRIPT_DIR / "shear_deformable_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
    csv_files = sorted(glob.glob(pattern))
    by_base: dict[int, list[tuple[int, Path, str]]] = {}
    for csv_path in csv_files:
        csv_file = Path(csv_path)
        job_dir = csv_file.parent.parent.parent
        name = job_dir.name
        m = re.match(r"job_(?P<id>\d+)_n(?P<n>\d+)_[\d\-_]+_pid\d+_[a-f0-9]+", name)
        if m:
            base_id = int(m.group("id"))
            n = int(m.group("n"))
            job_input_name = f"job_{base_id:04d}_n{n}"
        else:
            m2 = re.match(r"job_(?P<id>\d+)_[\d\-_]+_pid\d+_[a-f0-9]+", name)
            if m2:
                base_id = int(m2.group("id"))
                n = 0
                job_input_name = f"job_{base_id:04d}"
            else:
                continue
        if base_id not in (3, 4):
            continue
        by_base.setdefault(base_id, []).append((n, csv_file, job_input_name))
    # Keep best n per base
    job_to_result: dict[int, tuple[Path, str]] = {}
    for base_id in (3, 4):
        if base_id not in by_base:
            continue
        candidates = sorted(by_base[base_id], key=lambda t: -t[0])
        _, csv_file, job_input_name = candidates[0]
        job_to_result[base_id] = (csv_file, job_input_name)

    if not job_to_result:
        print("No Timoshenko/Levinson result dirs (job_0003_n*, job_0004_n*) found. Run jobs first.")
        return

    theta_analytical_rad = _analytical_tip_rotation()
    theta_analytical_deg = np.degrees(theta_analytical_rad)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle("Shear-deformable verification: tip u_y and θ_z vs analytical", fontsize=12, fontweight="bold")

    csv_rows: list[list[float]] = []
    labels = []
    fem_uy = []
    ana_uy = []
    fem_theta = []
    ana_theta = []

    cases = [
        (3, "Timoshenko", analytical_timoshenko_tip_deflection(P, E, I_z, L, G, A, KAPPA)[0]),
        (4, "Levinson", analytical_levinson_tip_deflection(P, E, I_z, L, G, A)[0]),
    ]
    for job_id, title, ana_deflection_m in cases:
        ana_uy_mm = ana_deflection_m * 1000
        ana_theta_deg = theta_analytical_deg
        if job_id not in job_to_result:
            print(f"WARNING: No result for job_{job_id:04d}, skipping '{title}'.")
            continue
        csv_file, job_input_name = job_to_result[job_id]
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
        tip_idx = int(np.argmax(x))
        uy_fem_mm = U[tip_idx, 1] * 1000
        theta_fem_deg = np.degrees(U[tip_idx, 5])
        err_uy = uy_fem_mm - ana_uy_mm
        err_theta = theta_fem_deg - ana_theta_deg
        labels.append(title)
        fem_uy.append(uy_fem_mm)
        ana_uy.append(ana_uy_mm)
        fem_theta.append(theta_fem_deg)
        ana_theta.append(ana_theta_deg)
        csv_rows.append([float(job_id), uy_fem_mm, ana_uy_mm, err_uy, theta_fem_deg, ana_theta_deg, err_theta])
        print(f"  job_{job_id:04d} ({title}): u_y FEM={uy_fem_mm:.6f} mm, analytical={ana_uy_mm:.6f} mm, error={err_uy:.6f} mm")
        print(f"    θ_z FEM={theta_fem_deg:.6f} deg, analytical={ana_theta_deg:.6f} deg, error={err_theta:.6f} deg")

    if not labels:
        plt.close(fig)
        return

    x_pos = np.arange(len(labels))
    width = 0.35
    axes[0].bar(x_pos - width / 2, fem_uy, width, label="FEM", color="#4F81BD")
    axes[0].bar(x_pos + width / 2, ana_uy, width, label="Analytical", color="#C0504D")
    axes[0].set_ylabel(r"Tip $u_y$ [mm]")
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(labels)
    axes[0].legend()
    axes[0].set_title("Tip deflection")
    axes[0].grid(axis="y", ls="--", alpha=0.6)

    axes[1].bar(x_pos - width / 2, fem_theta, width, label="FEM", color="#4F81BD")
    axes[1].bar(x_pos + width / 2, ana_theta, width, label="Analytical", color="#C0504D")
    axes[1].set_ylabel(r"Tip $\theta_z$ [deg]")
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(labels)
    axes[1].legend()
    axes[1].set_title("Tip rotation")
    axes[1].grid(axis="y", ls="--", alpha=0.6)

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)
    plot_path = out_dir / "shear_deformable_tip_verification.png"
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {plot_path}")

    if csv_rows:
        csv_path = out_dir / "shear_deformable_verification_data.csv"
        header = "job_id,uy_fem_mm,uy_analytical_mm,error_uy_mm,theta_z_fem_deg,theta_z_analytical_deg,error_theta_deg"
        np.savetxt(csv_path, csv_rows, delimiter=",", header=header, comments="")
        print(f"Saved: {csv_path}")
    print("Shear-deformable verification done.")


if __name__ == "__main__":
    run_shear_deformable_verification()
