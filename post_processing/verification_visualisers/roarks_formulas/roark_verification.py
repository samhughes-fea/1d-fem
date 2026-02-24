# post_processing/verification_visualisers/roarks_formulas/roark_verification.py
"""
FEM vs Roark verification: compare displacement/rotation (and optionally V, M) from
FEM results with Roark's formulas. Uses same job discovery and grid as
deflection_tables/deformation_convergence.py. Aligned with jobs/README_JOBS.md:
  - Point loads (Euler–Bernoulli): job_0000 end, job_0001 midspan, job_0002 quarter.
  - Distributed loads: job_0005 UDL, job_0006 triangular, job_0007 parabolic.
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
sys.path.insert(0, str(SCRIPT_DIR))

from pre_processing.parsing.grid_parser import GridParser  # type: ignore

from roarks_formulas_point import roark_point_load_response  # type: ignore
from roarks_formulas_distributed import roark_distributed_load_response  # type: ignore


# Beam and load parameters (match deflection_tables/deformation_convergence for point;
# distributed magnitude as in visualiser)
E: Final[float] = 2.0e11   # Pa
I_z: Final[float] = 2.08769e-06  # m^4
P_point: Final[float] = -500.0  # N (match verification convention: same sign as deformation_convergence F)
w_dist: Final[float] = 500.0  # N/m; must match jobs/job_0005, 0006, 0007 distributed_load.txt (F_y magnitude)

# Analytical curve: much finer grid than FEM for plotting
_ANALYTICAL_GRID_FACTOR: Final[int] = 25
_ANALYTICAL_GRID_MIN: Final[int] = 500

# Job layout per README_JOBS.md
POINT_LOAD_CASES: Final[list[tuple[int, str, float]]] = [
    (0, "End", 1.0),
    (1, "Mid-span", 0.5),
    (2, "Quarter-span", 0.25),
]
DIST_LOAD_CASES: Final[list[tuple[int, str, str]]] = [
    (5, "UDL", "udl"),
    (6, "Triangular", "triangular"),
    (7, "Parabolic", "parabolic"),
]

COLORS: Final[list[str]] = ["#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#4BACC6"]


def _analytical_grid(n_fem: int, L: float) -> np.ndarray:
    n = max(n_fem * _ANALYTICAL_GRID_FACTOR, _ANALYTICAL_GRID_MIN)
    return np.linspace(0.0, L, int(n))


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
            raise ValueError("DOF count not divisible by 6")
        return vals.reshape(-1, 6)
    except Exception as exc:
        print(f"Error reading {file}: {exc}")
        return None


def run_roark_verification(scale: float = 1.0) -> None:
    results_dir: Path = PROJECT_ROOT / "post_processing" / "results"
    jobs_dir: Path = PROJECT_ROOT / "jobs"
    out_dir: Path = SCRIPT_DIR / "verification"
    out_dir.mkdir(exist_ok=True)

    pattern = str(results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
    csv_files = sorted(glob.glob(pattern))
    # (base_id, n) -> (csv_file, job_dir, job_input_name); keep largest n per base_id
    by_base: dict[int, list[tuple[int, Path, Path, str]]] = {}
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
        by_base.setdefault(base_id, []).append((n, csv_file, job_dir, job_input_name))
    job_to_result: dict[int, tuple[Path, Path, str]] = {}
    for base_id, candidates in by_base.items():
        candidates.sort(key=lambda t: -t[0])
        n, csv_file, job_dir, job_input_name = candidates[0]
        job_to_result[base_id] = (csv_file, job_dir, job_input_name)

    if not job_to_result:
        print("No FEM result directories found. Run jobs first.")
        return

    # ----- Point loads: overlay u_y, θ_z (FEM vs Roark) -----
    fig_point, axes_point = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
    fig_point.suptitle("Roark verification: point loads (Euler–Bernoulli)", fontsize=14, fontweight="bold")

    csv_rows: list[list[float]] = []

    for row, (job_id, title, a_frac) in enumerate(POINT_LOAD_CASES):
        if job_id not in job_to_result:
            print(f"WARNING: No result for job_{job_id:04d}, skipping '{title}'.")
            continue
        csv_file, job_dir, job_input_name = job_to_result[job_id]
        grid_file = jobs_dir / job_input_name / "grid.txt"
        if not grid_file.is_file():
            print(f"WARNING: Missing {grid_file}, skipping job_{job_id:04d}.")
            continue
        U = _read_U_global(csv_file)
        if U is None:
            continue
        grid = GridParser(str(grid_file), str(jobs_dir / job_input_name)).parse()
        try:
            node_coords = _get_node_coordinates(grid)
        except Exception as exc:
            print(f"WARNING: {exc} job_{job_id:04d}")
            continue
        x = node_coords[:, 0]
        if x.shape[0] != U.shape[0]:
            print(f"ERROR: x/U length mismatch job_{job_id:04d}")
            continue
        L = float(np.max(x))
        a = a_frac * L
        load_type = "end" if a_frac >= 1.0 else ("mid" if a_frac >= 0.5 else "quarter")

        uy_fem = U[:, 1] * 1000 * scale  # mm
        th_fem = np.degrees(U[:, 5]) * scale  # deg

        # Roark at FEM nodes (for errors and CSV)
        roark_at_nodes = roark_point_load_response(x, L, E, I_z, P_point, load_type)
        uy_roark = roark_at_nodes["deflection"] * 1000   # m -> mm
        th_roark = np.degrees(roark_at_nodes["rotation"])

        # Fine grid for smooth analytical curve
        x_ana = _analytical_grid(x.shape[0], L)
        roark_fine = roark_point_load_response(x_ana, L, E, I_z, P_point, load_type)
        uy_ana = roark_fine["deflection"] * 1000
        th_ana = np.degrees(roark_fine["rotation"])

        axes_point[row, 0].plot(x, uy_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_point[row, 0].plot(x_ana, uy_ana, "k--", label="Roark")
        axes_point[row, 0].set_ylabel(r"$u_y$ [mm]")
        axes_point[row, 0].set_title(title, fontweight="bold")
        axes_point[row, 0].grid(ls="--", alpha=0.6)
        axes_point[row, 0].legend(loc="upper right", fontsize="small")

        axes_point[row, 1].plot(x, th_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_point[row, 1].plot(x_ana, th_ana, "k--", label="Roark")
        axes_point[row, 1].set_ylabel(r"$\theta_z$ [deg]")
        axes_point[row, 1].set_title(title, fontweight="bold")
        axes_point[row, 1].grid(ls="--", alpha=0.6)
        axes_point[row, 1].legend(loc="upper right", fontsize="small")

        err_uy = uy_fem - uy_roark
        err_th = th_fem - th_roark
        for i in range(len(x)):
            csv_rows.append([
                float(job_id), x[i], uy_fem[i], uy_roark[i], err_uy[i],
                th_fem[i], th_roark[i], err_th[i],
            ])
        max_uy = np.max(np.abs(err_uy))
        rms_uy = np.sqrt(np.mean(err_uy ** 2))
        max_th = np.max(np.abs(err_th))
        rms_th = np.sqrt(np.mean(err_th ** 2))
        print(f"  job_{job_id:04d} ({title}): u_y max|err|={max_uy:.6f} mm, RMS={rms_uy:.6f} mm; "
              f"theta_z max|err|={max_th:.6f} deg, RMS={rms_th:.6f} deg")

    axes_point[-1, 0].set_xlabel(r"$x$ [m]")
    axes_point[-1, 1].set_xlabel(r"$x$ [m]")
    fig_point.tight_layout()
    fig_point.subplots_adjust(top=0.92)
    fig_point.savefig(out_dir / "roark_verification_point_loads.png", dpi=300)
    plt.close(fig_point)
    print(f"Saved: {out_dir / 'roark_verification_point_loads.png'}")

    # ----- Distributed loads: overlay u_y, θ_z -----
    fig_dist, axes_dist = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
    fig_dist.suptitle("Roark verification: distributed loads (Euler–Bernoulli)", fontsize=14, fontweight="bold")

    for row, (job_id, title, roark_type) in enumerate(DIST_LOAD_CASES):
        if job_id not in job_to_result:
            print(f"WARNING: No result for job_{job_id:04d}, skipping '{title}'.")
            continue
        csv_file, job_dir, job_input_name = job_to_result[job_id]
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
        # Distributed x must be sorted for Roark's cumulative_trapezoid
        order = np.argsort(x)
        x_sorted = x[order]
        roark_at_nodes = roark_distributed_load_response(x_sorted, L, E, I_z, w_dist, roark_type)
        uy_roark = np.interp(x, x_sorted, roark_at_nodes["deflection"] * 1000)
        th_roark = np.interp(x, x_sorted, np.degrees(roark_at_nodes["rotation"]))

        uy_fem = U[:, 1] * 1000 * scale
        th_fem = np.degrees(U[:, 5]) * scale

        x_ana = _analytical_grid(x.shape[0], L)
        x_ana_sorted = np.sort(x_ana)
        roark_fine = roark_distributed_load_response(x_ana_sorted, L, E, I_z, w_dist, roark_type)
        uy_ana = np.interp(x_ana, x_ana_sorted, roark_fine["deflection"] * 1000)
        th_ana = np.interp(x_ana, x_ana_sorted, np.degrees(roark_fine["rotation"]))

        axes_dist[row, 0].plot(x, uy_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_dist[row, 0].plot(x_ana, uy_ana, "k--", label="Roark")
        axes_dist[row, 0].set_ylabel(r"$u_y$ [mm]")
        axes_dist[row, 0].set_title(title, fontweight="bold")
        axes_dist[row, 0].grid(ls="--", alpha=0.6)
        axes_dist[row, 0].legend(loc="upper right", fontsize="small")

        axes_dist[row, 1].plot(x, th_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_dist[row, 1].plot(x_ana, th_ana, "k--", label="Roark")
        axes_dist[row, 1].set_ylabel(r"$\theta_z$ [deg]")
        axes_dist[row, 1].set_title(title, fontweight="bold")
        axes_dist[row, 1].grid(ls="--", alpha=0.6)
        axes_dist[row, 1].legend(loc="upper right", fontsize="small")

        err_uy = uy_fem - uy_roark
        err_th = th_fem - th_roark
        for i in range(len(x)):
            csv_rows.append([
                float(job_id), x[i], uy_fem[i], uy_roark[i], err_uy[i],
                th_fem[i], th_roark[i], err_th[i],
            ])
        max_uy = np.max(np.abs(err_uy))
        max_th = np.max(np.abs(err_th))
        print(f"  job_{job_id:04d} ({title}): u_y max|err|={max_uy:.6f} mm; theta_z max|err|={max_th:.6f} deg")

    axes_dist[-1, 0].set_xlabel(r"$x$ [m]")
    axes_dist[-1, 1].set_xlabel(r"$x$ [m]")
    fig_dist.tight_layout()
    fig_dist.subplots_adjust(top=0.92)
    fig_dist.savefig(out_dir / "roark_verification_distributed_loads.png", dpi=300)
    plt.close(fig_dist)
    print(f"Saved: {out_dir / 'roark_verification_distributed_loads.png'}")

    if csv_rows:
        csv_path = out_dir / "roark_verification_data.csv"
        header = "job_id,x,uy_fem_mm,uy_roark_mm,error_uy_mm,theta_z_fem_deg,theta_z_roark_deg,error_theta_deg"
        np.savetxt(csv_path, csv_rows, delimiter=",", header=header, comments="")
        print(f"Saved: {csv_path}")
    print("Roark verification done.")


if __name__ == "__main__":
    run_roark_verification()
