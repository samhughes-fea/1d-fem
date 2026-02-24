# post_processing/verification_visualisers/roarks_formulas/roark_section_forces_verification.py
"""
FEM vs Roark section forces (V, M) verification.

Loads FEM section forces from tertiary_results (nodal or gaussian) and compares
Vy (shear) and Mz (bending moment) to Roark's V(x), M(x) for the same jobs as
roark_verification.py: point job_0000–0002, distributed job_0005–0007.
Output: overlay plots and CSV (FEM vs Roark, errors).
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
from pre_processing.parsing.element_parser import ElementParser  # type: ignore

from roarks_formulas_point import roark_point_load_response  # type: ignore
from roarks_formulas_distributed import roark_distributed_load_response  # type: ignore

# Beam and load parameters (match roark_verification.py)
E: Final[float] = 2.0e11   # Pa
I_z: Final[float] = 2.08769e-06  # m^4
P_point: Final[float] = -500.0  # N
w_dist: Final[float] = 500.0  # N/m

# Column order: N, Vy, Vz, T, My, Mz
IDX_VY: Final[int] = 1
IDX_MZ: Final[int] = 5

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

_ANALYTICAL_GRID_FACTOR: Final[int] = 25
_ANALYTICAL_GRID_MIN: Final[int] = 500


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


def _analytical_grid(n_fem: int, L: float) -> np.ndarray:
    n = max(n_fem * _ANALYTICAL_GRID_FACTOR, _ANALYTICAL_GRID_MIN)
    return np.linspace(0.0, L, int(n))


def _read_nodal_section_forces(csv_path: Path) -> Optional[np.ndarray]:
    """Read nodal section forces [N, Vy, Vz, T, My, Mz]. Returns (n_nodes, 6) or None."""
    try:
        with open(csv_path, encoding="utf-8") as f:
            first_line = f.readline()
        skip = 2 if "column_order=resultant" in first_line else 1
        data = np.genfromtxt(csv_path, delimiter=",", skip_header=skip)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] != 6:
            return None
        return data
    except Exception as exc:
        print(f"Error reading {csv_path}: {exc}")
        return None


def _gather_section_forces_from_gaussian(
    job_dir: Path,
    element_dictionary: dict,
    grid_dictionary: dict,
    n_nodes: int,
) -> Optional[tuple[np.ndarray, np.ndarray]]:
    """Build (x, forces) at nodes from gaussian section_forces_elem_*.csv. forces shape (n_nodes, 6)."""
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
    ids = element_dictionary.get("ids", np.arange(len(connectivity)))
    forces_nodal = np.zeros((n_nodes, 6))
    weight = np.zeros(n_nodes)
    for i, csv_path in enumerate(files):
        try:
            with open(csv_path, encoding="utf-8") as f:
                first_line = f.readline()
                skip = 2 if "column_order=resultant" in first_line else 1
                if skip == 2:
                    second = f.readline()
                    if second.strip().startswith("# xi_per_row="):
                        skip = 3
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
    nonzero = weight > 0
    if not np.any(nonzero):
        return None
    forces_nodal[nonzero] /= weight[nonzero, np.newaxis]
    return (x_all, forces_nodal)


def run_section_forces_verification() -> None:
    results_dir: Path = PROJECT_ROOT / "post_processing" / "results"
    jobs_dir: Path = PROJECT_ROOT / "jobs"
    out_dir: Path = SCRIPT_DIR / "verification"
    out_dir.mkdir(exist_ok=True)

    pattern = str(results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
    csv_files = sorted(glob.glob(pattern))
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
        _, csv_file, job_dir, job_input_name = candidates[0]
        job_to_result[base_id] = (csv_file, job_dir, job_input_name)

    if not job_to_result:
        print("No FEM result directories found. Run jobs first.")
        return

    csv_rows: list[list[float]] = []

    # ----- Point loads: V and M (FEM vs Roark) -----
    fig_point, axes_point = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
    fig_point.suptitle(
        "Roark section forces verification: point loads (Euler–Bernoulli)",
        fontsize=14,
        fontweight="bold",
    )

    for row, (job_id, title, a_frac) in enumerate(POINT_LOAD_CASES):
        if job_id not in job_to_result:
            print(f"WARNING: No result for job_{job_id:04d}, skipping '{title}'.")
            continue
        _, job_dir, job_input_name = job_to_result[job_id]
        grid_file = jobs_dir / job_input_name / "grid.txt"
        if not grid_file.is_file():
            print(f"WARNING: Missing {grid_file}, skipping job_{job_id:04d}.")
            continue
        grid = GridParser(str(grid_file), str(jobs_dir / job_input_name)).parse()
        try:
            node_coords = _get_node_coordinates(grid)
        except Exception as exc:
            print(f"WARNING: {exc} job_{job_id:04d}")
            continue
        x = node_coords[:, 0]
        n_nodes = x.shape[0]
        L = float(np.max(x))
        a = a_frac * L
        load_type = "end" if a_frac >= 1.0 else ("mid" if a_frac >= 0.5 else "quarter")

        nodal_csv = job_dir / "tertiary_results" / "nodal" / "nodal_section_forces.csv"
        if nodal_csv.is_file():
            forces = _read_nodal_section_forces(nodal_csv)
            if forces is None or forces.shape[0] != n_nodes:
                forces = None
            else:
                x_fem, forces_fem = x, forces
        else:
            forces = None
        if forces is None:
            element_file = jobs_dir / job_input_name / "element.txt"
            if element_file.is_file():
                elem_parsed = ElementParser(str(element_file), str(jobs_dir / job_input_name)).parse()
                grid_dict = grid["grid_dictionary"] if isinstance(grid, dict) else {}
                out = _gather_section_forces_from_gaussian(
                    job_dir,
                    elem_parsed["element_dictionary"],
                    grid_dict,
                    n_nodes,
                )
                if out is not None:
                    x_fem, forces_fem = out
                else:
                    print(f"WARNING: No section forces for job_{job_id:04d}, skipping.")
                    continue
            else:
                print(f"WARNING: No section forces for job_{job_id:04d}, skipping.")
                continue
        else:
            x_fem, forces_fem = x, forces

        V_fem = forces_fem[:, IDX_VY]
        M_fem = forces_fem[:, IDX_MZ]

        roark_at_nodes = roark_point_load_response(x_fem, L, E, I_z, P_point, load_type)
        V_roark = roark_at_nodes["shear"]
        M_roark = roark_at_nodes["moment"]

        x_ana = _analytical_grid(x_fem.shape[0], L)
        roark_fine = roark_point_load_response(x_ana, L, E, I_z, P_point, load_type)
        V_ana = roark_fine["shear"]
        M_ana = roark_fine["moment"]

        axes_point[row, 0].plot(x_fem, V_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_point[row, 0].plot(x_ana, V_ana, "k--", label="Roark")
        axes_point[row, 0].set_ylabel(r"$V_y$ [N]")
        axes_point[row, 0].set_title(title, fontweight="bold")
        axes_point[row, 0].grid(ls="--", alpha=0.6)
        axes_point[row, 0].legend(loc="upper right", fontsize="small")

        axes_point[row, 1].plot(x_fem, M_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_point[row, 1].plot(x_ana, M_ana, "k--", label="Roark")
        axes_point[row, 1].set_ylabel(r"$M_z$ [N·m]")
        axes_point[row, 1].set_title(title, fontweight="bold")
        axes_point[row, 1].grid(ls="--", alpha=0.6)
        axes_point[row, 1].legend(loc="upper right", fontsize="small")

        err_V = V_fem - V_roark
        err_M = M_fem - M_roark
        for i in range(len(x_fem)):
            csv_rows.append([
                float(job_id), x_fem[i], V_fem[i], V_roark[i], err_V[i],
                M_fem[i], M_roark[i], err_M[i],
            ])
        print(f"  job_{job_id:04d} ({title}): V_y max|err|={np.max(np.abs(err_V)):.4f} N; "
              f"M_z max|err|={np.max(np.abs(err_M)):.4f} N·m")

    axes_point[-1, 0].set_xlabel(r"$x$ [m]")
    axes_point[-1, 1].set_xlabel(r"$x$ [m]")
    fig_point.tight_layout()
    fig_point.subplots_adjust(top=0.92)
    fig_point.savefig(out_dir / "roark_section_forces_point_loads.png", dpi=300)
    plt.close(fig_point)
    print(f"Saved: {out_dir / 'roark_section_forces_point_loads.png'}")

    # ----- Distributed loads: V and M -----
    fig_dist, axes_dist = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
    fig_dist.suptitle(
        "Roark section forces verification: distributed loads (Euler–Bernoulli)",
        fontsize=14,
        fontweight="bold",
    )

    for row, (job_id, title, roark_type) in enumerate(DIST_LOAD_CASES):
        if job_id not in job_to_result:
            print(f"WARNING: No result for job_{job_id:04d}, skipping '{title}'.")
            continue
        _, job_dir, job_input_name = job_to_result[job_id]
        grid_file = jobs_dir / job_input_name / "grid.txt"
        if not grid_file.is_file():
            continue
        grid = GridParser(str(grid_file), str(jobs_dir / job_input_name)).parse()
        try:
            node_coords = _get_node_coordinates(grid)
        except Exception:
            continue
        x = node_coords[:, 0]
        n_nodes = x.shape[0]
        L = float(np.max(x))

        nodal_csv = job_dir / "tertiary_results" / "nodal" / "nodal_section_forces.csv"
        if nodal_csv.is_file():
            forces = _read_nodal_section_forces(nodal_csv)
            if forces is not None and forces.shape[0] == n_nodes:
                x_fem, forces_fem = x, forces
            else:
                forces = None
        else:
            forces = None
        if forces is None:
            element_file = jobs_dir / job_input_name / "element.txt"
            if element_file.is_file():
                elem_parsed = ElementParser(str(element_file), str(jobs_dir / job_input_name)).parse()
                grid_dict = grid["grid_dictionary"] if isinstance(grid, dict) else {}
                out = _gather_section_forces_from_gaussian(
                    job_dir,
                    elem_parsed["element_dictionary"],
                    grid_dict,
                    n_nodes,
                )
                if out is not None:
                    x_fem, forces_fem = out
                else:
                    continue
            else:
                continue
        else:
            x_fem, forces_fem = x, forces

        V_fem = forces_fem[:, IDX_VY]
        M_fem = forces_fem[:, IDX_MZ]

        order = np.argsort(x_fem)
        x_sorted = x_fem[order]
        roark_at_nodes = roark_distributed_load_response(x_sorted, L, E, I_z, w_dist, roark_type)
        V_roark = np.interp(x_fem, x_sorted, roark_at_nodes["shear"])
        M_roark = np.interp(x_fem, x_sorted, roark_at_nodes["moment"])

        x_ana = _analytical_grid(x_fem.shape[0], L)
        x_ana_sorted = np.sort(x_ana)
        roark_fine = roark_distributed_load_response(x_ana_sorted, L, E, I_z, w_dist, roark_type)
        V_ana = np.interp(x_ana, x_ana_sorted, roark_fine["shear"])
        M_ana = np.interp(x_ana, x_ana_sorted, roark_fine["moment"])

        axes_dist[row, 0].plot(x_fem, V_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_dist[row, 0].plot(x_ana, V_ana, "k--", label="Roark")
        axes_dist[row, 0].set_ylabel(r"$V_y$ [N]")
        axes_dist[row, 0].set_title(title, fontweight="bold")
        axes_dist[row, 0].grid(ls="--", alpha=0.6)
        axes_dist[row, 0].legend(loc="upper right", fontsize="small")

        axes_dist[row, 1].plot(x_fem, M_fem, color=COLORS[row], linestyle="-", label=f"FEM job_{job_id:04d}")
        axes_dist[row, 1].plot(x_ana, M_ana, "k--", label="Roark")
        axes_dist[row, 1].set_ylabel(r"$M_z$ [N·m]")
        axes_dist[row, 1].set_title(title, fontweight="bold")
        axes_dist[row, 1].grid(ls="--", alpha=0.6)
        axes_dist[row, 1].legend(loc="upper right", fontsize="small")

        err_V = V_fem - V_roark
        err_M = M_fem - M_roark
        for i in range(len(x_fem)):
            csv_rows.append([
                float(job_id), x_fem[i], V_fem[i], V_roark[i], err_V[i],
                M_fem[i], M_roark[i], err_M[i],
            ])
        print(f"  job_{job_id:04d} ({title}): V_y max|err|={np.max(np.abs(err_V)):.4f} N; "
              f"M_z max|err|={np.max(np.abs(err_M)):.4f} N·m")

    axes_dist[-1, 0].set_xlabel(r"$x$ [m]")
    axes_dist[-1, 1].set_xlabel(r"$x$ [m]")
    fig_dist.tight_layout()
    fig_dist.subplots_adjust(top=0.92)
    fig_dist.savefig(out_dir / "roark_section_forces_distributed_loads.png", dpi=300)
    plt.close(fig_dist)
    print(f"Saved: {out_dir / 'roark_section_forces_distributed_loads.png'}")

    if csv_rows:
        csv_path = out_dir / "roark_section_forces_verification_data.csv"
        header = "job_id,x,V_fem_N,V_roark_N,error_V_N,M_fem_Nm,M_roark_Nm,error_M_Nm"
        np.savetxt(csv_path, csv_rows, delimiter=",", header=header, comments="")
        print(f"Saved: {csv_path}")
    print("Roark section forces verification done.")


if __name__ == "__main__":
    run_section_forces_verification()
