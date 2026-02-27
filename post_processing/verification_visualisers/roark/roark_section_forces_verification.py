# post_processing/verification_visualisers/roark_section_forces_verification.py
"""
FEM vs Roark section forces (V, M) verification.

Loads FEM section forces from tertiary_results (nodal or gaussian) and compares
Vy (shear) and Mz (bending moment) to Roark's V(x), M(x) for the same jobs as
roark_verification.py: point job_0000–0002, distributed job_0003–0005.
Output: overlay plots and CSV (FEM vs Roark, errors).

Sign convention (SFD and BMD)
----------------------------
We use the standard structural convention:

- **Positive shear** → tends to rotate a small element clockwise.
- **Positive bending moment** → causes sagging (bottom fibre in tension).
  **Negative bending moment** → hogging (top fibre in tension).

For a cantilever fixed at the left (x=0) and free at x=L with **downward** transverse
load (positive q(x) or F_y negative in global y-up):
- **SFD (V_y)**: negative; magnitude grows from tip toward the fixed end.
- **BMD (M_z)**: zero at the tip, negative (hogging) toward the fixed end; maximum
  magnitude at the support. Same applies for point loads, UDL, triangular, or parabolic q(x).

Roark point-load formulas (V = -P, M = -P*(a-x)) use the opposite sign for our load
(F_y = -500). We apply ROARK_POINT_LOAD_SIGN = -1 to Roark V and M for point-load
plots and CSV so they match FEM and the convention above.
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

try:
    from labellines import labelLines
except ImportError:
    labelLines = None  # optional: pip install matplotlib-label-lines

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = next(
    (p for p in SCRIPT_DIR.parents if (p / "pre_processing").is_dir()),
    SCRIPT_DIR.parents[4],
)
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR / "roark_utilities"))

from pre_processing.parsing.grid_parser import GridParser  # type: ignore
from pre_processing.parsing.element_parser import ElementParser  # type: ignore

from roarks_formulas_euler_bernoulli_point import roark_point_load_response  # type: ignore
from roarks_formulas_euler_bernoulli_distributed import roark_distributed_load_response  # type: ignore

# Beam and load parameters (match roark_verification.py)
E: Final[float] = 2.0e11   # Pa
I_z: Final[float] = 2.08769e-06  # m^4
P_point: Final[float] = -500.0  # N (matches job F_y: positive transverse load = downward in global y)
w_dist: Final[float] = 500.0  # N/m

# Sign convention: for positive transverse (downward) load, structural convention gives V_y < 0 and M_z < 0
# (shear and hogging moment). Roark point-load formulas use V = -P, M = -P*(a-x); with P = F_y = -500
# they yield V, M positive. We negate Roark V and M for point-load plots/CSV so they match FEM.
ROARK_POINT_LOAD_SIGN: Final[float] = -1.0  # apply to Roark V, M for point loads only
ROARK_FILL_ALPHA: Final[float] = 0.25
ROARK_FILL_COLOR: Final[str] = "#808080"  # neutral gray; fill between FEM curve and zero axis

# Column order: N, Vy, Vz, T, My, Mz
IDX_VY: Final[int] = 1
IDX_MZ: Final[int] = 5

POINT_LOAD_CASES: Final[list[tuple[int, str, float]]] = [
    (0, "End", 1.0),
    (1, "Mid-span", 0.5),
    (2, "Quarter-span", 0.25),
]
DIST_LOAD_CASES: Final[list[tuple[int, str, str]]] = [
    (3, "UDL", "udl"),
    (4, "Triangular", "triangular"),
    (5, "Parabolic", "parabolic"),
]
# Timoshenko jobs 6–11 (V, M same as E–B Roark)
POINT_LOAD_CASES_TIMS: Final[list[tuple[int, str, str]]] = [
    (6, "End", "end"),
    (7, "Mid-span", "mid"),
    (8, "Quarter-span", "quarter"),
]
DIST_LOAD_CASES_TIMS: Final[list[tuple[int, str, str]]] = [
    (9, "UDL", "udl"),
    (10, "Triangular", "triangular"),
    (11, "Parabolic", "parabolic"),
]

# Analytical V_y(x) and M_z(x) for legend and labelLines. Key: (job_id, col) with col 0 = V_y, 1 = M_z.
FORMULA_V_M_POINT: Final[dict[tuple[int, int], str]] = {
    (0, 0): r"$V_y = -P$",
    (0, 1): r"$M_z = -P(L-x)$",
    (1, 0): r"$V_y = -P$, $x<a$",
    (1, 1): r"$M_z = -P(a-x)$, $x<a$",
    (2, 0): r"$V_y = -P$, $x<a$",
    (2, 1): r"$M_z = -P(a-x)$, $x<a$",
    (6, 0): r"$V_y = -P$",
    (6, 1): r"$M_z = -P(L-x)$",
    (7, 0): r"$V_y = -P$, $x<a$",
    (7, 1): r"$M_z = -P(a-x)$, $x<a$",
    (8, 0): r"$V_y = -P$, $x<a$",
    (8, 1): r"$M_z = -P(a-x)$, $x<a$",
}
FORMULA_V_M_DIST: Final[dict[tuple[int, int], str]] = {
    (3, 0): r"$V_y = -w(L-x)$",
    (3, 1): r"$M_z = -\frac{w}{2}(L-x)^2$",
    (4, 0): r"$V_y = -\frac{w(L^2-x^2)}{2L}$",
    (4, 1): r"$M_z = -\frac{w(L-x)^2(2L+x)}{6L}$",
    (5, 0): r"$V_y = -\frac{w(L^3-x^3)}{3L^2}$",
    (5, 1): r"$M_z = -\frac{w(3L^4-4L^3x+x^4)}{12L^2}$",
    (9, 0): r"$V_y = -w(L-x)$",
    (9, 1): r"$M_z = -\frac{w}{2}(L-x)^2$",
    (10, 0): r"$V_y = -\frac{w(L^2-x^2)}{2L}$",
    (10, 1): r"$M_z = -\frac{w(L-x)^2(2L+x)}{6L}$",
    (11, 0): r"$V_y = -\frac{w(L^3-x^3)}{3L^2}$",
    (11, 1): r"$M_z = -\frac{w(3L^4-4L^3x+x^4)}{12L^2}$",
}

COLORS: Final[list[str]] = ["#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#4BACC6"]
COLORS_6: Final[list[str]] = ["#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#4BACC6", "#F79646"]


def _legend_num_anal(handles: list, labels: list) -> tuple[list, list]:
    """Interleave (FEM, Roark) so each legend row is (num_i, anal_i); for 1×2 all-loads figures."""
    num = [(h, l) for h, l in zip(handles, labels) if l.startswith("job_")]
    anal = [(h, l) for h, l in zip(handles, labels) if not l.startswith("job_")]
    h_inter = [h for p in zip([n[0] for n in num], [a[0] for a in anal]) for h in p]
    l_inter = [l for p in zip([n[1] for n in num], [a[1] for a in anal]) for l in p]
    return h_inter, l_inter


def _section_forces_col_legend(axes_2col: np.ndarray, col: int, nrows: int = 3) -> tuple[list, list]:
    """Build interleaved (FEM, Roark) handles/labels for a single column of the 3×2 grid.

    With ``ncol=2`` in the legend the left column shows FEM entries and
    the right column shows analytical formulas, matching the deformation
    convergence plot layout.
    """
    handles, labels = [], []
    for r in range(nrows):
        h, l = axes_2col[r, col].get_legend_handles_labels()
        num = [(hi, li) for hi, li in zip(h, l) if "FEM" in li]
        anal = [(hi, li) for hi, li in zip(h, l) if "FEM" not in li]
        if num and anal:
            handles.extend([num[0][0], anal[0][0]])
            labels.extend([num[0][1], anal[0][1]])
    return handles, labels


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
    out_dir: Path = SCRIPT_DIR / "section_forces_plots"
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
    fig_point, axes_point = plt.subplots(3, 2, figsize=(14, 15), sharex=True)
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
        V_roark = ROARK_POINT_LOAD_SIGN * roark_at_nodes["shear"]
        M_roark = ROARK_POINT_LOAD_SIGN * roark_at_nodes["moment"]

        x_ana = _analytical_grid(x_fem.shape[0], L)
        roark_fine = roark_point_load_response(x_ana, L, E, I_z, P_point, load_type)
        V_ana = ROARK_POINT_LOAD_SIGN * roark_fine["shear"]
        M_ana = ROARK_POINT_LOAD_SIGN * roark_fine["moment"]

        formula_v = FORMULA_V_M_POINT.get((job_id, 0), "Roark")
        formula_m = FORMULA_V_M_POINT.get((job_id, 1), "Roark")
        axes_point[row, 0].fill_between(x_fem, V_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
        axes_point[row, 0].plot(x_fem, V_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
        (line_v,) = axes_point[row, 0].plot(x_ana, V_ana, "k-", linewidth=0.5, label=formula_v)
        axes_point[row, 0].set_ylabel(r"$V_y$ [N]")
        axes_point[row, 0].set_title(title, fontweight="bold")
        axes_point[row, 0].grid(ls="--", alpha=0.6)
        if labelLines is not None:
            y_mid_v = (np.min(V_ana) + np.max(V_ana)) / 2
            # Point-load shear is a step; place label on the jump segment at x-mid of jump so y ≈ (V_min+V_max)/2
            if len(np.unique(V_ana)) == 2:
                idx_r = np.searchsorted(x_ana, a, side="left")
                idx_r = min(idx_r, len(x_ana) - 1)
                idx_l = max(0, idx_r - 1)
                x_mid_v = float((x_ana[idx_l] + x_ana[idx_r]) / 2)
            else:
                x_mid_v = float(x_ana[np.argmin(np.abs(V_ana - y_mid_v))])
            labelLines([line_v], align=True, fontsize="small", xvals=[x_mid_v])

        axes_point[row, 1].fill_between(x_fem, M_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
        axes_point[row, 1].plot(x_fem, M_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
        (line_m,) = axes_point[row, 1].plot(x_ana, M_ana, "k-", linewidth=0.5, label=formula_m)
        axes_point[row, 1].set_ylabel(r"$M_z$ [N·m]")
        axes_point[row, 1].set_title(title, fontweight="bold")
        axes_point[row, 1].grid(ls="--", alpha=0.6)
        if labelLines is not None:
            y_mid_m = (np.min(M_ana) + np.max(M_ana)) / 2
            x_mid_m = float(x_ana[np.argmin(np.abs(M_ana - y_mid_m))])
            labelLines([line_m], align=True, fontsize="small", xvals=[x_mid_m])

        err_V = V_fem - V_roark
        err_M = M_fem - M_roark
        for i in range(len(x_fem)):
            csv_rows.append([
                float(job_id), x_fem[i], V_fem[i], V_roark[i], err_V[i],
                M_fem[i], M_roark[i], err_M[i],
            ])
        print(f"  job_{job_id:04d} ({title}): V_y max|err|={np.max(np.abs(err_V)):.4f} N; "
              f"M_z max|err|={np.max(np.abs(err_M)):.4f} N·m")

    for c in range(2):
        axes_point[2, c].set_xlabel(r"$x$ [m]")
    h_v, l_v = _section_forces_col_legend(axes_point, col=0)
    h_m, l_m = _section_forces_col_legend(axes_point, col=1)
    fig_point.tight_layout()
    fig_point.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.28, wspace=0.22, hspace=0.46)
    if h_v:
        axes_point[2, 0].legend(h_v, l_v, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
    if h_m:
        axes_point[2, 1].legend(h_m, l_m, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
    fig_point.savefig(out_dir / "roark_section_forces_point_loads_euler_bernoulli.png", dpi=300, bbox_inches="tight")
    plt.close(fig_point)
    print(f"Saved: {out_dir / 'roark_section_forces_point_loads_euler_bernoulli.png'}")

    # ----- Distributed loads: V and M -----
    fig_dist, axes_dist = plt.subplots(3, 2, figsize=(14, 15), sharex=True)
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

        formula_v = FORMULA_V_M_DIST.get((job_id, 0), "Roark")
        formula_m = FORMULA_V_M_DIST.get((job_id, 1), "Roark")
        axes_dist[row, 0].fill_between(x_fem, V_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
        axes_dist[row, 0].plot(x_fem, V_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
        (line_v,) = axes_dist[row, 0].plot(x_ana, V_ana, "k-", linewidth=0.5, label=formula_v)
        axes_dist[row, 0].set_ylabel(r"$V_y$ [N]")
        axes_dist[row, 0].set_title(title, fontweight="bold")
        axes_dist[row, 0].grid(ls="--", alpha=0.6)
        if labelLines is not None:
            y_mid_v = (np.min(V_ana) + np.max(V_ana)) / 2
            x_mid_v = float(x_ana[np.argmin(np.abs(V_ana - y_mid_v))])
            labelLines([line_v], align=True, fontsize="x-small", xvals=[x_mid_v])

        axes_dist[row, 1].fill_between(x_fem, M_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
        axes_dist[row, 1].plot(x_fem, M_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
        (line_m,) = axes_dist[row, 1].plot(x_ana, M_ana, "k-", linewidth=0.5, label=formula_m)
        axes_dist[row, 1].set_ylabel(r"$M_z$ [N·m]")
        axes_dist[row, 1].set_title(title, fontweight="bold")
        axes_dist[row, 1].grid(ls="--", alpha=0.6)
        if labelLines is not None:
            y_mid_m = (np.min(M_ana) + np.max(M_ana)) / 2
            x_mid_m = float(x_ana[np.argmin(np.abs(M_ana - y_mid_m))])
            labelLines([line_m], align=True, fontsize="x-small", xvals=[x_mid_m])

        err_V = V_fem - V_roark
        err_M = M_fem - M_roark
        for i in range(len(x_fem)):
            csv_rows.append([
                float(job_id), x_fem[i], V_fem[i], V_roark[i], err_V[i],
                M_fem[i], M_roark[i], err_M[i],
            ])
        print(f"  job_{job_id:04d} ({title}): V_y max|err|={np.max(np.abs(err_V)):.4f} N; "
              f"M_z max|err|={np.max(np.abs(err_M)):.4f} N·m")

    for c in range(2):
        axes_dist[2, c].set_xlabel(r"$x$ [m]")
    h_v, l_v = _section_forces_col_legend(axes_dist, col=0)
    h_m, l_m = _section_forces_col_legend(axes_dist, col=1)
    fig_dist.tight_layout()
    fig_dist.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.28, wspace=0.22, hspace=0.46)
    if h_v:
        axes_dist[2, 0].legend(h_v, l_v, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
    if h_m:
        axes_dist[2, 1].legend(h_m, l_m, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
    fig_dist.savefig(out_dir / "roark_section_forces_distributed_loads_euler_bernoulli.png", dpi=300, bbox_inches="tight")
    plt.close(fig_dist)
    print(f"Saved: {out_dir / 'roark_section_forces_distributed_loads_euler_bernoulli.png'}")

    # ----- Timoshenko jobs 6–11: section forces (V, M same as Roark E–B) -----
    if any(job_id in job_to_result for job_id, _, _ in POINT_LOAD_CASES_TIMS + DIST_LOAD_CASES_TIMS):
        fig_pt, axes_pt = plt.subplots(3, 2, figsize=(14, 15), sharex=True)
        fig_pt.suptitle("Roark section forces: point loads (Timoshenko jobs 6–8)", fontsize=14, fontweight="bold")
        for row, (job_id, title, load_type) in enumerate(POINT_LOAD_CASES_TIMS):
            if job_id not in job_to_result:
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
                if forces is None or forces.shape[0] != n_nodes:
                    forces = None
                else:
                    x_fem, forces_fem = x, forces
            else:
                forces = None
            if forces is None:
                element_file = jobs_dir / job_input_name / "element.txt"
                if not element_file.is_file():
                    continue
                elem_parsed = ElementParser(str(element_file), str(jobs_dir / job_input_name)).parse()
                grid_dict = grid["grid_dictionary"] if isinstance(grid, dict) else {}
                out = _gather_section_forces_from_gaussian(job_dir, elem_parsed["element_dictionary"], grid_dict, n_nodes)
                if out is None:
                    continue
                x_fem, forces_fem = out
            else:
                x_fem, forces_fem = x, forces
            V_fem = forces_fem[:, IDX_VY]
            M_fem = forces_fem[:, IDX_MZ]
            roark_at_nodes = roark_point_load_response(x_fem, L, E, I_z, P_point, load_type)
            x_ana = _analytical_grid(x_fem.shape[0], L)
            roark_fine = roark_point_load_response(x_ana, L, E, I_z, P_point, load_type)
            V_ana_pt = ROARK_POINT_LOAD_SIGN * roark_fine["shear"]
            M_ana_pt = ROARK_POINT_LOAD_SIGN * roark_fine["moment"]
            formula_v = FORMULA_V_M_POINT.get((job_id, 0), "Roark")
            formula_m = FORMULA_V_M_POINT.get((job_id, 1), "Roark")
            axes_pt[row, 0].fill_between(x_fem, V_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
            axes_pt[row, 0].plot(x_fem, V_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
            (line_v,) = axes_pt[row, 0].plot(x_ana, V_ana_pt, "k-", linewidth=0.5, label=formula_v)
            axes_pt[row, 0].set_ylabel(r"$V_y$ [N]")
            axes_pt[row, 0].set_title(title, fontweight="bold")
            axes_pt[row, 0].grid(ls="--", alpha=0.6)
            if labelLines is not None:
                y_mid_v = (np.min(V_ana_pt) + np.max(V_ana_pt)) / 2
                a_pt = L if load_type == "end" else (L / 2 if load_type == "mid" else L / 4)
                if len(np.unique(V_ana_pt)) == 2:
                    idx_r = np.searchsorted(x_ana, a_pt, side="left")
                    idx_r = min(idx_r, len(x_ana) - 1)
                    idx_l = max(0, idx_r - 1)
                    x_mid_v = float((x_ana[idx_l] + x_ana[idx_r]) / 2)
                else:
                    x_mid_v = float(x_ana[np.argmin(np.abs(V_ana_pt - y_mid_v))])
                labelLines([line_v], align=True, fontsize="small", xvals=[x_mid_v])
            axes_pt[row, 1].fill_between(x_fem, M_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
            axes_pt[row, 1].plot(x_fem, M_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
            (line_m,) = axes_pt[row, 1].plot(x_ana, M_ana_pt, "k-", linewidth=0.5, label=formula_m)
            axes_pt[row, 1].set_ylabel(r"$M_z$ [N·m]")
            axes_pt[row, 1].set_title(title, fontweight="bold")
            axes_pt[row, 1].grid(ls="--", alpha=0.6)
            if labelLines is not None:
                y_mid_m = (np.min(M_ana_pt) + np.max(M_ana_pt)) / 2
                x_mid_m = float(x_ana[np.argmin(np.abs(M_ana_pt - y_mid_m))])
                labelLines([line_m], align=True, fontsize="small", xvals=[x_mid_m])
        for c in range(2):
            axes_pt[2, c].set_xlabel(r"$x$ [m]")
        h_v, l_v = _section_forces_col_legend(axes_pt, col=0)
        h_m, l_m = _section_forces_col_legend(axes_pt, col=1)
        fig_pt.tight_layout()
        fig_pt.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.28, wspace=0.22, hspace=0.46)
        if h_v:
            axes_pt[2, 0].legend(h_v, l_v, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        if h_m:
            axes_pt[2, 1].legend(h_m, l_m, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        fig_pt.savefig(out_dir / "roark_section_forces_point_loads_timoshenko.png", dpi=300, bbox_inches="tight")
        plt.close(fig_pt)
        print(f"Saved: {out_dir / 'roark_section_forces_point_loads_timoshenko.png'}")

        fig_dt, axes_dt = plt.subplots(3, 2, figsize=(14, 15), sharex=True)
        fig_dt.suptitle("Roark section forces: distributed loads (Timoshenko jobs 9–11)", fontsize=14, fontweight="bold")
        for row, (job_id, title, roark_type) in enumerate(DIST_LOAD_CASES_TIMS):
            if job_id not in job_to_result:
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
                if not element_file.is_file():
                    continue
                elem_parsed = ElementParser(str(element_file), str(jobs_dir / job_input_name)).parse()
                grid_dict = grid["grid_dictionary"] if isinstance(grid, dict) else {}
                out = _gather_section_forces_from_gaussian(job_dir, elem_parsed["element_dictionary"], grid_dict, n_nodes)
                if out is None:
                    continue
                x_fem, forces_fem = out
            else:
                x_fem, forces_fem = x, forces
            V_fem = forces_fem[:, IDX_VY]
            M_fem = forces_fem[:, IDX_MZ]
            order = np.argsort(x_fem)
            x_sorted = x_fem[order]
            roark_at_nodes = roark_distributed_load_response(x_sorted, L, E, I_z, w_dist, roark_type)
            x_ana = _analytical_grid(x_fem.shape[0], L)
            x_ana_sorted = np.sort(x_ana)
            roark_fine = roark_distributed_load_response(x_ana_sorted, L, E, I_z, w_dist, roark_type)
            formula_v = FORMULA_V_M_DIST.get((job_id, 0), "Roark")
            formula_m = FORMULA_V_M_DIST.get((job_id, 1), "Roark")
            V_ana_dt = np.interp(x_ana, x_ana_sorted, roark_fine["shear"])
            M_ana_dt = np.interp(x_ana, x_ana_sorted, roark_fine["moment"])
            axes_dt[row, 0].fill_between(x_fem, V_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
            axes_dt[row, 0].plot(x_fem, V_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
            (line_v,) = axes_dt[row, 0].plot(x_ana, V_ana_dt, "k-", linewidth=0.5, label=formula_v)
            axes_dt[row, 0].set_ylabel(r"$V_y$ [N]")
            axes_dt[row, 0].set_title(title, fontweight="bold")
            axes_dt[row, 0].grid(ls="--", alpha=0.6)
            if labelLines is not None:
                y_mid_v = (np.min(V_ana_dt) + np.max(V_ana_dt)) / 2
                x_mid_v = float(x_ana[np.argmin(np.abs(V_ana_dt - y_mid_v))])
                labelLines([line_v], align=True, fontsize="x-small", xvals=[x_mid_v])
            axes_dt[row, 1].fill_between(x_fem, M_fem, 0, color=ROARK_FILL_COLOR, alpha=ROARK_FILL_ALPHA)
            axes_dt[row, 1].plot(x_fem, M_fem, color=COLORS[row], linestyle="--", linewidth=1.25, label=f"FEM job_{job_id:04d}")
            (line_m,) = axes_dt[row, 1].plot(x_ana, M_ana_dt, "k-", linewidth=0.5, label=formula_m)
            axes_dt[row, 1].set_ylabel(r"$M_z$ [N·m]")
            axes_dt[row, 1].set_title(title, fontweight="bold")
            axes_dt[row, 1].grid(ls="--", alpha=0.6)
            if labelLines is not None:
                y_mid_m = (np.min(M_ana_dt) + np.max(M_ana_dt)) / 2
                x_mid_m = float(x_ana[np.argmin(np.abs(M_ana_dt - y_mid_m))])
                labelLines([line_m], align=True, fontsize="x-small", xvals=[x_mid_m])
        for c in range(2):
            axes_dt[2, c].set_xlabel(r"$x$ [m]")
        h_v, l_v = _section_forces_col_legend(axes_dt, col=0)
        h_m, l_m = _section_forces_col_legend(axes_dt, col=1)
        fig_dt.tight_layout()
        fig_dt.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.28, wspace=0.22, hspace=0.46)
        if h_v:
            axes_dt[2, 0].legend(h_v, l_v, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        if h_m:
            axes_dt[2, 1].legend(h_m, l_m, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        fig_dt.savefig(out_dir / "roark_section_forces_distributed_loads_timoshenko.png", dpi=300, bbox_inches="tight")
        plt.close(fig_dt)
        print(f"Saved: {out_dir / 'roark_section_forces_distributed_loads_timoshenko.png'}")

    # ----- 1×2 Euler–Bernoulli: V_y and M_z, all six load cases overlaid -----
    all_load_cases_eb: list[tuple[int, str]] = [
        (0, "End"),
        (1, "Mid-span"),
        (2, "Quarter-span"),
        (3, "UDL"),
        (4, "Triangular"),
        (5, "Parabolic"),
    ]
    if any(job_id in job_to_result for job_id, _ in all_load_cases_eb):
        fig_eb, axes_eb = plt.subplots(1, 2, figsize=(12, 8), sharex=True)
        fig_eb.suptitle(
            r"$V_y$ and $M_z$: all load cases (Euler–Bernoulli)",
            fontsize=14,
            fontweight="bold",
        )
        for col, ylabel in enumerate([r"$V_y$ [N]", r"$M_z$ [N·m]"]):
            ax = axes_eb[col]
            roark_lines: list = []
            for load_idx, (job_id, title) in enumerate(all_load_cases_eb):
                if job_id not in job_to_result:
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
                            continue
                    else:
                        continue
                else:
                    x_fem, forces_fem = x, forces

                V_fem = forces_fem[:, IDX_VY]
                M_fem = forces_fem[:, IDX_MZ]
                if job_id <= 2:
                    a_frac = POINT_LOAD_CASES[job_id][2]
                    load_type = "end" if a_frac >= 1.0 else ("mid" if a_frac >= 0.5 else "quarter")
                    roark_at_nodes = roark_point_load_response(x_fem, L, E, I_z, P_point, load_type)
                    V_roark = ROARK_POINT_LOAD_SIGN * roark_at_nodes["shear"]
                    M_roark = ROARK_POINT_LOAD_SIGN * roark_at_nodes["moment"]
                    x_ana = _analytical_grid(x_fem.shape[0], L)
                    roark_fine = roark_point_load_response(x_ana, L, E, I_z, P_point, load_type)
                    V_ana = ROARK_POINT_LOAD_SIGN * roark_fine["shear"]
                    M_ana = ROARK_POINT_LOAD_SIGN * roark_fine["moment"]
                else:
                    roark_type = DIST_LOAD_CASES[job_id - 3][2]
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

                formula_v = FORMULA_V_M_POINT.get((job_id, 0)) or FORMULA_V_M_DIST.get((job_id, 0), title)
                formula_m = FORMULA_V_M_POINT.get((job_id, 1)) or FORMULA_V_M_DIST.get((job_id, 1), title)
                color = COLORS_6[load_idx % len(COLORS_6)]
                vals_fem = V_fem if col == 0 else M_fem
                vals_ana = V_ana if col == 0 else M_ana
                formula = formula_v if col == 0 else formula_m
                ax.plot(x_fem, vals_fem, color=color, linestyle="--", linewidth=1.25, label=f"job_{job_id:04d} ({title})")
                (line,) = ax.plot(x_ana, vals_ana, "k-", linewidth=0.5, label=formula)
                roark_lines.append(line)
            ax.set_ylabel(ylabel)
            ax.grid(ls="--", alpha=0.6)
            if labelLines is not None and roark_lines:
                L_ax = float(np.max(roark_lines[0].get_xdata()))
                n_lines = len(roark_lines)
                xvals = np.linspace(0.2 * L_ax, 0.85 * L_ax, n_lines).tolist()
                labelLines(roark_lines, align=True, fontsize="x-small", xvals=xvals)
        h_left, l_left = _legend_num_anal(*axes_eb[0].get_legend_handles_labels())
        h_right, l_right = _legend_num_anal(*axes_eb[1].get_legend_handles_labels())
        fig_eb.tight_layout()
        fig_eb.subplots_adjust(top=0.88, bottom=0.32)
        axes_eb[0].legend(h_left, l_left, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        axes_eb[1].legend(h_right, l_right, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        axes_eb[0].set_xlabel(r"$x$ [m]")
        axes_eb[1].set_xlabel(r"$x$ [m]")
        path_eb = out_dir / "roark_section_forces_convergence_all_loads_euler_bernoulli.png"
        fig_eb.savefig(path_eb, dpi=300, bbox_inches="tight")
        plt.close(fig_eb)
        print(f"Saved: {path_eb}")

    # ----- 1×2 Timoshenko: V_y and M_z, all six load cases overlaid -----
    all_load_cases_tims: list[tuple[int, str]] = [
        (6, "End"),
        (7, "Mid-span"),
        (8, "Quarter-span"),
        (9, "UDL"),
        (10, "Triangular"),
        (11, "Parabolic"),
    ]
    if any(job_id in job_to_result for job_id, _ in all_load_cases_tims):
        fig_tims, axes_tims = plt.subplots(1, 2, figsize=(12, 8), sharex=True)
        fig_tims.suptitle(
            r"$V_y$ and $M_z$: all load cases (Timoshenko jobs 6–11)",
            fontsize=14,
            fontweight="bold",
        )
        for col, ylabel in enumerate([r"$V_y$ [N]", r"$M_z$ [N·m]"]):
            ax = axes_tims[col]
            roark_lines_t: list = []
            for load_idx, (job_id, title) in enumerate(all_load_cases_tims):
                if job_id not in job_to_result:
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
                            continue
                    else:
                        continue
                else:
                    x_fem, forces_fem = x, forces

                V_fem = forces_fem[:, IDX_VY]
                M_fem = forces_fem[:, IDX_MZ]
                if job_id <= 8:
                    load_type = POINT_LOAD_CASES_TIMS[job_id - 6][2]
                    roark_at_nodes = roark_point_load_response(x_fem, L, E, I_z, P_point, load_type)
                    V_roark = ROARK_POINT_LOAD_SIGN * roark_at_nodes["shear"]
                    M_roark = ROARK_POINT_LOAD_SIGN * roark_at_nodes["moment"]
                    x_ana = _analytical_grid(x_fem.shape[0], L)
                    roark_fine = roark_point_load_response(x_ana, L, E, I_z, P_point, load_type)
                    V_ana = ROARK_POINT_LOAD_SIGN * roark_fine["shear"]
                    M_ana = ROARK_POINT_LOAD_SIGN * roark_fine["moment"]
                else:
                    roark_type = DIST_LOAD_CASES_TIMS[job_id - 9][2]
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

                formula_v = FORMULA_V_M_POINT.get((job_id, 0)) or FORMULA_V_M_DIST.get((job_id, 0), title)
                formula_m = FORMULA_V_M_POINT.get((job_id, 1)) or FORMULA_V_M_DIST.get((job_id, 1), title)
                color = COLORS_6[load_idx % len(COLORS_6)]
                vals_fem = V_fem if col == 0 else M_fem
                vals_ana = V_ana if col == 0 else M_ana
                formula = formula_v if col == 0 else formula_m
                ax.plot(x_fem, vals_fem, color=color, linestyle="--", linewidth=1.25, label=f"job_{job_id:04d} ({title})")
                (line,) = ax.plot(x_ana, vals_ana, "k-", linewidth=0.5, label=formula)
                roark_lines_t.append(line)
            ax.set_ylabel(ylabel)
            ax.grid(ls="--", alpha=0.6)
            if labelLines is not None and roark_lines_t:
                L_ax_t = float(np.max(roark_lines_t[0].get_xdata()))
                n_lines_t = len(roark_lines_t)
                xvals_t = np.linspace(0.2 * L_ax_t, 0.85 * L_ax_t, n_lines_t).tolist()
                labelLines(roark_lines_t, align=True, fontsize="x-small", xvals=xvals_t)
        h_left_t, l_left_t = _legend_num_anal(*axes_tims[0].get_legend_handles_labels())
        h_right_t, l_right_t = _legend_num_anal(*axes_tims[1].get_legend_handles_labels())
        fig_tims.tight_layout()
        fig_tims.subplots_adjust(top=0.88, bottom=0.32)
        axes_tims[0].legend(h_left_t, l_left_t, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        axes_tims[1].legend(h_right_t, l_right_t, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize="small")
        axes_tims[0].set_xlabel(r"$x$ [m]")
        axes_tims[1].set_xlabel(r"$x$ [m]")
        path_tims = out_dir / "roark_section_forces_convergence_all_loads_timoshenko.png"
        fig_tims.savefig(path_tims, dpi=300, bbox_inches="tight")
        plt.close(fig_tims)
        print(f"Saved: {path_tims}")

    if csv_rows:
        csv_path = out_dir / "roark_section_forces_verification_data.csv"
        header = "job_id,x,V_fem_N,V_roark_N,error_V_N,M_fem_Nm,M_roark_Nm,error_M_Nm"
        np.savetxt(csv_path, csv_rows, delimiter=",", header=header, comments="")
        print(f"Saved: {csv_path}")
    print("Roark section forces verification done.")


if __name__ == "__main__":
    run_section_forces_verification()
