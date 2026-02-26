# post_processing/verification_visualisers/deformation_convergence.py

import glob
import re
import sys
from pathlib import Path
from typing import Final, Optional

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend so script saves and exits without blocking
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------#
#  Project paths
# ---------------------------------------------------------------------------#
SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent

PROJECT_ROOT: Final[Path] = next(
    (p for p in SCRIPT_DIR.parents if (p / "pre_processing").is_dir()),
    SCRIPT_DIR.parents[4],
)
sys.path.append(str(PROJECT_ROOT))
sys.path.insert(0, str(SCRIPT_DIR / "roark_utilities"))

from pre_processing.parsing.grid_parser import GridParser  # type: ignore
from roarks_formulas_euler_bernoulli_point import roark_point_load_response  # type: ignore
from roarks_formulas_euler_bernoulli_distributed import roark_distributed_load_response  # type: ignore
from roarks_formulas_timoshenko_point import timoshenko_point_load_response  # type: ignore
from roarks_formulas_timoshenko_distributed import timoshenko_distributed_load_response  # type: ignore

# Beam parameters (must match job material/section and gci_richardson_roark_report)
E: Final[float] = 2.1e11  # Pa
I_z: Final[float] = 2.08769e-06  # m^4
P_point: Final[float] = -500.0  # N
w_dist: Final[float] = 500.0  # N/m
# Timoshenko (match shear_deformable_verification / mesh section)
A: Final[float] = 0.00131  # m^2
G: Final[float] = 8.1e10  # Pa
K_S: Final[float] = 5.0 / 6.0
ROARK_LOAD_TYPE: Final[dict[int, str]] = {0: "end", 1: "mid", 2: "quarter"}
DIST_JOBS_ROARK_TYPE: Final[dict[int, str]] = {3: "udl", 4: "triangular", 5: "parabolic"}
# Timoshenko job 6–11: point 6,7,8 and distributed 9,10,11
ROARK_LOAD_TYPE_TIMS: Final[dict[int, str]] = {6: "end", 7: "mid", 8: "quarter"}
DIST_JOBS_ROARK_TYPE_TIMS: Final[dict[int, str]] = {9: "udl", 10: "triangular", 11: "parabolic"}


def _roark_uy_theta_at_x(x: np.ndarray, L: float, job_id: int) -> tuple[np.ndarray, np.ndarray]:
    """Roark deflection [m] and rotation [rad] at node x for point-load job_id (0,1,2)."""
    load_type = ROARK_LOAD_TYPE[job_id]
    roark = roark_point_load_response(x, L, E, I_z, P_point, load_type)
    return roark["deflection"], roark["rotation"]


def _roark_uy_theta_distributed_at_x(x: np.ndarray, L: float, job_id: int) -> tuple[np.ndarray, np.ndarray]:
    """Roark deflection [m] and rotation [rad] at node x for distributed-load job_id (3,4,5)."""
    roark_type = DIST_JOBS_ROARK_TYPE[job_id]
    order = np.argsort(x)
    x_sorted = x[order]
    roark = roark_distributed_load_response(x_sorted, L, E, I_z, w_dist, roark_type)
    uy = np.interp(x, x_sorted, roark["deflection"])
    th = np.interp(x, x_sorted, roark["rotation"])
    return uy, th


def _timoshenko_uy_theta_point_at_x(x: np.ndarray, L: float, job_id: int) -> tuple[np.ndarray, np.ndarray]:
    """Timoshenko Roark deflection [m] and rotation [rad] at x for point-load job_id (6,7,8)."""
    load_type = ROARK_LOAD_TYPE_TIMS[job_id]
    roark = timoshenko_point_load_response(x, L, E, I_z, P_point, A, G, K_S, load_type)
    return roark["deflection"], roark["rotation"]


def _timoshenko_uy_theta_distributed_at_x(x: np.ndarray, L: float, job_id: int) -> tuple[np.ndarray, np.ndarray]:
    """Timoshenko Roark deflection [m] and rotation [rad] at x for distributed-load job_id (9,10,11)."""
    roark_type = DIST_JOBS_ROARK_TYPE_TIMS[job_id]
    order = np.argsort(x)
    x_sorted = x[order]
    roark = timoshenko_distributed_load_response(x_sorted, L, E, I_z, w_dist, A, G, roark_type, K_S)
    uy = np.interp(x, x_sorted, roark["deflection"])
    th = np.interp(x, x_sorted, roark["rotation"])
    return uy, th


class VisualiseDeformationConvergence:
    """
    Overlay translation / rotation profiles from *_U_global.csv files
    on a single figure, assuming each belongs to a convergence study.
    """

    _LINESTYLES: Final[list[str]] = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
    _COLORS: Final[list[str]] = ["#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#4BACC6"]
    _COLORS_6: Final[list[str]] = ["#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#4BACC6", "#F79646"]

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "deformation_plots"
        self.figure_output_dir.mkdir(exist_ok=True)

        self.results_dir: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
        self.jobs_dir: Final[Path] = PROJECT_ROOT / "jobs"

    @staticmethod
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

    @staticmethod
    def _read_U_global(file: Path) -> Optional[np.ndarray]:
        try:
            vals = np.genfromtxt(file, delimiter=",", skip_header=1, usecols=1)
            if vals.size % 6:
                raise ValueError("DOF count not divisible by 6")
            return vals.reshape(-1, 6)
        except Exception as exc:
            print(f"Error reading {file}: {exc}")
            return None

    def process_convergence_plot(self, scale: float = 1.0) -> None:
        pattern = str(self.results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No deformation files found.")
            return

        # Map base_id -> (csv_path, job_dir, job_input_name); keep result with largest n per base
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

        fig, axes = plt.subplots(3, 2, figsize=(15, 10), sharex=True)
        fig.suptitle(
            "Deformation convergence: point load at end / mid-span / quarter-span",
            fontsize=16,
            fontweight="bold",
        )

        # (job_id, a as fraction of L, row title)
        load_cases = [
            (0, 1.0, "End load"),
            (1, 0.5, "Mid-span load"),
            (2, 0.25, "Quarter-span load"),
        ]

        for row, (job_id, _a_frac, title) in enumerate(load_cases):
            if job_id not in job_to_result:
                print(f"WARNING: No result for job_{job_id:04d}, skipping row '{title}'.")
                axes[row, 0].set_title(title, fontweight="bold")
                axes[row, 1].set_visible(True)
                continue

            csv_file, job_dir, job_input_name = job_to_result[job_id]
            grid_file = self.jobs_dir / job_input_name / "grid.txt"
            if not grid_file.is_file():
                print(f"WARNING: Missing grid file for job_{job_id:04d}, skipping row '{title}'.")
                continue

            U = self._read_U_global(csv_file)
            if U is None:
                print(f"WARNING: Failed to read displacements for job_{job_id:04d}, skipping.")
                continue

            grid = GridParser(str(grid_file), str(self.jobs_dir / job_input_name)).parse()
            try:
                node_coords = self._get_node_coordinates(grid)
            except Exception as exc:
                print(f"WARNING: {exc} for job_{job_id:04d}, skipping.")
                continue

            x = node_coords[:, 0]
            if x.shape[0] != U.shape[0]:
                print(f"ERROR: Mismatch x/U for job_{job_id:04d}, skipping.")
                continue

            L = float(np.max(x))

            uy_mm = U[:, 1] * 1000 * scale
            thetaz_deg = np.degrees(U[:, 5]) * scale

            uy_roark_m, thetaz_roark_rad = _roark_uy_theta_at_x(x, L, job_id)
            uy_roark_mm = uy_roark_m * 1000
            thetaz_roark_deg = np.degrees(thetaz_roark_rad)

            axes[row, 0].plot(x, uy_mm, color=self._COLORS[row], linestyle="-", label=f"FEM (job_{job_id:04d})")
            axes[row, 0].plot(x, uy_roark_mm, "k--", label="Roark")
            axes[row, 1].plot(x, thetaz_deg, color=self._COLORS[row], linestyle="-", label=f"FEM (job_{job_id:04d})")
            axes[row, 1].plot(x, thetaz_roark_deg, "k--", label="Roark")

            axes[row, 0].set_ylabel(r"$u_y$ [mm]")
            axes[row, 1].set_ylabel(r"$\theta_z$ [deg]")
            axes[row, 0].set_title(title, fontweight="bold")
            axes[row, 1].set_title(title, fontweight="bold")
            axes[row, 0].grid(ls="--", alpha=0.6)
            axes[row, 1].grid(ls="--", alpha=0.6)
            axes[row, 0].legend(loc="upper right", fontsize="small")
            axes[row, 1].legend(loc="upper right", fontsize="small")

        axes[-1, 0].set_xlabel(r"$x$ [m]")
        axes[-1, 1].set_xlabel(r"$x$ [m]")

        fig.tight_layout()
        fig.subplots_adjust(top=0.92)

        plot_path = self.figure_output_dir / "deformation_convergence_overlay.png"
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)
        print(f"Saved: {plot_path}")

        # Single 2×3 subplot: all 6 deformation DOFs, all point-load jobs overlaid
        fig_2x3, axes_2x3 = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
        fig_2x3.suptitle(
            "Deformation DOFs: all point-load jobs (End / Mid-span / Quarter)",
            fontsize=14,
            fontweight="bold",
        )
        dof_specs = [
            (0, r"$u_x$ [mm]", "disp"),
            (1, r"$u_y$ [mm]", "disp"),
            (2, r"$u_z$ [mm]", "disp"),
            (3, r"$\theta_x$ [deg]", "rot"),
            (4, r"$\theta_y$ [deg]", "rot"),
            (5, r"$\theta_z$ [deg]", "rot"),
        ]
        for col, (dof_idx, ylabel, kind) in enumerate(dof_specs):
            row_2x3, col_2x3 = col // 3, col % 3
            ax = axes_2x3[row_2x3, col_2x3]
            for load_idx, (job_id, _a_frac, title) in enumerate(load_cases):
                if job_id not in job_to_result:
                    continue
                csv_file, _job_dir, job_input_name = job_to_result[job_id]
                grid_file = self.jobs_dir / job_input_name / "grid.txt"
                if not grid_file.is_file():
                    continue
                U = self._read_U_global(csv_file)
                if U is None:
                    continue
                grid = GridParser(str(grid_file), str(self.jobs_dir / job_input_name)).parse()
                try:
                    node_coords = self._get_node_coordinates(grid)
                except Exception:
                    continue
                x = node_coords[:, 0]
                if x.shape[0] != U.shape[0]:
                    continue
                L = float(np.max(x))
                if kind == "rot":
                    vals = np.degrees(U[:, dof_idx]) * scale
                else:
                    vals = U[:, dof_idx] * 1000 * scale
                ax.plot(x, vals, color=self._COLORS[load_idx], linestyle="-", label=f"job_{job_id:04d} ({title})")
                if dof_idx == 1:
                    uy_roark_m, _ = _roark_uy_theta_at_x(x, L, job_id)
                    ax.plot(x, uy_roark_m * 1000, "k--", label="Roark" if load_idx == 0 else None)
                elif dof_idx == 5:
                    _, th_roark_rad = _roark_uy_theta_at_x(x, L, job_id)
                    ax.plot(x, np.degrees(th_roark_rad), "k--", label="Roark" if load_idx == 0 else None)
            ax.set_ylabel(ylabel)
            ax.grid(ls="--", alpha=0.6)
            ax.legend(loc="upper right", fontsize="small")
        for ax in axes_2x3[1, :]:
            ax.set_xlabel(r"$x$ [m]")
        fig_2x3.tight_layout()
        fig_2x3.subplots_adjust(top=0.90)
        path_2x3 = self.figure_output_dir / "deformation_convergence_all_dofs_2x3.png"
        fig_2x3.savefig(path_2x3, dpi=300)
        plt.close(fig_2x3)
        print(f"Saved: {path_2x3}")

        # 1×2: u_y and θ_z only, all six load cases (End, Mid, Quarter, UDL, Triangular, Parabolic) overlaid
        all_load_cases = [
            (0, "End"),
            (1, "Mid-span"),
            (2, "Quarter"),
            (3, "UDL"),
            (4, "Triangular"),
            (5, "Parabolic"),
        ]
        fig_1x2, axes_1x2 = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
        fig_1x2.suptitle(
            r"$u_y$ and $\theta_z$: all load cases (point + distributed)",
            fontsize=14,
            fontweight="bold",
        )
        for col, (dof_idx, ylabel) in enumerate([(1, r"$u_y$ [mm]"), (5, r"$\theta_z$ [deg]")]):
            ax = axes_1x2[col]
            for load_idx, (job_id, title) in enumerate(all_load_cases):
                if job_id not in job_to_result:
                    continue
                csv_file, _job_dir, job_input_name = job_to_result[job_id]
                grid_file = self.jobs_dir / job_input_name / "grid.txt"
                if not grid_file.is_file():
                    continue
                U = self._read_U_global(csv_file)
                if U is None:
                    continue
                grid = GridParser(str(grid_file), str(self.jobs_dir / job_input_name)).parse()
                try:
                    node_coords = self._get_node_coordinates(grid)
                except Exception:
                    continue
                x = node_coords[:, 0]
                if x.shape[0] != U.shape[0]:
                    continue
                L = float(np.max(x))
                color = self._COLORS_6[load_idx % len(self._COLORS_6)]
                if dof_idx == 1:
                    vals = U[:, 1] * 1000 * scale
                else:
                    vals = np.degrees(U[:, 5]) * scale
                ax.plot(x, vals, color=color, linestyle="-", label=f"job_{job_id:04d} ({title})")
                if job_id in ROARK_LOAD_TYPE:
                    uy_r, th_r = _roark_uy_theta_at_x(x, L, job_id)
                else:
                    uy_r, th_r = _roark_uy_theta_distributed_at_x(x, L, job_id)
                if dof_idx == 1:
                    ax.plot(x, uy_r * 1000, "k--", label="Roark" if load_idx == 0 else None)
                else:
                    ax.plot(x, np.degrees(th_r), "k--", label="Roark" if load_idx == 0 else None)
            ax.set_ylabel(ylabel)
            ax.grid(ls="--", alpha=0.6)
            ax.legend(loc="lower left", fontsize="small")
        axes_1x2[0].set_xlabel(r"$x$ [m]")
        axes_1x2[1].set_xlabel(r"$x$ [m]")
        fig_1x2.tight_layout()
        fig_1x2.subplots_adjust(top=0.88)
        path_1x2 = self.figure_output_dir / "deformation_convergence_uy_theta_all_loads.png"
        fig_1x2.savefig(path_1x2, dpi=300)
        plt.close(fig_1x2)
        print(f"Saved: {path_1x2}")

        # 1×2 Timoshenko: u_y and θ_z, all six Timoshenko load cases (6–11) overlaid vs Roark Timoshenko
        all_load_cases_tims = [
            (6, "End"),
            (7, "Mid-span"),
            (8, "Quarter"),
            (9, "UDL"),
            (10, "Triangular"),
            (11, "Parabolic"),
        ]
        if any(job_id in job_to_result for job_id, _ in all_load_cases_tims):
            fig_tims, axes_tims = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
            fig_tims.suptitle(
                r"$u_y$ and $\theta_z$: Timoshenko jobs 6–11 vs Roark Timoshenko",
                fontsize=14,
                fontweight="bold",
            )
            for col, (dof_idx, ylabel) in enumerate([(1, r"$u_y$ [mm]"), (5, r"$\theta_z$ [deg]")]):
                ax = axes_tims[col]
                for load_idx, (job_id, title) in enumerate(all_load_cases_tims):
                    if job_id not in job_to_result:
                        continue
                    csv_file, _job_dir, job_input_name = job_to_result[job_id]
                    grid_file = self.jobs_dir / job_input_name / "grid.txt"
                    if not grid_file.is_file():
                        continue
                    U = self._read_U_global(csv_file)
                    if U is None:
                        continue
                    grid = GridParser(str(grid_file), str(self.jobs_dir / job_input_name)).parse()
                    try:
                        node_coords = self._get_node_coordinates(grid)
                    except Exception:
                        continue
                    x = node_coords[:, 0]
                    if x.shape[0] != U.shape[0]:
                        continue
                    L = float(np.max(x))
                    color = self._COLORS_6[load_idx % len(self._COLORS_6)]
                    if dof_idx == 1:
                        vals = U[:, 1] * 1000 * scale
                    else:
                        vals = np.degrees(U[:, 5]) * scale
                    ax.plot(x, vals, color=color, linestyle="-", label=f"job_{job_id:04d} ({title})")
                    if job_id in ROARK_LOAD_TYPE_TIMS:
                        uy_r, th_r = _timoshenko_uy_theta_point_at_x(x, L, job_id)
                    else:
                        uy_r, th_r = _timoshenko_uy_theta_distributed_at_x(x, L, job_id)
                    if dof_idx == 1:
                        ax.plot(x, uy_r * 1000, "k--", label="Roark Timoshenko" if load_idx == 0 else None)
                    else:
                        ax.plot(x, np.degrees(th_r), "k--", label="Roark Timoshenko" if load_idx == 0 else None)
                ax.set_ylabel(ylabel)
                ax.grid(ls="--", alpha=0.6)
                ax.legend(loc="lower left", fontsize="small")
            axes_tims[0].set_xlabel(r"$x$ [m]")
            axes_tims[1].set_xlabel(r"$x$ [m]")
            fig_tims.tight_layout()
            fig_tims.subplots_adjust(top=0.88)
            path_tims = self.figure_output_dir / "deformation_convergence_uy_theta_all_loads_timoshenko.png"
            fig_tims.savefig(path_tims, dpi=300)
            plt.close(fig_tims)
            print(f"Saved: {path_tims}")

        # CSV export: FEM vs Roark and divergence analysis
        csv_rows: list[list[float]] = []
        for load_idx, (job_id, _a_frac, title) in enumerate(load_cases):
            if job_id not in job_to_result:
                continue
            csv_file, job_dir, job_input_name = job_to_result[job_id]
            grid_file = self.jobs_dir / job_input_name / "grid.txt"
            if not grid_file.is_file():
                continue
            U = self._read_U_global(csv_file)
            if U is None:
                continue
            grid = GridParser(str(grid_file), str(self.jobs_dir / job_input_name)).parse()
            try:
                node_coords = self._get_node_coordinates(grid)
            except Exception:
                continue
            x = node_coords[:, 0]
            if x.shape[0] != U.shape[0]:
                continue
            L = float(np.max(x))
            uy_fem_mm = U[:, 1] * 1000 * scale
            thetaz_fem_deg = np.degrees(U[:, 5]) * scale
            uy_roark_m, thetaz_roark_rad = _roark_uy_theta_at_x(x, L, job_id)
            uy_roark_mm = uy_roark_m * 1000
            thetaz_roark_deg = np.degrees(thetaz_roark_rad)
            err_uy = uy_fem_mm - uy_roark_mm
            err_theta = thetaz_fem_deg - thetaz_roark_deg
            for i in range(len(x)):
                csv_rows.append([
                    float(job_id), x[i], uy_fem_mm[i], uy_roark_mm[i], err_uy[i],
                    thetaz_fem_deg[i], thetaz_roark_deg[i], err_theta[i],
                ])
        if csv_rows:
            csv_path = self.figure_output_dir / "deformation_convergence_data.csv"
            header = "job_id,x,uy_fem_mm,uy_roark_mm,error_uy_mm,theta_z_fem_deg,theta_z_roark_deg,error_theta_z_deg"
            np.savetxt(csv_path, csv_rows, delimiter=",", header=header, comments="")
            print(f"Saved: {csv_path}")
            # Divergence summary per job
            arr = np.array(csv_rows)
            print("\n--- FEM vs Roark (point load) ---")
            for job_id in (0, 1, 2):
                mask = arr[:, 0] == job_id
                if not np.any(mask):
                    continue
                j = arr[mask]
                max_uy = np.max(np.abs(j[:, 4]))
                rms_uy = np.sqrt(np.mean(j[:, 4] ** 2))
                max_th = np.max(np.abs(j[:, 7]))
                rms_th = np.sqrt(np.mean(j[:, 7] ** 2))
                tip = j[-1]
                tip_uy_err = tip[4]
                tip_th_err = tip[7]
                names = ["End", "Mid-span", "Quarter"]
                print(f"  job_{job_id:04d} ({names[job_id]}):")
                print(f"    u_y:     max|error| = {max_uy:.6f} mm,  RMS = {rms_uy:.6f} mm,  tip error = {tip_uy_err:.6f} mm")
                print(f"    theta_z: max|error| = {max_th:.6f} deg, RMS = {rms_th:.6f} deg, tip error = {tip_th_err:.6f} deg")
            print("")

        # Three separate 3×2 "all DOF" figures (one per load case)
        pairs = [
            (0, r"$u_x$ [mm]", 3, r"$\theta_x$ [deg]"),
            (1, r"$u_y$ [mm]", 5, r"$\theta_z$ [deg]"),
            (2, r"$u_z$ [mm]", 4, r"$\theta_y$ [deg]"),
        ]
        suffixes = ["end", "midspan", "quarter"]
        for load_idx, (job_id, _a_frac, title) in enumerate(load_cases):
            if job_id not in job_to_result:
                continue
            csv_file, job_dir, job_input_name = job_to_result[job_id]
            grid_file = self.jobs_dir / job_input_name / "grid.txt"
            if not grid_file.is_file():
                continue
            U = self._read_U_global(csv_file)
            if U is None:
                continue
            grid = GridParser(str(grid_file), str(self.jobs_dir / job_input_name)).parse()
            try:
                node_coords = self._get_node_coordinates(grid)
            except Exception:
                continue
            x = node_coords[:, 0]
            if x.shape[0] != U.shape[0]:
                continue
            L = float(np.max(x))
            uy_roark_m, thetaz_roark_rad = _roark_uy_theta_at_x(x, L, job_id)
            uy_roark_mm = uy_roark_m * 1000
            thetaz_roark_deg = np.degrees(thetaz_roark_rad)

            fig_all, axes_all = plt.subplots(3, 2, figsize=(15, 10), sharex=True)
            fig_all.suptitle(f"All DOFs: {title} (job_{job_id:04d})", fontsize=16, fontweight="bold")
            color = self._COLORS[load_idx]
            for i, (disp_idx, disp_lbl, rot_idx, rot_lbl) in enumerate(pairs):
                axes_all[i, 0].plot(
                    x, U[:, disp_idx] * 1000 * scale, color=color, linestyle="-", label=f"FEM (job_{job_id:04d})"
                )
                axes_all[i, 1].plot(
                    x, np.degrees(U[:, rot_idx]) * scale, color=color, linestyle="-", label=f"FEM (job_{job_id:04d})"
                )
                if i == 1:
                    axes_all[i, 0].plot(x, uy_roark_mm, "k--", label="Roark")
                    axes_all[i, 1].plot(x, thetaz_roark_deg, "k--", label="Roark")
                axes_all[i, 0].set_ylabel(disp_lbl)
                axes_all[i, 1].set_ylabel(rot_lbl)
                axes_all[i, 0].set_title("Translation" if i == 0 else "", fontweight="bold")
                axes_all[i, 1].set_title("Rotation" if i == 0 else "", fontweight="bold")
                axes_all[i, 0].grid(ls="--", alpha=0.6)
                axes_all[i, 1].grid(ls="--", alpha=0.6)
                axes_all[i, 0].legend(loc="upper right", fontsize="small")
                axes_all[i, 1].legend(loc="upper right", fontsize="small")
            axes_all[-1, 0].set_xlabel(r"$x$ [m]")
            axes_all[-1, 1].set_xlabel(r"$x$ [m]")
            fig_all.tight_layout()
            fig_all.subplots_adjust(top=0.92)
            all_dof_path = self.figure_output_dir / f"deformation_convergence_all_dofs_{suffixes[load_idx]}.png"
            fig_all.savefig(all_dof_path, dpi=300)
            plt.close(fig_all)
            print(f"Saved: {all_dof_path}")


if __name__ == "__main__":
    VisualiseDeformationConvergence().process_convergence_plot()
