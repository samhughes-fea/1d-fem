# post_processing/validation_visualisers/deflection_tables/deformation_convergence.py

import glob
import re
import sys
from pathlib import Path
from typing import Final, Optional

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

from pre_processing.parsing.grid_parser import GridParser  # type: ignore


def _uy_theta_point_load(
    x: np.ndarray, a: float, L: float, F: float, E: float, I_z: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Euler–Bernoulli cantilever: point load F at x = a.
    Returns u_y [m] and θ_z [rad] on the given x grid.

    For x > a: M(x)=0 (no load to the right), so κ=dθ/dx=0 ⇒ θ=θ(a)=const
    and u_y = u_y(a) + θ(a)*(x - a).
    """
    x = np.asarray(x, dtype=float)
    uy = np.empty_like(x)
    theta_z = np.empty_like(x)
    mask_left = x <= a
    mask_right = ~mask_left
    # x ≤ a
    uy[mask_left] = (F / (6 * E * I_z)) * (x[mask_left] ** 2) * (3 * a - x[mask_left])
    theta_z[mask_left] = (F / (2 * E * I_z)) * x[mask_left] * (2 * a - x[mask_left])
    # x > a: M=0 ⇒ θ constant = θ(a), u_y = u_y(a) + θ(a)*(x - a)
    theta_a = (F / (2 * E * I_z)) * (a ** 2)
    uy_a = (F / (6 * E * I_z)) * (a ** 2) * (3 * a - a)  # = (F/(6*E*I_z))*a^2*2*a
    uy[mask_right] = uy_a + theta_a * (x[mask_right] - a)
    theta_z[mask_right] = theta_a
    return uy, theta_z


class VisualiseDeformationConvergence:
    """
    Overlay translation / rotation profiles from *_U_global.csv files
    on a single figure, assuming each belongs to a convergence study.
    """

    _LINESTYLES: Final[list[str]] = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
    _COLORS: Final[list[str]] = ["#4F81BD", "#C0504D", "#9BBB59", "#8064A2", "#4BACC6"]
    # Analytical curves use a much finer grid than FEM for smooth reference lines
    _ANALYTICAL_GRID_FACTOR: Final[int] = 25
    _ANALYTICAL_GRID_MIN_POINTS: Final[int] = 500

    @staticmethod
    def _analytical_grid(n_fem: int, L: float) -> np.ndarray:
        """Evaluation grid for analytical formulas: much finer than FEM for plotting."""
        n = max(n_fem * VisualiseDeformationConvergence._ANALYTICAL_GRID_FACTOR,
                VisualiseDeformationConvergence._ANALYTICAL_GRID_MIN_POINTS)
        return np.linspace(0.0, L, int(n))

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

        # Map job_id -> (csv_path, job_dir); keep first result per job
        job_to_result: dict[int, tuple[Path, Path]] = {}
        for csv_path in csv_files:
            csv_file = Path(csv_path)
            job_dir = csv_file.parent.parent.parent
            m = re.match(r"job_(?P<id>\d+)_(?P<ts>[\d\-_]+_pid\d+_[a-f0-9]+)", job_dir.name)
            if m:
                job_id = int(m.group("id"))
                if job_id not in job_to_result:
                    job_to_result[job_id] = (csv_file, job_dir)

        F = -500.0  # Load [N]
        E = 2e11  # Young's modulus [Pa]
        I_z = 2.08769e-06  # Moment of inertia [m^4]

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

        for row, (job_id, a_frac, title) in enumerate(load_cases):
            if job_id not in job_to_result:
                print(f"WARNING: No result for job_{job_id:04d}, skipping row '{title}'.")
                axes[row, 0].set_title(title, fontweight="bold")
                axes[row, 1].set_visible(True)
                continue

            csv_file, job_dir = job_to_result[job_id]
            grid_file = self.jobs_dir / f"job_{job_id:04d}" / "grid.txt"
            if not grid_file.is_file():
                print(f"WARNING: Missing grid file for job_{job_id:04d}, skipping row '{title}'.")
                continue

            U = self._read_U_global(csv_file)
            if U is None:
                print(f"WARNING: Failed to read displacements for job_{job_id:04d}, skipping.")
                continue

            grid = GridParser(str(grid_file), str(job_dir)).parse()
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
            a = a_frac * L

            uy_mm = U[:, 1] * 1000 * scale
            thetaz_deg = np.degrees(U[:, 5]) * scale

            x_ana = self._analytical_grid(x.shape[0], L)
            uy_analytical, thetaz_analytical_rad = _uy_theta_point_load(x_ana, a, L, F, E, I_z)
            thetaz_analytical_deg = np.degrees(thetaz_analytical_rad)

            axes[row, 0].plot(x, uy_mm, color=self._COLORS[row], linestyle="-", label=f"FEM (job_{job_id:04d})")
            axes[row, 0].plot(x_ana, uy_analytical * 1000, "k--", label="Analytical")
            axes[row, 1].plot(x, thetaz_deg, color=self._COLORS[row], linestyle="-", label=f"FEM (job_{job_id:04d})")
            axes[row, 1].plot(x_ana, thetaz_analytical_deg, "k--", label="Analytical")

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

        # CSV export: FEM vs analytical (correct formula per load case) and divergence analysis
        csv_rows: list[list[float]] = []
        for load_idx, (job_id, a_frac, title) in enumerate(load_cases):
            if job_id not in job_to_result:
                continue
            csv_file, job_dir = job_to_result[job_id]
            grid_file = self.jobs_dir / f"job_{job_id:04d}" / "grid.txt"
            if not grid_file.is_file():
                continue
            U = self._read_U_global(csv_file)
            if U is None:
                continue
            grid = GridParser(str(grid_file), str(job_dir)).parse()
            try:
                node_coords = self._get_node_coordinates(grid)
            except Exception:
                continue
            x = node_coords[:, 0]
            if x.shape[0] != U.shape[0]:
                continue
            L = float(np.max(x))
            a = a_frac * L
            uy_fem_mm = U[:, 1] * 1000 * scale
            thetaz_fem_deg = np.degrees(U[:, 5]) * scale
            uy_analytical, thetaz_analytical_rad = _uy_theta_point_load(x, a, L, F, E, I_z)
            uy_analytical_mm = uy_analytical * 1000
            thetaz_analytical_deg = np.degrees(thetaz_analytical_rad)
            err_uy = uy_fem_mm - uy_analytical_mm
            err_theta = thetaz_fem_deg - thetaz_analytical_deg
            for i in range(len(x)):
                csv_rows.append([
                    float(job_id), x[i], uy_fem_mm[i], uy_analytical_mm[i], err_uy[i],
                    thetaz_fem_deg[i], thetaz_analytical_deg[i], err_theta[i],
                ])
        if csv_rows:
            csv_path = self.figure_output_dir / "deformation_convergence_data.csv"
            header = "job_id,x,uy_fem_mm,uy_analytical_mm,error_uy_mm,theta_z_fem_deg,theta_z_analytical_deg,error_theta_z_deg"
            np.savetxt(csv_path, csv_rows, delimiter=",", header=header, comments="")
            print(f"Saved: {csv_path}")
            # Divergence summary per job
            arr = np.array(csv_rows)
            print("\n--- Divergence from theory (Euler–Bernoulli point load) ---")
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
        for load_idx, (job_id, a_frac, title) in enumerate(load_cases):
            if job_id not in job_to_result:
                continue
            csv_file, job_dir = job_to_result[job_id]
            grid_file = self.jobs_dir / f"job_{job_id:04d}" / "grid.txt"
            if not grid_file.is_file():
                continue
            U = self._read_U_global(csv_file)
            if U is None:
                continue
            grid = GridParser(str(grid_file), str(job_dir)).parse()
            try:
                node_coords = self._get_node_coordinates(grid)
            except Exception:
                continue
            x = node_coords[:, 0]
            if x.shape[0] != U.shape[0]:
                continue
            L = float(np.max(x))
            a = a_frac * L
            x_ana = self._analytical_grid(x.shape[0], L)
            uy_analytical, thetaz_analytical_rad = _uy_theta_point_load(x_ana, a, L, F, E, I_z)
            thetaz_analytical_deg = np.degrees(thetaz_analytical_rad)

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
                    axes_all[i, 0].plot(x_ana, uy_analytical * 1000, "k--", label="Analytical")
                    axes_all[i, 1].plot(x_ana, thetaz_analytical_deg, "k--", label="Analytical")
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