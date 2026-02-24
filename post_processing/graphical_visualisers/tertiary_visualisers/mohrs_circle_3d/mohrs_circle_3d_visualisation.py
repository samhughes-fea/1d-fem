"""
Mohr's circle 3D visualisation (tertiary results, Gaussian resolution).

Reads job_*/tertiary_results/gaussian/principal_stress/principal_stress_elem_*.csv
(σ1, σ2, σ3 at Gauss points) and plots the three Mohr's circles for a representative
stress state per job. The plotted state is the Gauss point with maximum shear
(σ1 − σ3)/2. Axes: normal stress σ (horizontal), shear stress τ (vertical).
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path
from typing import Final, Optional

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = next(
    (p for p in SCRIPT_DIR.parents if (p / "pre_processing").is_dir()),
    SCRIPT_DIR.parents[4],
)
sys.path.append(str(PROJECT_ROOT))

from post_processing.graphical_visualisers.job_discovery_utils import parse_job_result_dir_name


def _mohr_circle_xy(
    sigma_i: float,
    sigma_j: float,
    n_theta: int = 128,
) -> tuple[np.ndarray, np.ndarray]:
    """
    (σ, τ) points for one Mohr's circle between principal stresses σi and σj.
    Convention: σi ≥ σj; center ((σi+σj)/2, 0), radius (σi-σj)/2.
    """
    center = (sigma_i + sigma_j) / 2.0
    radius = abs(sigma_i - sigma_j) / 2.0
    theta = np.linspace(0, 2 * np.pi, n_theta)
    sigma = center + radius * np.cos(theta)
    tau = radius * np.sin(theta)
    return sigma, tau


def _draw_mohrs_circle_3d(
    ax: plt.Axes,
    s1: float,
    s2: float,
    s3: float,
    *,
    outer_color: str = "#4F81BD",
    mid_color: str = "#9BBB59",
    inner_color: str = "#C0504D",
) -> None:
    """
    Draw the three Mohr's circles for principal stresses s1 ≥ s2 ≥ s3.
    Outer (s1-s3), middle (s1-s2), inner (s2-s3).
    """
    # Ensure descending order (CSV is already sigma1, sigma2, sigma3)
    sigmas = np.sort([s1, s2, s3])[::-1]
    s1, s2, s3 = float(sigmas[0]), float(sigmas[1]), float(sigmas[2])

    # Outer circle: s1 - s3
    sigma_outer, tau_outer = _mohr_circle_xy(s1, s3)
    ax.plot(sigma_outer, tau_outer, color=outer_color, linewidth=2, label=r"$\sigma_1$–$\sigma_3$")
    # Middle: s1 - s2
    sigma_mid, tau_mid = _mohr_circle_xy(s1, s2)
    ax.plot(sigma_mid, tau_mid, color=mid_color, linewidth=1.5, label=r"$\sigma_1$–$\sigma_2$")
    # Inner: s2 - s3
    sigma_inner, tau_inner = _mohr_circle_xy(s2, s3)
    ax.plot(sigma_inner, tau_inner, color=inner_color, linewidth=1.5, label=r"$\sigma_2$–$\sigma_3$")

    # Principal stress points on σ-axis
    ax.scatter([s1, s2, s3], [0, 0, 0], color="k", s=40, zorder=5)
    ax.axhline(0, color="k", linewidth=0.8, linestyle="-")
    ax.axvline(0, color="k", linewidth=0.8, linestyle="--", alpha=0.7)


class VisualiseMohrsCircle3D:
    """Produce Mohr's circle 3D plots from tertiary principal_stress results."""

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "mohrs_circle_3d_plots"
        self.figure_output_dir.mkdir(exist_ok=True)
        self.results_dir: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
        self.jobs_dir: Final[Path] = PROJECT_ROOT / "jobs"

    def _plot(
        self,
        s1: float,
        s2: float,
        s3: float,
        *,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        fig, ax = plt.subplots(figsize=(7, 6))
        _draw_mohrs_circle_3d(ax, s1, s2, s3)
        ax.set_xlabel(r"Normal stress $\sigma$ [Pa]")
        ax.set_ylabel(r"Shear stress $\tau$ [Pa]")
        ax.set_title(f"Mohr's circle (3D) {title_suffix}".strip())
        ax.legend(loc="upper right")
        ax.grid(ls="--", alpha=0.6)
        ax.set_aspect("equal")
        fig.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=300)
            plt.close(fig)
        else:
            plt.show()

    def process_all(self) -> None:
        pattern = str(
            self.results_dir / "job_*" / "tertiary_results" / "gaussian" / "principal_stress" / "principal_stress_elem_*.csv"
        )
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No principal_stress files found.")
            return

        by_job: dict[str, list[Path]] = {}
        for p in csv_files:
            path = Path(p)
            job_dir = path.parent.parent.parent.parent
            by_job.setdefault(job_dir.name, []).append(path)

        for job_dirname, files in by_job.items():
            parsed = parse_job_result_dir_name(job_dirname)
            if not parsed:
                print(f"Skipping unrecognised folder '{job_dirname}'")
                continue

            job_id, job_folder_name, timestamp = parsed
            job_dir = self.results_dir / job_dirname
            grid_file = self.jobs_dir / job_folder_name / "grid.txt"
            element_file = self.jobs_dir / job_folder_name / "element.txt"

            if not grid_file.is_file() or not element_file.is_file():
                print(f"WARNING: Grid or element file missing for job {job_id}, skipping.")
                continue

            def elem_index(path: Path) -> int:
                stem = path.stem
                num = stem.split("_")[-1]
                return int(num)

            files_sorted = sorted(files, key=elem_index)

            # Load all principal stresses and pick the state with max (σ1 - σ3)/2
            best_s1, best_s2, best_s3 = None, None, None
            max_shear = -np.inf

            for csv_path in files_sorted:
                try:
                    data = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
                except Exception:
                    continue
                if data.ndim == 1:
                    data = data.reshape(1, -1)
                if data.shape[1] != 3:
                    continue
                for row in data:
                    s1, s2, s3 = row[0], row[1], row[2]
                    shear = abs(s1 - s3) / 2.0
                    if shear > max_shear:
                        max_shear = shear
                        best_s1, best_s2, best_s3 = s1, s2, s3

            if best_s1 is None:
                print(f"No valid principal stress data for job {job_id}, skipping.")
                continue

            fig_name = f"mohrs_circle_3d_{job_folder_name}_{timestamp}.png"
            save_path = self.figure_output_dir / fig_name
            self._plot(
                best_s1,
                best_s2,
                best_s3,
                title_suffix=f"{job_folder_name}_{timestamp} (max τ)",
                save_path=save_path,
            )
            print(f"Saved: {save_path}")


if __name__ == "__main__":
    VisualiseMohrsCircle3D().process_all()
