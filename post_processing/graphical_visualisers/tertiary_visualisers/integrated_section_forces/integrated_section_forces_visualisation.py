"""
Integrated section forces visualisation (tertiary results, elemental resolution only).

Reads job_*/tertiary_results/elemental/integrated_section_forces.csv and plots
N, Vy, Vz, T, My, Mz vs position along structure (element midpoint x).
Elemental-level resolution only: B2 convention uses only filled square markers
at element midpoints (no nodal, no Gauss, no interpolant).
"""

from __future__ import annotations

import glob
import re
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

from pre_processing.parsing.grid_parser import GridParser  # type: ignore
from pre_processing.parsing.element_parser import ElementParser  # type: ignore
from post_processing.graphical_visualisers.resolution_plotting_utils import (
    get_element_node_coords,
    plot_elemental_points,
    ELEMENTAL_MARKER_SIZE,
)


def _get_node_coordinates(grid_obj: object) -> np.ndarray:
    if isinstance(grid_obj, dict) and "grid_dictionary" in grid_obj:
        inner = grid_obj["grid_dictionary"]
        if isinstance(inner, dict) and "coordinates" in inner:
            return inner["coordinates"]  # type: ignore[index]
    if isinstance(grid_obj, dict) and "node_coordinates" in grid_obj:
        return grid_obj["node_coordinates"]  # type: ignore[index]
    if hasattr(grid_obj, "node_coordinates"):
        return getattr(grid_obj, "node_coordinates")  # type: ignore[arg-type]
    raise KeyError("grid data does not contain 'grid_dictionary' -> 'coordinates'")


COMPONENTS = ["N", "Vy", "Vz", "T", "My", "Mz"]


class VisualiseIntegratedSectionForces:
    """Produce integrated section forces per element from tertiary_results."""

    _BLUE: Final[str] = "#4F81BD"

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "integrated_section_forces_plots"
        self.figure_output_dir.mkdir(exist_ok=True)
        self.results_dir: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
        self.jobs_dir: Final[Path] = PROJECT_ROOT / "jobs"

    def _plot(
        self,
        x_midpoints: np.ndarray,
        forces: np.ndarray,
        *,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        # forces shape (n_elements, 6) -> N, Vy, Vz, T, My, Mz
        fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
        axes = axes.ravel()
        for k, (ax, name) in enumerate(zip(axes, COMPONENTS)):
            plot_elemental_points(
                ax, x_midpoints, forces[:, k],
                color=self._BLUE, size=ELEMENTAL_MARKER_SIZE,
                label=name,
            )
            ax.set_ylabel(name)
            ax.grid(ls="--", alpha=0.6)
            ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
        axes[0].set_title("Axial N")
        axes[1].set_title("Shear Vy")
        axes[2].set_title("Shear Vz")
        axes[3].set_title("Torque T")
        axes[4].set_title("Moment My")
        axes[5].set_title("Moment Mz")
        for ax in axes[3:]:
            ax.set_xlabel(r"$x$ [m]")
        fig.suptitle(f"Integrated section forces {title_suffix}".strip(), fontweight="bold")
        fig.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=300)
            plt.close(fig)
        else:
            plt.show()

    def process_all(self) -> None:
        pattern = str(self.results_dir / "job_*" / "tertiary_results" / "elemental" / "integrated_section_forces.csv")
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No integrated_section_forces files found.")
            return

        for csv_path in csv_files:
            csv_file = Path(csv_path)
            job_dir = csv_file.parent.parent.parent  # elemental -> tertiary_results -> job_XXX
            m = re.match(r"job_(?P<id>\d+)_(?P<ts>[\d\-_]+_pid\d+_[a-f0-9]+)", job_dir.name)
            if not m:
                print(f"Skipping unrecognised folder '{job_dir.name}'")
                continue

            job_id, timestamp = m.group("id"), m.group("ts")
            grid_file = self.jobs_dir / f"job_{job_id}" / "grid.txt"
            element_file = self.jobs_dir / f"job_{job_id}" / "element.txt"

            if not grid_file.is_file() or not element_file.is_file():
                print(f"WARNING: Grid or element file missing for job {job_id}, skipping.")
                continue

            try:
                forces = np.genfromtxt(csv_file, delimiter=",", skip_header=1)
            except Exception as exc:
                print(f"Error reading {csv_file}: {exc}")
                continue

            if forces.ndim == 1:
                forces = forces.reshape(1, -1)
            if forces.shape[1] != 6:
                print(f"Expected 6 columns in {csv_file}, got {forces.shape[1]}, skipping.")
                continue

            grid = GridParser(str(grid_file), str(job_dir)).parse()
            try:
                _get_node_coordinates(grid)
            except (AttributeError, KeyError) as exc:
                print(f"WARNING: {exc} for job {job_id}, skipping.")
                continue

            elem_parsed = ElementParser(str(element_file), str(job_dir)).parse()
            element_dictionary = elem_parsed["element_dictionary"]
            grid_dictionary = grid["grid_dictionary"] if isinstance(grid, dict) else {}

            ids = element_dictionary.get("ids", np.arange(len(forces)))
            if len(ids) != len(forces):
                print(f"Element count mismatch for job {job_id}, skipping.")
                continue

            x_midpoints = []
            for i in range(len(forces)):
                elem_id = int(ids[i]) if hasattr(ids, "__getitem__") else i
                try:
                    node_coords = get_element_node_coords(
                        elem_id, element_dictionary, grid_dictionary
                    )
                    x_mid = float(np.mean(node_coords[:, 0]))
                    x_midpoints.append(x_mid)
                except Exception:
                    x_midpoints.append(np.nan)
            x_midpoints = np.array(x_midpoints)
            valid = np.isfinite(x_midpoints)
            if not np.any(valid):
                print(f"Could not compute midpoints for job {job_id}, skipping.")
                continue
            x_midpoints = x_midpoints[valid]
            forces = forces[valid]

            fig_name = f"integrated_section_forces_job_{job_id}_{timestamp}.png"
            self._plot(
                x_midpoints,
                forces,
                title_suffix=f"job_{job_id}_{timestamp}",
                save_path=self.figure_output_dir / fig_name,
            )
            print(f"Saved: {self.figure_output_dir / fig_name}")


if __name__ == "__main__":
    VisualiseIntegratedSectionForces().process_all()
