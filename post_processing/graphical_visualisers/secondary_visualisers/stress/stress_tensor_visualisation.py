"""
Stress visualisation utility.

Reads stress results from secondary_results and produces stress component profiles with:
- Nodal markers at node locations
- Continuous interpolated fields using element shape functions

Stress components: [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]
"""

from __future__ import annotations

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

# Walk upwards until we find the repo root (must contain pre_processing)
PROJECT_ROOT: Final[Path] = next(
    (p for p in SCRIPT_DIR.parents if (p / "pre_processing").is_dir()),
    SCRIPT_DIR.parents[4],
)
sys.path.append(str(PROJECT_ROOT))

# --- parsers ---------------------------------------------------------------#
from pre_processing.parsing.grid_parser import GridParser  # type: ignore
from pre_processing.parsing.element_parser import ElementParser  # type: ignore

# --- resolution plotting utilities -----------------------------------------#
from post_processing.graphical_visualisers.resolution_plotting_utils import (
    get_element_node_coords,
    plot_nodal_points,
    plot_interpolated_field,
    interpolate_field_nodal_to_continuous,
)


# ---------------------------------------------------------------------------#
#  Visualiser
# ---------------------------------------------------------------------------#
class VisualiseStress:
    """Produce stress component profiles from nodal_stress.csv files."""

    _BLUE: Final[str] = "#4F81BD"

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "stress_plots"
        self.figure_output_dir.mkdir(exist_ok=True)

        self.results_dir: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
        self.jobs_dir: Final[Path] = PROJECT_ROOT / "jobs"

    # ------------------------------------------------------------------#
    #  Internal helpers
    # ------------------------------------------------------------------#
    @staticmethod
    def _get_node_coordinates(grid_obj: object) -> np.ndarray:
        """
        Extract the (N, 3) array of node coordinates from the object returned
        by ``GridParser.parse()``.
        """
        # 1️⃣ official / nested layout
        if isinstance(grid_obj, dict) and "grid_dictionary" in grid_obj:
            inner = grid_obj["grid_dictionary"]
            if isinstance(inner, dict) and "coordinates" in inner:
                return inner["coordinates"]  # type: ignore[index]

        # 2️⃣ optional flat / attribute fall-backs
        if isinstance(grid_obj, dict) and "node_coordinates" in grid_obj:
            return grid_obj["node_coordinates"]  # type: ignore[index]
        if hasattr(grid_obj, "node_coordinates"):
            return getattr(grid_obj, "node_coordinates")  # type: ignore[arg-type]

        raise KeyError(
            "grid data does not contain 'grid_dictionary' -> 'coordinates'"
        )

    # ------------------------------------------------------------------#
    #  Plotting
    # ------------------------------------------------------------------#
    def _plot(
        self,
        stress: np.ndarray,
        node_positions: np.ndarray,
        element_dictionary: dict,
        grid_dictionary: dict,
        *,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        """
        Plot stress component profiles with nodal markers and interpolated continuous fields.
        
        Parameters
        ----------
        stress : np.ndarray
            Stress tensor components, shape (n_nodes, 6) with columns [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_xz]
        node_positions : np.ndarray
            Node x-coordinates, shape (n_nodes,)
        element_dictionary : dict
            Element dictionary with connectivity and types
        grid_dictionary : dict
            Grid dictionary with node coordinates
        title_suffix : str
            Suffix for plot title
        save_path : Optional[Path]
            Path to save figure
        """
        if stress.shape[1] != 6:
            raise ValueError("Stress must be shaped (n_nodes, 6)")

        x_min, x_max = float(node_positions.min()), float(node_positions.max())

        fig, axes = plt.subplots(3, 2, figsize=(15, 10), sharex=True)
        fig.suptitle(
            rf"Stress profiles${' – ' + title_suffix if title_suffix else ''}",
            fontsize=16,
            fontweight="bold",
        )

        # Stress component pairs: (normal_stress_index, normal_stress_label, shear_stress_index, shear_stress_label)
        pairs = [
            (0, r"$\sigma_{xx}$ [Pa]", 3, r"$\tau_{xy}$ [Pa]"),
            (1, r"$\sigma_{yy}$ [Pa]", 4, r"$\tau_{yz}$ [Pa]"),
            (2, r"$\sigma_{zz}$ [Pa]", 5, r"$\tau_{xz}$ [Pa]"),
        ]

        # Check if we have element data for interpolation
        has_elements = (
            element_dictionary and 
            "ids" in element_dictionary and 
            "connectivity" in element_dictionary and
            grid_dictionary and
            "coordinates" in grid_dictionary
        )

        # Process each component pair
        for i, (ax_l, ax_r, (norm_idx, norm_lbl, shear_idx, shear_lbl)) in enumerate(
            zip(axes[:, 0], axes[:, 1], pairs)
        ):
            # Get stress components
            norm_stress = stress[:, norm_idx]
            shear_stress = stress[:, shear_idx]

            # Plot interpolated continuous fields for each element (if available)
            if has_elements:
                element_ids = element_dictionary["ids"]
                for elem_id in element_ids:
                    try:
                        # Get element node coordinates
                        elem_node_coords = get_element_node_coords(
                            elem_id, element_dictionary, grid_dictionary
                        )
                        
                        # Get node IDs for this element
                        elem_idx = int(np.where(element_dictionary["ids"] == elem_id)[0][0])
                        node_ids = element_dictionary["connectivity"][elem_idx]
                        
                        # Get stresses at element nodes
                        elem_norm_stress = norm_stress[node_ids]
                        elem_shear_stress = shear_stress[node_ids]
                        
                        # Interpolate normal stress field
                        x_interp_norm, norm_interp = interpolate_field_nodal_to_continuous(
                            nodal_values=elem_norm_stress,
                            element_id=elem_id,
                            element_dictionary=element_dictionary,
                            grid_dictionary=grid_dictionary,
                            n_points=50
                        )
                        
                        # Interpolate shear stress field
                        x_interp_shear, shear_interp = interpolate_field_nodal_to_continuous(
                            nodal_values=elem_shear_stress,
                            element_id=elem_id,
                            element_dictionary=element_dictionary,
                            grid_dictionary=grid_dictionary,
                            n_points=50
                        )
                        
                        # Plot interpolated fields (label only once per subplot)
                        plot_interpolated_field(
                            ax_l, x_interp_norm, norm_interp,
                            linestyle='-', linewidth=2.0, alpha=0.7,
                            color=self._BLUE,
                            label="Interpolated (shape functions)" if i == 0 and elem_id == element_ids[0] else None
                        )
                        plot_interpolated_field(
                            ax_r, x_interp_shear, shear_interp,
                            linestyle='-', linewidth=2.0, alpha=0.7,
                            color=self._BLUE, label=None
                        )
                    except Exception:
                        # Skip elements that fail
                        continue

            # Plot nodal markers - label only once
            plot_nodal_points(
                ax_l, node_positions.reshape(-1, 1), norm_stress,
                marker='s', color=self._BLUE, size=70.0, alpha=0.9,
                edgecolors='black', linewidths=1.0,
                label="Nodes" if i == 0 else None
            )
            plot_nodal_points(
                ax_r, node_positions.reshape(-1, 1), shear_stress,
                marker='s', color=self._BLUE, size=70.0, alpha=0.9,
                edgecolors='black', linewidths=1.0,
                label=None  # Only label on left subplot
            )

            # Beam-end anchors + baseline
            for ax in (ax_l, ax_r):
                ax.plot([x_min], [0], marker="o", color="k", ms=3, zorder=3)
                ax.plot([x_max], [0], marker="o", color="k", ms=3, zorder=3)
                ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
                ax.grid(ls="--", alpha=0.6)

            ax_l.set_ylabel(norm_lbl)
            ax_r.set_ylabel(shear_lbl)

            if i == 0:
                ax_l.set_title("Normal stress components", fontweight="bold")
                ax_r.set_title("Shear stress components", fontweight="bold")

        axes[-1, 0].set_xlabel(r"$x$ [m]")
        axes[-1, 1].set_xlabel(r"$x$ [m]")
        
        # Add unified legend at bottom of figure
        if has_elements:
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], linestyle='-', linewidth=2.0, color=self._BLUE, 
                       label='Interpolated (shape functions)'),
                Line2D([0], [0], marker='s', linestyle='None', markersize=8, 
                       color=self._BLUE, markeredgecolor='black', markeredgewidth=1.0,
                       label='Nodes'),
                Line2D([0], [0], marker='o', linestyle='None', markersize=5, 
                       color='red', label='Gauss points'),
            ]
            fig.legend(handles=legend_elements, loc='lower center', ncol=3, 
                      fontsize=9, frameon=True, bbox_to_anchor=(0.5, 0.02))
        
        fig.tight_layout()
        fig.subplots_adjust(top=0.9, bottom=0.08 if has_elements else 0.1)

        if save_path:
            fig.savefig(save_path, dpi=300)
            plt.close(fig)
        else:
            plt.show()

    # ------------------------------------------------------------------#
    #  CSV helper
    # ------------------------------------------------------------------#
    @staticmethod
    def _read_nodal_stress(file: Path) -> Optional[np.ndarray]:
        """Read nodal stress CSV file."""
        try:
            stress = np.genfromtxt(file, delimiter=",", skip_header=1)
            if stress.ndim == 1:
                stress = stress.reshape(-1, 6)
            if stress.shape[1] != 6:
                raise ValueError("Stress must have 6 components")
            return stress
        except Exception as exc:
            print(f"Error reading {file}: {exc}")
            return None

    # ------------------------------------------------------------------#
    #  Driver
    # ------------------------------------------------------------------#
    def process_all(self) -> None:
        pattern = str(self.results_dir / "job_*" / "secondary_results" / "nodal" / "nodal_stress.csv")
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No stress files found.")
            return

        for csv_path in csv_files:
            csv_file = Path(csv_path)
            job_dir = csv_file.parent.parent.parent
            m = re.match(r"job_(?P<id>\d+)_(?P<ts>[\d\-_]+_pid\d+_[a-f0-9]+)", job_dir.name)
            if not m:
                print(f"Skipping unrecognised folder '{job_dir.name}'")
                continue

            job_id, timestamp = m.group("id"), m.group("ts")
            grid_file = self.jobs_dir / f"job_{job_id}" / "grid.txt"
            element_file = self.jobs_dir / f"job_{job_id}" / "element.txt"

            print(f"> Processing job {job_id} ({timestamp})")

            # ---- Stresses ----------------------------------------------- #
            stress = self._read_nodal_stress(csv_file)
            if stress is None:
                print(f"WARNING: Could not read stresses for job {job_id}, skipping.")
                continue

            # ---- Geometry (grid) ---------------------------------------- #
            if not grid_file.is_file():
                print(f"WARNING: Grid file missing for job {job_id}, skipping.")
                continue

            grid = GridParser(str(grid_file), str(job_dir)).parse()
            try:
                node_coords = self._get_node_coordinates(grid)
            except (AttributeError, KeyError) as exc:
                print(f"WARNING: {exc} for job {job_id}, skipping.")
                continue

            # ---- Element data (for interpolation) ---------------------- #
            if not element_file.is_file():
                print(f"WARNING: Element file missing for job {job_id}, skipping interpolation.")
                element_dictionary = None
                grid_dictionary = None
            else:
                try:
                    element_parsed = ElementParser(str(element_file), str(job_dir)).parse()
                    element_dictionary = element_parsed["element_dictionary"]
                    grid_dictionary = grid["grid_dictionary"]
                except Exception as exc:
                    print(f"WARNING: Could not parse element file for job {job_id}: {exc}")
                    print("   Plotting nodal points only (no interpolation).")
                    element_dictionary = None
                    grid_dictionary = None

            # ---- Plot ---------------------------------------------------- #
            fig_name = f"stress_job_{job_id}_{timestamp}.png"
            self._plot(
                stress,
                node_coords[:, 0],
                element_dictionary=element_dictionary if element_dictionary else {},
                grid_dictionary=grid_dictionary if grid_dictionary else {},
                title_suffix=f"job_{job_id}_{timestamp}",
                save_path=self.figure_output_dir / fig_name,
            )


if __name__ == "__main__":
    VisualiseStress().process_all()

