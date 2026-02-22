"""
Deformation visualisation utility – July 2025
---------------------------------------------

Reads every *_U_global.csv under post_processing/results/**/primary_results
and produces translation / rotation profiles with:
- Nodal markers at node locations
- Continuous interpolated fields using element shape functions

Mesh geometry is supplied by ``grid.txt`` and ``element.txt`` files parsed
with ``GridParser`` and ``ElementParser``.
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
    INTERPOLANT_LINEWIDTH,
    NODAL_MARKER_SIZE,
    LEGEND_MARKER_SIZE,
    LEGEND_MARKER_SIZE_SECONDARY,
)


# ---------------------------------------------------------------------------#
#  Visualiser
# ---------------------------------------------------------------------------#
class VisualiseDeformation:
    """Produce translation / rotation profiles from *_U_global.csv files."""

    _BLUE: Final[str] = "#4F81BD"

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "deformation_plots"
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

        Expected shape (dict only):

            grid_obj["grid_dictionary"]["coordinates"]

        Falls back to ``.node_coordinates`` or ``["node_coordinates"]`` if they
        ever appear in a future refactor.
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
        U: np.ndarray,
        node_positions: np.ndarray,
        element_dictionary: dict,
        grid_dictionary: dict,
        *,
        scale: float = 1.0,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        """
        Plot displacement/rotation profiles with nodal markers and interpolated continuous fields.
        
        Parameters
        ----------
        U : np.ndarray
            Global displacements, shape (n_nodes, 6) with columns [u_x, u_y, u_z, θ_x, θ_y, θ_z]
        node_positions : np.ndarray
            Node x-coordinates, shape (n_nodes,)
        element_dictionary : dict
            Element dictionary with connectivity and types
        grid_dictionary : dict
            Grid dictionary with node coordinates
        scale : float
            Scale factor for displacements
        title_suffix : str
            Suffix for plot title
        save_path : Optional[Path]
            Path to save figure
        """
        if U.shape[1] != 6:
            raise ValueError("U must be shaped (n_nodes, 6)")

        x_min, x_max = float(node_positions.min()), float(node_positions.max())

        fig, axes = plt.subplots(3, 2, figsize=(15, 10), sharex=True)
        fig.suptitle(
            rf"Deformation profiles${' – ' + title_suffix if title_suffix else ''}",
            fontsize=16,
            fontweight="bold",
        )

        # Component pairs: (displacement_index, displacement_label, rotation_index, rotation_label)
        pairs = [
            (0, r"$u_x(x)\ \mathrm{[mm]}$", 3, r"$\theta_x(x)\ [^\circ]$"),
            (1, r"$u_y(x)\ \mathrm{[mm]}$", 4, r"$\theta_y(x)\ [^\circ]$"),
            (2, r"$u_z(x)\ \mathrm{[mm]}$", 5, r"$\theta_z(x)\ [^\circ]$"),
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
        for i, (ax_l, ax_r, (disp_idx, disp_lbl, rot_idx, rot_lbl)) in enumerate(
            zip(axes[:, 0], axes[:, 1], pairs)
        ):
            # Get displacement and rotation components
            disp_values = U[:, disp_idx] * 1000 * scale  # Convert to mm
            rot_values = np.degrees(U[:, rot_idx]) * scale  # Convert to degrees

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
                        
                        # Get displacements at element nodes
                        # node_ids are indices into the node array
                        elem_disp = disp_values[node_ids]
                        elem_rot = rot_values[node_ids]
                        
                        # Interpolate displacement field
                        x_interp_disp, disp_interp = interpolate_field_nodal_to_continuous(
                            nodal_values=elem_disp,
                            element_id=elem_id,
                            element_dictionary=element_dictionary,
                            grid_dictionary=grid_dictionary,
                            n_points=50  # Fine interpolation for smooth curves
                        )
                        
                        # Interpolate rotation field
                        x_interp_rot, rot_interp = interpolate_field_nodal_to_continuous(
                            nodal_values=elem_rot,
                            element_id=elem_id,
                            element_dictionary=element_dictionary,
                            grid_dictionary=grid_dictionary,
                            n_points=50
                        )
                        
                        # Plot interpolated fields (B2: thin black line)
                        plot_interpolated_field(
                            ax_l, x_interp_disp, disp_interp,
                            linestyle='-', alpha=0.7,
                            label="Interpolated (shape functions)" if i == 0 and elem_id == element_ids[0] else None
                        )
                        plot_interpolated_field(
                            ax_r, x_interp_rot, rot_interp,
                            linestyle='-', alpha=0.7, label=None
                        )
                    except Exception as exc:
                        # Skip elements that fail (e.g. wrong/missing shape functions for element type)
                        import warnings
                        warnings.warn(
                            f"Deformation interpolation skipped for element {elem_id}: {exc}",
                            UserWarning,
                            stacklevel=2,
                        )
                        continue

            # Plot nodal markers (B2: explicit = solid triangle)
            plot_nodal_points(
                ax_l, node_positions.reshape(-1, 1), disp_values,
                color=self._BLUE, size=NODAL_MARKER_SIZE, alpha=0.9,
                nodal_data_type='explicit',
                label="Nodes" if i == 0 else None
            )
            plot_nodal_points(
                ax_r, node_positions.reshape(-1, 1), rot_values,
                color=self._BLUE, size=NODAL_MARKER_SIZE, alpha=0.9,
                nodal_data_type='explicit',
                label=None  # Only label on left subplot
            )

            # Beam-end anchors + baseline
            for ax in (ax_l, ax_r):
                ax.plot([x_min], [0], marker="o", color="k", ms=3, zorder=3)
                ax.plot([x_max], [0], marker="o", color="k", ms=3, zorder=3)
                ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
                ax.grid(ls="--", alpha=0.6)

            ax_l.set_ylabel(disp_lbl)
            ax_r.set_ylabel(rot_lbl)

            if i == 0:
                ax_l.set_title("Translation profiles", fontweight="bold")
                ax_r.set_title("Rotation profiles", fontweight="bold")

        axes[-1, 0].set_xlabel(r"$x$ [m]")
        axes[-1, 1].set_xlabel(r"$x$ [m]")
        
        # Add unified legend at bottom of figure (B2 convention; only resolution levels present)
        from matplotlib.lines import Line2D
        legend_elements = []
        if has_elements:
            legend_elements.append(
                Line2D([0], [0], linestyle='-', linewidth=INTERPOLANT_LINEWIDTH, color='black',
                       label='Interpolated (shape functions)'))
        legend_elements.append(
            Line2D([0], [0], marker='^', linestyle='None', markersize=LEGEND_MARKER_SIZE,
                   color=self._BLUE, markeredgecolor='black', markeredgewidth=1.0,
                   label='Nodes'))
        if legend_elements:
            fig.legend(handles=legend_elements, loc='lower center', ncol=len(legend_elements),
                      fontsize=9, frameon=True, bbox_to_anchor=(0.5, 0.06))

        fig.tight_layout()
        fig.subplots_adjust(top=0.9, bottom=0.14 if legend_elements else 0.1)

        if save_path:
            fig.savefig(save_path, dpi=300)
            plt.close(fig)
        else:
            plt.show()

    # ------------------------------------------------------------------#
    #  CSV helper
    # ------------------------------------------------------------------#
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

    # ------------------------------------------------------------------#
    #  Driver
    # ------------------------------------------------------------------#
    def process_all(self) -> None:
        pattern = str(self.results_dir / "job_*" / "primary_results" / "global" / "U_global.csv")
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No deformation files found.")
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

            # ---- Displacements ------------------------------------------ #
            U = self._read_U_global(csv_file)
            if U is None:
                print(f"WARNING: Could not read displacements for job {job_id}, skipping.")
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
            fig_name = f"deformation_job_{job_id}_{timestamp}.png"
            self._plot(
                U,
                node_coords[:, 0],
                element_dictionary=element_dictionary if element_dictionary else {},
                grid_dictionary=grid_dictionary if grid_dictionary else {},
                title_suffix=f"job_{job_id}_{timestamp}",
                save_path=self.figure_output_dir / fig_name,
            )


if __name__ == "__main__":
    VisualiseDeformation().process_all()