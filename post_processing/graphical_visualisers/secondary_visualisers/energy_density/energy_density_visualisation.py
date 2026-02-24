"""
Strain energy density visualisation utility.

Reads energy density results from secondary_results and produces energy density profiles with:
- Nodal markers at node locations (B2: projected = hollow circle)
- Gauss point markers when gaussian energy density is saved (B2: small solid circle)
- Continuous interpolated fields using element shape functions
"""

from __future__ import annotations

import glob
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
from post_processing.graphical_visualisers.job_discovery_utils import parse_job_result_dir_name
from post_processing.graphical_visualisers.resolution_plotting_utils import (
    get_element_node_coords,
    plot_nodal_points,
    plot_gauss_points,
    plot_interpolated_field,
    interpolate_field_nodal_to_continuous,
    natural_to_physical_coords,
    INTERPOLANT_LINEWIDTH,
    NODAL_MARKER_SIZE,
    LEGEND_MARKER_SIZE,
    LEGEND_MARKER_SIZE_SECONDARY,
    GAUSS_MARKER_SIZE,
)


# ---------------------------------------------------------------------------#
#  Visualiser
# ---------------------------------------------------------------------------#
class VisualiseEnergyDensity:
    """Produce strain energy density profiles from nodal_strain_energy_density.csv files."""

    _BLUE: Final[str] = "#4F81BD"

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "energy_density_plots"
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
        energy_density: np.ndarray,
        node_positions: np.ndarray,
        element_dictionary: dict,
        grid_dictionary: dict,
        *,
        x_gauss: Optional[np.ndarray] = None,
        energy_gauss: Optional[np.ndarray] = None,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        """
        Plot strain energy density profile with nodal markers, optional Gauss point markers, and interpolated continuous field.
        
        Parameters
        ----------
        energy_density : np.ndarray
            Strain energy density, shape (n_nodes,)
        node_positions : np.ndarray
            Node x-coordinates, shape (n_nodes,)
        element_dictionary : dict
            Element dictionary with connectivity and types
        grid_dictionary : dict
            Grid dictionary with node coordinates
        x_gauss : Optional[np.ndarray]
            Gauss point x-coordinates (when gaussian data is loaded)
        energy_gauss : Optional[np.ndarray]
            Strain energy density at Gauss points, shape (n_gauss,)
        title_suffix : str
            Suffix for plot title
        save_path : Optional[Path]
            Path to save figure
        """
        if energy_density.ndim != 1:
            raise ValueError("Energy density must be 1D array, shape (n_nodes,)")

        x_min, x_max = float(node_positions.min()), float(node_positions.max())
        has_gauss = (
            x_gauss is not None and energy_gauss is not None
            and len(x_gauss) == len(energy_gauss) and len(x_gauss) > 0
        )

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        fig.suptitle(
            rf"Strain energy density${' – ' + title_suffix if title_suffix else ''}",
            fontsize=16,
            fontweight="bold",
        )

        # Check if we have element data for interpolation
        has_elements = (
            element_dictionary and 
            "ids" in element_dictionary and 
            "connectivity" in element_dictionary and
            grid_dictionary and
            "coordinates" in grid_dictionary
        )

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
                    
                    # Get energy density at element nodes
                    elem_energy = energy_density[node_ids]
                    
                    # Interpolate energy density field
                    x_interp, energy_interp = interpolate_field_nodal_to_continuous(
                        nodal_values=elem_energy,
                        element_id=elem_id,
                        element_dictionary=element_dictionary,
                        grid_dictionary=grid_dictionary,
                        n_points=50
                    )
                    
                    # Plot interpolated field (B2: thin black line)
                    plot_interpolated_field(
                        ax, x_interp, energy_interp,
                        linestyle='-', alpha=0.7,
                        label="Interpolated (shape functions)" if elem_id == element_ids[0] else None
                    )
                except Exception:
                    # Skip elements that fail
                    continue

        # Plot nodal markers (B2: projected = hollow circle)
        plot_nodal_points(
            ax, node_positions.reshape(-1, 1), energy_density,
            color=self._BLUE, size=NODAL_MARKER_SIZE, alpha=0.9,
            nodal_data_type='projected',
            label="Nodes"
        )

        # Plot Gauss point markers (B2: small solid circle) when gaussian data is available
        if has_gauss:
            plot_gauss_points(
                ax, x_gauss, energy_gauss,
                color="red", size=GAUSS_MARKER_SIZE, alpha=0.9,
                label="Gauss points",
            )

        # Beam-end anchors + baseline
        ax.plot([x_min], [0], marker="o", color="k", ms=3, zorder=3)
        ax.plot([x_max], [0], marker="o", color="k", ms=3, zorder=3)
        ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
        ax.grid(ls="--", alpha=0.6)
        ax.set_ylabel(r"$w$ [J/m³]")
        ax.set_xlabel(r"$x$ [m]")
        
        # Add unified legend at bottom of figure (B2 convention; only resolution levels present)
        from matplotlib.lines import Line2D
        legend_elements = []
        if has_elements:
            legend_elements.append(
                Line2D([0], [0], linestyle='-', linewidth=INTERPOLANT_LINEWIDTH, color='black',
                       label='Interpolated (shape functions)'))
        legend_elements.append(
            Line2D([0], [0], marker='o', linestyle='None', markersize=LEGEND_MARKER_SIZE,
                   markerfacecolor='none', markeredgecolor=self._BLUE, markeredgewidth=1.0,
                   label='Nodes'))
        if has_gauss:
            legend_elements.append(
                Line2D([0], [0], marker='o', linestyle='None', markersize=LEGEND_MARKER_SIZE_SECONDARY,
                       color='red', label='Gauss points'))
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
    def _read_nodal_energy_density(file: Path) -> Optional[np.ndarray]:
        """Read nodal energy density CSV file."""
        try:
            energy = np.genfromtxt(file, delimiter=",", skip_header=1)
            if energy.ndim == 0:
                energy = np.array([energy])
            elif energy.ndim > 1:
                energy = energy.flatten()
            return energy
        except Exception as exc:
            print(f"Error reading {file}: {exc}")
            return None

    def _load_gaussian_energy_density(
        self,
        job_dir: Path,
        element_dictionary: dict,
        grid_dictionary: dict,
    ) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Load strain energy density at Gauss points and compute their x-coordinates.
        Returns (x_gauss, energy_gauss) or (None, None) if not available.
        """
        energy_dir = job_dir / "secondary_results" / "gaussian" / "energy_density"
        if not energy_dir.is_dir():
            return None, None
        files = sorted(
            energy_dir.glob("energy_density_elem_*.csv"),
            key=lambda p: int(p.stem.split("_")[-1]),
        )
        if not files:
            return None, None
        ids = element_dictionary.get("ids", np.arange(len(files)))
        if len(ids) != len(files):
            return None, None
        x_list = []
        energy_list = []
        for elem_idx, csv_path in enumerate(files):
            try:
                data = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
            except Exception:
                continue
            data = np.atleast_1d(data).flatten()
            elem_id = int(ids[elem_idx]) if hasattr(ids, "__getitem__") else elem_idx
            try:
                node_coords = get_element_node_coords(
                    elem_id, element_dictionary, grid_dictionary
                )
            except Exception:
                continue
            n_gp = len(data)
            xi = (
                np.polynomial.legendre.leggauss(3)[0]
                if n_gp == 3
                else np.polynomial.legendre.leggauss(n_gp)[0]
            )
            x_gp = natural_to_physical_coords(xi, node_coords)
            x_list.append(x_gp)
            energy_list.append(data)
        if not x_list:
            return None, None
        x_gauss = np.concatenate(x_list)
        energy_gauss = np.concatenate(energy_list)
        return x_gauss, energy_gauss

    # ------------------------------------------------------------------#
    #  Driver
    # ------------------------------------------------------------------#
    def process_all(self) -> None:
        pattern = str(self.results_dir / "job_*" / "secondary_results" / "nodal" / "nodal_strain_energy_density.csv")
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No energy density files found.")
            return

        for csv_path in csv_files:
            csv_file = Path(csv_path)
            job_dir = csv_file.parent.parent.parent
            parsed = parse_job_result_dir_name(job_dir.name)
            if not parsed:
                print(f"Skipping unrecognised folder '{job_dir.name}'")
                continue

            job_id, job_folder_name, timestamp = parsed
            grid_file = self.jobs_dir / job_folder_name / "grid.txt"
            element_file = self.jobs_dir / job_folder_name / "element.txt"

            print(f"> Processing job {job_id} ({timestamp})")

            # ---- Energy density ----------------------------------------- #
            energy_density = self._read_nodal_energy_density(csv_file)
            if energy_density is None:
                print(f"WARNING: Could not read energy density for job {job_id}, skipping.")
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

            # ---- Gaussian energy density (optional) ------------------------ #
            x_gauss, energy_gauss = None, None
            if element_dictionary and grid_dictionary:
                x_gauss, energy_gauss = self._load_gaussian_energy_density(
                    job_dir, element_dictionary, grid_dictionary
                )

            # ---- Plot ---------------------------------------------------- #
            fig_name = f"energy_density_{job_folder_name}_{timestamp}.png"
            self._plot(
                energy_density,
                node_coords[:, 0],
                element_dictionary=element_dictionary if element_dictionary else {},
                grid_dictionary=grid_dictionary if grid_dictionary else {},
                x_gauss=x_gauss,
                energy_gauss=energy_gauss,
                title_suffix=f"{job_folder_name}_{timestamp}",
                save_path=self.figure_output_dir / fig_name,
            )


if __name__ == "__main__":
    VisualiseEnergyDensity().process_all()

