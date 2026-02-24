"""
Stress visualisation utility.

Reads stress results from secondary_results and produces stress resultant profiles.
For beam elements the stored "stress" is the 6 resultants [N, Vy, Vz, T, My, Mz].
The bending stress normal to the section (σ_xx) is recovered from N, My, Mz and
section properties; here we plot the resultants so that e.g. transverse load in y
activates Vy and Mz.

- Nodal markers at node locations (B2: projected = hollow circle)
- Gauss point markers when gaussian stress is saved (B2: small solid circle)
- Continuous interpolated fields using element shape functions

Components: [N, Vy, Vz, T, My, Mz] (axial, shear y/z, torsion, bending y/z).

Euler-Bernoulli elements yield Vy = Vz = 0 (no shear strain); shear-deformable
elements (e.g. Timoshenko, Levinson) yield non-zero shear resultants.
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
        x_gauss: Optional[np.ndarray] = None,
        stress_gauss: Optional[np.ndarray] = None,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        """
        Plot stress component profiles with nodal markers, optional Gauss point markers, and interpolated continuous fields.
        
        Parameters
        ----------
        stress : np.ndarray
            Stress resultants, shape (n_nodes, 6) with columns [N, Vy, Vz, T, My, Mz]
        node_positions : np.ndarray
            Node x-coordinates, shape (n_nodes,)
        element_dictionary : dict
            Element dictionary with connectivity and types
        grid_dictionary : dict
            Grid dictionary with node coordinates
        x_gauss : Optional[np.ndarray]
            Gauss point x-coordinates (when gaussian data is loaded)
        stress_gauss : Optional[np.ndarray]
            Stress at Gauss points, shape (n_gauss, 6)
        title_suffix : str
            Suffix for plot title
        save_path : Optional[Path]
            Path to save figure
        """
        if stress.shape[1] != 6:
            raise ValueError("Stress must be shaped (n_nodes, 6)")

        x_min, x_max = float(node_positions.min()), float(node_positions.max())
        has_gauss = (
            x_gauss is not None and stress_gauss is not None
            and stress_gauss.shape[1] == 6
            and len(x_gauss) == len(stress_gauss) and len(x_gauss) > 0
        )

        fig, axes = plt.subplots(3, 2, figsize=(15, 10), sharex=True)
        fig.suptitle(
            rf"Stress profiles${' – ' + title_suffix if title_suffix else ''}",
            fontsize=16,
            fontweight="bold",
        )

        # Beam stress resultants: (left_idx, left_label, right_idx, right_label)
        # Order in data: [N, Vy, Vz, T, My, Mz]
        pairs = [
            (0, r"$N$ [N]", 1, r"$V_y$ [N]"),
            (2, r"$V_z$ [N]", 3, r"$T$ [N·m]"),
            (4, r"$M_y$ [N·m]", 5, r"$M_z$ [N·m]"),
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
                        
                        # Plot interpolated fields (B2: thin black line)
                        plot_interpolated_field(
                            ax_l, x_interp_norm, norm_interp,
                            linestyle='-', alpha=0.7,
                            label="Interpolated (shape functions)" if i == 0 and elem_id == element_ids[0] else None
                        )
                        plot_interpolated_field(
                            ax_r, x_interp_shear, shear_interp,
                            linestyle='-', alpha=0.7, label=None
                        )
                    except Exception:
                        # Skip elements that fail
                        continue

            # Plot nodal markers (B2: projected = hollow circle)
            plot_nodal_points(
                ax_l, node_positions.reshape(-1, 1), norm_stress,
                color=self._BLUE, size=NODAL_MARKER_SIZE, alpha=0.9,
                nodal_data_type='projected',
                label="Nodes" if i == 0 else None
            )
            plot_nodal_points(
                ax_r, node_positions.reshape(-1, 1), shear_stress,
                color=self._BLUE, size=NODAL_MARKER_SIZE, alpha=0.9,
                nodal_data_type='projected',
                label=None  # Only label on left subplot
            )

            # Plot Gauss point markers (B2: small solid circle) when available
            if has_gauss:
                norm_g = stress_gauss[:, norm_idx]
                shear_g = stress_gauss[:, shear_idx]
                plot_gauss_points(
                    ax_l, x_gauss, norm_g,
                    color="red", size=GAUSS_MARKER_SIZE, alpha=0.9,
                    label="Gauss points" if i == 0 else None
                )
                plot_gauss_points(
                    ax_r, x_gauss, shear_g,
                    color="red", size=GAUSS_MARKER_SIZE, alpha=0.9,
                    label=None
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
                ax_l.set_title(r"Resultants $N$, $V_z$, $M_y$", fontweight="bold")
                ax_r.set_title(r"Resultants $V_y$, $T$, $M_z$", fontweight="bold")

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
    _RESULTANT_ORDER = (0, 3, 4, 5, 1, 2)  # formulation [N,My,Mz,Vy,Vz,T] -> [N,Vy,Vz,T,My,Mz]

    @classmethod
    def _read_nodal_stress(cls, file: Path) -> Optional[np.ndarray]:
        """Read nodal stress CSV. If header is old (σ_xx, ...), reorder to resultants [N,Vy,Vz,T,My,Mz]."""
        try:
            with open(file, encoding="utf-8") as f:
                first_line = f.readline().strip()
            stress = np.genfromtxt(file, delimiter=",", skip_header=1)
            if stress.ndim == 1:
                stress = stress.reshape(-1, 6)
            if stress.shape[1] != 6:
                raise ValueError("Stress must have 6 components")
            # Old CSVs used 3D Voigt header; data was in formulation order
            if "σ_xx" in first_line or "sigma" in first_line.lower():
                stress = stress[:, cls._RESULTANT_ORDER]
            return stress
        except Exception as exc:
            print(f"Error reading {file}: {exc}")
            return None

    def _load_gaussian_stress(
        self,
        job_dir: Path,
        element_dictionary: dict,
        grid_dictionary: dict,
    ) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Load stress at Gauss points and compute their x-coordinates.
        Returns (x_gauss, stress_gauss) or (None, None) if not available.
        """
        stress_dir = job_dir / "secondary_results" / "gaussian" / "stress"
        if not stress_dir.is_dir():
            return None, None
        files = sorted(
            stress_dir.glob("stress_elem_*.csv"),
            key=lambda p: int(p.stem.split("_")[-1]),
        )
        if not files:
            return None, None
        ids = element_dictionary.get("ids", np.arange(len(files)))
        if len(ids) != len(files):
            return None, None
        x_list = []
        stress_list = []
        for elem_idx, csv_path in enumerate(files):
            try:
                with open(csv_path, encoding="utf-8") as f:
                    first_line = f.readline().strip()
                data = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
            except Exception:
                continue
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] != 6:
                continue
            if "σ_xx" in first_line or "sigma" in first_line.lower():
                data = data[:, self._RESULTANT_ORDER]
            elem_id = int(ids[elem_idx]) if hasattr(ids, "__getitem__") else elem_idx
            try:
                node_coords = get_element_node_coords(
                    elem_id, element_dictionary, grid_dictionary
                )
            except Exception:
                continue
            n_gp = data.shape[0]
            xi = (
                np.polynomial.legendre.leggauss(3)[0]
                if n_gp == 3
                else np.polynomial.legendre.leggauss(n_gp)[0]
            )
            x_gp = natural_to_physical_coords(xi, node_coords)
            x_list.append(x_gp)
            stress_list.append(data)
        if not x_list:
            return None, None
        x_gauss = np.concatenate(x_list)
        stress_gauss = np.vstack(stress_list)
        return x_gauss, stress_gauss

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
            parsed = parse_job_result_dir_name(job_dir.name)
            if not parsed:
                print(f"Skipping unrecognised folder '{job_dir.name}'")
                continue

            job_id, job_folder_name, timestamp = parsed
            grid_file = self.jobs_dir / job_folder_name / "grid.txt"
            element_file = self.jobs_dir / job_folder_name / "element.txt"

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

            # ---- Gaussian stress (optional) ------------------------------ #
            x_gauss, stress_gauss = None, None
            if element_dictionary and grid_dictionary:
                x_gauss, stress_gauss = self._load_gaussian_stress(
                    job_dir, element_dictionary, grid_dictionary
                )

            # ---- Plot ---------------------------------------------------- #
            fig_name = f"stress_{job_folder_name}_{timestamp}.png"
            self._plot(
                stress,
                node_coords[:, 0],
                element_dictionary=element_dictionary if element_dictionary else {},
                grid_dictionary=grid_dictionary if grid_dictionary else {},
                x_gauss=x_gauss,
                stress_gauss=stress_gauss,
                title_suffix=f"{job_folder_name}_{timestamp}",
                save_path=self.figure_output_dir / fig_name,
            )


if __name__ == "__main__":
    VisualiseStress().process_all()

