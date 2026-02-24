"""
Principal stress visualisation (tertiary results, Gaussian resolution).

Reads job_*/tertiary_results/gaussian/principal_stress/principal_stress_elem_*.csv
and plots σ1, σ2, σ3 vs position along structure.
- B2: red small solid circle at Gauss points; projected nodal = hollow circle;
  interpolated = thin black line (shape functions).
GP positions from element geometry and 3-point Gauss-Legendre rule.
Nodal values extrapolated from Gauss point data for consistency with other plots.
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

from pre_processing.parsing.grid_parser import GridParser  # type: ignore
from pre_processing.parsing.element_parser import ElementParser  # type: ignore
from post_processing.graphical_visualisers.job_discovery_utils import parse_job_result_dir_name
from post_processing.graphical_visualisers.resolution_plotting_utils import (
    get_element_node_coords,
    natural_to_physical_coords,
    plot_nodal_points,
    plot_gauss_points,
    plot_interpolated_field,
    interpolate_field_nodal_to_continuous,
    INTERPOLANT_LINEWIDTH,
    NODAL_MARKER_SIZE,
    LEGEND_MARKER_SIZE,
    LEGEND_MARKER_SIZE_SECONDARY,
    GAUSS_MARKER_SIZE,
)

# 3-point Gauss-Legendre (same as EB 3D default quadrature)
GAUSS_3PT_XI: Final[np.ndarray] = np.polynomial.legendre.leggauss(3)[0]

COMPONENTS = [r"$\sigma_1$", r"$\sigma_2$", r"$\sigma_3$"]


def _nodal_shape_matrix_at_xi(xi: np.ndarray, n_nodes: int) -> np.ndarray:
    """Shape function matrix N at natural coords xi. Shape (len(xi), n_nodes)."""
    if n_nodes == 2:
        N = np.zeros((len(xi), 2))
        N[:, 0] = (1 - xi) / 2
        N[:, 1] = (1 + xi) / 2
        return N
    if n_nodes == 3:
        N = np.zeros((len(xi), 3))
        N[:, 0] = xi * (xi - 1) / 2
        N[:, 1] = 1 - xi**2
        N[:, 2] = xi * (xi + 1) / 2
        return N
    raise ValueError(f"Unsupported n_nodes={n_nodes}")


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


def _gauss_point_x_for_element(
    node_coords: np.ndarray,
    xi: np.ndarray = GAUSS_3PT_XI,
) -> np.ndarray:
    """Physical x-coordinates of Gauss points for one element."""
    return natural_to_physical_coords(xi, node_coords)


class VisualisePrincipalStress:
    """Produce principal stress from tertiary_results (Gauss + projected nodal + interpolated)."""

    _BLUE: Final[str] = "#4F81BD"

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "principal_stress_plots"
        self.figure_output_dir.mkdir(exist_ok=True)
        self.results_dir: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
        self.jobs_dir: Final[Path] = PROJECT_ROOT / "jobs"

    def _plot(
        self,
        x_gauss: np.ndarray,
        stress_gauss: np.ndarray,
        node_positions: np.ndarray,
        stress_nodal: np.ndarray,
        element_dictionary: dict,
        grid_dictionary: dict,
        *,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        # stress_gauss (n_gp, 3), stress_nodal (n_nodes, 3)
        has_elements = (
            element_dictionary
            and "ids" in element_dictionary
            and "connectivity" in element_dictionary
            and grid_dictionary
            and "coordinates" in grid_dictionary
        )
        has_gauss = x_gauss is not None and len(x_gauss) > 0

        fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
        for k, (ax, name) in enumerate(zip(axes, COMPONENTS)):
            if has_elements:
                element_ids = element_dictionary["ids"]
                for elem_id in element_ids:
                    try:
                        elem_idx = int(np.where(element_dictionary["ids"] == elem_id)[0][0])
                        node_ids = element_dictionary["connectivity"][elem_idx]
                        elem_vals = stress_nodal[node_ids, k]
                        x_interp, val_interp = interpolate_field_nodal_to_continuous(
                            nodal_values=elem_vals,
                            element_id=elem_id,
                            element_dictionary=element_dictionary,
                            grid_dictionary=grid_dictionary,
                            n_points=50,
                        )
                        plot_interpolated_field(
                            ax, x_interp, val_interp,
                            linestyle="-", alpha=0.7,
                            label="Interpolated (shape functions)" if k == 0 and elem_id == element_ids[0] else None,
                        )
                    except Exception:
                        continue

            plot_nodal_points(
                ax, node_positions.reshape(-1, 1), stress_nodal[:, k],
                color=self._BLUE, size=NODAL_MARKER_SIZE, alpha=0.9,
                nodal_data_type="projected",
                label="Nodes" if k == 0 else None,
            )

            plot_gauss_points(
                ax, x_gauss, stress_gauss[:, k],
                color="red", size=GAUSS_MARKER_SIZE, alpha=0.9,
                label="Gauss points" if k == 0 else None,
            )

            ax.set_ylabel(name)
            ax.grid(ls="--", alpha=0.6)
            ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
        axes[-1].set_xlabel(r"$x$ [m]")
        fig.suptitle(f"Principal stresses {title_suffix}".strip(), fontweight="bold")

        # Legend: only resolution levels present (B2 convention)
        from matplotlib.lines import Line2D
        legend_elements = []
        if has_elements:
            legend_elements.append(
                Line2D([0], [0], linestyle="-", linewidth=INTERPOLANT_LINEWIDTH, color="black",
                       label="Interpolated (shape functions)"))
        legend_elements.append(
            Line2D([0], [0], marker="o", linestyle="None", markersize=LEGEND_MARKER_SIZE,
                   markerfacecolor="none", markeredgecolor=self._BLUE, markeredgewidth=1.0,
                   label="Nodes"))
        if has_gauss:
            legend_elements.append(
                Line2D([0], [0], marker="o", linestyle="None", markersize=LEGEND_MARKER_SIZE_SECONDARY,
                       color="red", label="Gauss points"))
        if legend_elements:
            fig.legend(handles=legend_elements, loc="lower center", ncol=len(legend_elements),
                      fontsize=9, frameon=True, bbox_to_anchor=(0.5, 0.06))

        fig.tight_layout()
        fig.subplots_adjust(top=0.9, bottom=0.14 if legend_elements else 0.1)
        if save_path:
            fig.savefig(save_path, dpi=300)
            plt.close(fig)
        else:
            plt.show()

    def process_all(self) -> None:
        pattern = str(self.results_dir / "job_*" / "tertiary_results" / "gaussian" / "principal_stress" / "principal_stress_elem_*.csv")
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No principal_stress files found.")
            return

        by_job: dict[str, list[Path]] = {}
        for p in csv_files:
            path = Path(p)
            job_dir = path.parent.parent.parent.parent  # principal_stress -> gaussian -> tertiary_results -> job_XXX
            key = job_dir.name
            by_job.setdefault(key, []).append(path)

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

            grid = GridParser(str(grid_file), str(job_dir)).parse()
            try:
                _get_node_coordinates(grid)
            except (AttributeError, KeyError) as exc:
                print(f"WARNING: {exc} for job {job_id}, skipping.")
                continue

            elem_parsed = ElementParser(str(element_file), str(job_dir)).parse()
            element_dictionary = elem_parsed["element_dictionary"]
            grid_dictionary = grid["grid_dictionary"] if isinstance(grid, dict) else {}

            ids = element_dictionary.get("ids", np.arange(len(files_sorted)))
            if len(ids) != len(files_sorted):
                print(f"Element count mismatch for job {job_id}, skipping.")
                continue

            coords = grid_dictionary.get("coordinates")
            n_nodes = coords.shape[0] if coords is not None else 0
            node_positions = coords[:, 0] if coords is not None else np.array([])
            stress_nodal = np.zeros((n_nodes, 3))
            weight = np.zeros(n_nodes)

            x_list = []
            stress_list = []
            for i, csv_path in enumerate(files_sorted):
                elem_id = int(ids[i]) if hasattr(ids, "__getitem__") else i
                node_ids = element_dictionary["connectivity"][i]
                n_nodes_elem = len(node_ids)
                try:
                    node_coords = get_element_node_coords(
                        elem_id, element_dictionary, grid_dictionary
                    )
                except Exception:
                    continue
                try:
                    data = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
                except Exception:
                    continue
                if data.ndim == 1:
                    data = data.reshape(1, -1)
                if data.shape[1] != 3:
                    continue
                n_gp = data.shape[0]
                xi_used = GAUSS_3PT_XI if n_gp == 3 else np.polynomial.legendre.leggauss(n_gp)[0]
                x_gp = _gauss_point_x_for_element(node_coords, xi_used)
                x_list.append(x_gp)
                stress_list.append(data)

                if n_nodes_elem not in (2, 3):
                    continue
                N_mat = _nodal_shape_matrix_at_xi(xi_used, n_nodes_elem)
                if n_gp == n_nodes_elem:
                    nodal_elem = np.linalg.solve(N_mat, data)
                else:
                    nodal_elem = np.linalg.lstsq(N_mat, data, rcond=None)[0]
                for node_idx in range(n_nodes_elem):
                    nid = node_ids[node_idx]
                    if nid < n_nodes:
                        stress_nodal[nid] += nodal_elem[node_idx]
                        weight[nid] += 1.0

            if not x_list:
                print(f"No valid principal stress data for job {job_id}, skipping.")
                continue

            nonzero = weight > 0
            if np.any(nonzero):
                stress_nodal[nonzero] /= weight[nonzero, np.newaxis]

            x_gauss = np.concatenate(x_list)
            stress_gauss = np.vstack(stress_list)

            fig_name = f"principal_stress_{job_folder_name}_{timestamp}.png"
            save_path = self.figure_output_dir / fig_name
            self._plot(
                x_gauss,
                stress_gauss,
                node_positions,
                stress_nodal,
                element_dictionary,
                grid_dictionary,
                title_suffix=f"{job_folder_name}_{timestamp}",
                save_path=save_path,
            )
            print(f"Saved: {save_path}")


if __name__ == "__main__":
    VisualisePrincipalStress().process_all()
