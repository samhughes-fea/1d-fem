"""
Section forces visualisation (tertiary results, Gaussian resolution).

Reads job_*/tertiary_results/gaussian/section_forces/section_forces_elem_*.csv
and plots N, Vy, Vz, T, My, Mz vs position along structure.

Column order: CSV and plot use [N, Vy, Vz, T, My, Mz]. The formulation (B, D)
outputs [N, M_y, M_z, V_y, V_z, T]; section force computation reorders before
save. Legacy CSVs (no "# column_order=resultant" first line) were saved in
formulation order; we reorder on read so existing plots correct without re-running.

Gauss points: Physical plot positions use natural coordinates **ξ** per CSV row.
When the saver writes ``# xi_per_row=...`` (ascending ξ, same order as GP rows),
those values are used. Otherwise ξ defaults to Gauss–Legendre points for ``n_gp``
(3-point: ξ ≈ −0.77, 0, 0.77). See docs/proofs/timoshenko/gauss_points_and_reduced_integration.md.

**Nodal markers:** If ``tertiary_results/nodal/nodal_section_forces.csv`` exists
(pipeline ``NodalSectionForcesProjector``: formulation-cache shape matrix,
least-squares / solve where applicable, **element mean at boundary nodes**), those
values are used. If it is missing, nodal markers fall back to the **element-wise
mean of GP** section forces averaged over elements meeting each node.

Euler-Bernoulli: Vy and Vz are zero (no shear strain); shear in EB comes from
equilibrium (V = dM/dx). Timoshenko/Levinson elements produce non-zero Vy, Vz.

Plot elements: red solid circles at Gauss points; hollow circles at nodal samples;
thin interpolated line where applicable.
"""

from __future__ import annotations

import glob
import os
import sys
import warnings
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

COMPONENTS = ["N", "Vy", "Vz", "T", "My", "Mz"]

# Legacy CSVs (no format marker) had formulation order [N, M_y, M_z, V_y, V_z, T]; reorder to [N,Vy,Vz,T,My,Mz]
_FORMULATION_TO_RESULTANT = (0, 3, 4, 5, 1, 2)


def _read_elem_section_forces_csv(
    csv_path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load one element section-forces CSV.

    Returns
    -------
    data : ndarray, shape (n_gp, 6)
        Resultant order [N, Vy, Vz, T, My, Mz].
    xi_used : ndarray, shape (n_gp,)
        Natural coordinates for each row: from ``# xi_per_row=`` when present
        and length matches ``n_gp``; else Gauss–Legendre points (3-point shortcut
        when ``n_gp == 3``).
    """
    with open(csv_path, encoding="utf-8") as f:
        first_line = f.readline()
    skip = 2 if "column_order=resultant" in first_line else 1
    xi_from_file: Optional[np.ndarray] = None
    if skip == 2:
        with open(csv_path, encoding="utf-8") as f:
            f.readline()
            second_line = f.readline()
        if second_line.strip().startswith("# xi_per_row="):
            rest = second_line.split("=", 1)[1].strip()
            xi_from_file = np.array([float(x) for x in rest.split(",")], dtype=np.float64)
            skip = 3

    data = np.genfromtxt(csv_path, delimiter=",", skip_header=skip)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 6:
        raise ValueError(f"Expected 6 columns, got {data.shape[1]} in {csv_path}")
    if skip == 1:
        data = np.asarray(data[:, _FORMULATION_TO_RESULTANT], dtype=np.float64)

    n_gp = int(data.shape[0])
    if xi_from_file is not None and xi_from_file.size == n_gp:
        xi_used = xi_from_file
    else:
        if xi_from_file is not None and xi_from_file.size != n_gp:
            warnings.warn(
                f"{csv_path.name}: xi_per_row length {xi_from_file.size} != n_gp {n_gp}; "
                "using default Gauss–Legendre ξ.",
                stacklevel=2,
            )
        xi_used = GAUSS_3PT_XI if n_gp == 3 else np.polynomial.legendre.leggauss(n_gp)[0]
    return data, xi_used


def _nodal_shape_matrix_at_xi(xi: np.ndarray, n_nodes: int) -> np.ndarray:
    """Shape function matrix N at natural coords xi. Shape (len(xi), n_nodes)."""
    if n_nodes == 2:
        # Linear: N1 = (1-xi)/2, N2 = (1+xi)/2
        N = np.zeros((len(xi), 2))
        N[:, 0] = (1 - xi) / 2
        N[:, 1] = (1 + xi) / 2
        return N
    if n_nodes == 3:
        # Quadratic: N1 = xi*(xi-1)/2, N2 = 1-xi^2, N3 = xi*(xi+1)/2
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


class VisualiseSectionForces:
    """Produce section forces from tertiary_results (Gauss + projected nodal + interpolated)."""

    _BLUE: Final[str] = "#4F81BD"

    def __init__(self) -> None:
        self.figure_output_dir: Final[Path] = SCRIPT_DIR / "section_forces_plots"
        self.figure_output_dir.mkdir(exist_ok=True)
        self.results_dir: Final[Path] = PROJECT_ROOT / "post_processing" / "results"
        self.jobs_dir: Final[Path] = PROJECT_ROOT / "jobs"

    def _plot(
        self,
        x_gauss: np.ndarray,
        forces_gauss: np.ndarray,
        node_positions: np.ndarray,
        forces_nodal: np.ndarray,
        element_dictionary: dict,
        grid_dictionary: dict,
        *,
        title_suffix: str = "",
        save_path: Optional[Path] = None,
    ) -> None:
        # forces_gauss (n_gp, 6), forces_nodal (n_nodes, 6)
        has_elements = (
            element_dictionary
            and "ids" in element_dictionary
            and "connectivity" in element_dictionary
            and grid_dictionary
            and "coordinates" in grid_dictionary
        )
        has_gauss = x_gauss is not None and len(x_gauss) > 0

        fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
        axes = axes.ravel()
        for k, (ax, name) in enumerate(zip(axes, COMPONENTS)):
            # Interpolated field per element (B2: thin black line)
            if has_elements:
                element_ids = element_dictionary["ids"]
                for elem_id in element_ids:
                    try:
                        elem_idx = int(np.where(element_dictionary["ids"] == elem_id)[0][0])
                        node_ids = element_dictionary["connectivity"][elem_idx]
                        elem_vals = forces_nodal[node_ids, k]
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

            # Nodal markers (B2: projected = hollow circle)
            plot_nodal_points(
                ax, node_positions.reshape(-1, 1), forces_nodal[:, k],
                color=self._BLUE, size=NODAL_MARKER_SIZE, alpha=0.9,
                nodal_data_type="projected",
                label="Nodes" if k == 0 else None,
            )

            # Gauss point markers (B2: red small solid circle)
            plot_gauss_points(
                ax, x_gauss, forces_gauss[:, k],
                color="red", size=GAUSS_MARKER_SIZE, alpha=0.9,
                label="Gauss points" if k == 0 else None,
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
        fig.suptitle(f"Section forces {title_suffix}".strip(), fontweight="bold")

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
        # Discover jobs that have section_forces dir
        pattern = str(self.results_dir / "job_*" / "tertiary_results" / "gaussian" / "section_forces" / "section_forces_elem_*.csv")
        csv_files = sorted(glob.glob(pattern))
        if not csv_files:
            print("No section_forces files found.")
            return

        # Group by job (job_dir is parent of gaussian/section_forces)
        by_job: dict[str, list[Path]] = {}
        for p in csv_files:
            path = Path(p)
            # .../job_XXX_ts/tertiary_results/gaussian/section_forces/section_forces_elem_000000.csv
            job_dir = path.parent.parent.parent.parent  # section_forces -> gaussian -> tertiary_results -> job_XXX
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

            # Sort by element index
            def elem_index(path: Path) -> int:
                stem = path.stem  # section_forces_elem_000000
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
            forces_nodal = np.zeros((n_nodes, 6))
            weight = np.zeros(n_nodes)

            # Prefer pipeline-projected nodal section forces when available (Option A)
            nodal_csv = job_dir / "tertiary_results" / "nodal" / "nodal_section_forces.csv"
            if nodal_csv.is_file():
                try:
                    with open(nodal_csv, encoding="utf-8") as f:
                        first_line = f.readline()
                    skip = 2 if "column_order=resultant" in first_line else 1
                    nodal_data = np.genfromtxt(nodal_csv, delimiter=",", skip_header=skip)
                    if nodal_data.ndim == 1:
                        nodal_data = nodal_data.reshape(1, -1)
                    if nodal_data.shape[0] == n_nodes and nodal_data.shape[1] == 6:
                        forces_nodal = nodal_data.copy()
                        weight = np.ones(n_nodes)  # mark as filled so we skip GP-based nodal below
                except Exception:
                    pass

            x_list = []
            forces_list = []
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
                    data, xi_used = _read_elem_section_forces_csv(Path(csv_path))
                except Exception:
                    continue
                x_gp = _gauss_point_x_for_element(node_coords, xi_used)
                x_list.append(x_gp)
                forces_list.append(data)

                if n_nodes_elem not in (2, 3):
                    continue
                # Nodal value = mean of GP values in this element (fallback when no nodal_section_forces.csv)
                if np.any(weight == 0):
                    elem_mean = data.mean(axis=0)
                    for nid in node_ids:
                        if nid < n_nodes:
                            forces_nodal[nid] += elem_mean
                            weight[nid] += 1.0

            if not x_list:
                print(f"No valid section force data for job {job_id}, skipping.")
                continue

            nonzero = weight > 0
            if np.any(nonzero):
                forces_nodal[nonzero] /= weight[nonzero, np.newaxis]

            x_gauss = np.concatenate(x_list)
            forces_gauss = np.vstack(forces_list)

            # Optional diagnostic: print min/max per component (e.g. for job_0003 expect |Vy| ~ 500 N)
            if os.environ.get("DEBUG_SECTION_FORCES"):
                print(f"  [job_{job_id}] Section force ranges (N or N·m):")
                for k, name in enumerate(COMPONENTS):
                    col = forces_gauss[:, k]
                    print(f"    {name}: min={col.min():.3g} max={col.max():.3g}")

            fig_name = f"section_forces_{job_folder_name}_{timestamp}.png"
            save_path = self.figure_output_dir / fig_name
            self._plot(
                x_gauss,
                forces_gauss,
                node_positions,
                forces_nodal,
                element_dictionary,
                grid_dictionary,
                title_suffix=f"{job_folder_name}_{timestamp}",
                save_path=save_path,
            )
            print(f"Saved: {save_path}")


if __name__ == "__main__":
    VisualiseSectionForces().process_all()
