# processing_OOP\static\results\compute_tertiary\nodal_section_forces_projector.py

"""
Nodal Section Forces Projector

Projects section force resultants [N, Vy, Vz, T, My, Mz] from Gauss points to nodes
using shape function extrapolation (same N_mat as NodalResultProjector). At boundary
nodes (one contributing element), uses the element mean instead of the extrapolate
to avoid spikes.
"""

import numpy as np
from typing import List, Dict, Any
import logging

logger = logging.getLogger(__name__)


class NodalSectionForcesProjector:
    """
    Projects Gaussian section forces to nodal resolution using the formulation
    cache shape matrix. Boundary nodes get the element mean; interior nodes get
    the average of extrapolated values from adjacent elements.
    """

    def __init__(
        self,
        section_forces_gauss: List[List[np.ndarray]],
        formulation_cache: Any,
        element_dictionary: Dict[str, Any],
        grid_dictionary: Dict[str, Any],
    ):
        """
        Parameters
        ----------
        section_forces_gauss : List[List[np.ndarray]]
            Section forces at Gauss points, same order as formulation_cache elements.
            Shape: List[element] -> List[gauss_point] -> np.ndarray(6,)
        formulation_cache
            Cached formulation data with element_objects and gauss_data
        element_dictionary : Dict
            Element connectivity and metadata (ids, connectivity)
        grid_dictionary : Dict
            Node/grid data (coordinates for n_nodes)
        """
        self.section_forces_gauss = section_forces_gauss
        self.formulation_cache = formulation_cache
        self.element_dictionary = element_dictionary
        self.grid_dictionary = grid_dictionary

        # Map element_id -> index into section_forces_gauss (same order as element_objects)
        self.elem_id_to_idx = {}
        for idx, elem_obj in enumerate(self.formulation_cache.element_objects):
            self.elem_id_to_idx[elem_obj.element_id] = idx

        # Number of nodes
        if "coordinates" in grid_dictionary:
            self.n_nodes = len(grid_dictionary["coordinates"])
        elif "grid_dictionary" in grid_dictionary and "coordinates" in grid_dictionary.get("grid_dictionary", {}):
            self.n_nodes = len(grid_dictionary["grid_dictionary"]["coordinates"])
        else:
            raise ValueError("Cannot determine number of nodes from grid_dictionary")

    def _element_has_shape_functions(self, elem_obj: Any, n_nodes_elem: int) -> bool:
        """True if every GP has shape_functions and n_gauss >= n_nodes_elem."""
        n_gauss = len(elem_obj.gauss_data)
        if n_gauss < n_nodes_elem:
            return False
        return all(gp.shape_functions is not None for gp in elem_obj.gauss_data)

    def _build_nodal_shape_matrix(self, elem_obj: Any, n_nodes_elem: int) -> np.ndarray:
        """
        Build (n_gauss, n_nodes_elem) matrix N such that values_g = N @ values_nodal.
        Mirrors NodalResultProjector logic for 1D or beam-like shape_functions.
        """
        rows = []
        for gp in elem_obj.gauss_data:
            sf = np.asarray(gp.shape_functions).squeeze()
            if sf.ndim == 1 and sf.size == n_nodes_elem:
                rows.append(sf.astype(float))
            elif sf.ndim >= 2 and sf.shape[0] >= n_nodes_elem:
                dofs_per_node = sf.shape[0] // n_nodes_elem
                row = np.array(
                    [sf[i * dofs_per_node, 0] for i in range(n_nodes_elem)],
                    dtype=float,
                )
                rows.append(row)
            else:
                flat = sf.flatten()
                rows.append(flat[:n_nodes_elem].astype(float))
        return np.array(rows)

    def _get_element_node_ids(self, elem_id: int) -> List[int]:
        """Get node IDs for an element by element_id."""
        elem_idx = np.where(self.element_dictionary["ids"] == elem_id)[0][0]
        node_ids = self.element_dictionary["connectivity"][elem_idx]
        return node_ids.tolist()

    def project(self) -> np.ndarray:
        """
        Project section forces from Gauss points to nodes.

        Returns
        -------
        np.ndarray
            Shape (n_nodes, 6), columns [N, Vy, Vz, T, My, Mz].
            Boundary nodes (one element) get the element mean; interior nodes
            get the average of extrapolated values from adjacent elements.
        """
        nodal_section_forces = np.zeros((self.n_nodes, 6))
        nodal_weight = np.zeros(self.n_nodes)
        # For boundary nodes (weight == 1), we replace with element mean
        nodal_elem_mean = np.full((self.n_nodes, 6), np.nan)

        for elem_obj in self.formulation_cache.element_objects:
            elem_id = elem_obj.element_id
            if elem_id not in self.elem_id_to_idx:
                logger.warning(
                    "Element %s not found in section_forces mapping, skipping",
                    elem_id,
                )
                continue

            gauss_idx = self.elem_id_to_idx[elem_id]
            node_ids = self._get_element_node_ids(elem_id)
            n_nodes_elem = len(node_ids)
            elem_forces = self.section_forces_gauss[gauss_idx]
            section_g = np.array(elem_forces)  # (n_gauss, 6)
            n_gauss = len(elem_obj.gauss_data)
            elem_mean = section_g.mean(axis=0)  # (6,)

            if self._element_has_shape_functions(elem_obj, n_nodes_elem):
                N_mat = self._build_nodal_shape_matrix(elem_obj, n_nodes_elem)
                nodal_elem = np.zeros((n_nodes_elem, 6))
                for j in range(6):
                    if n_gauss == n_nodes_elem:
                        nodal_elem[:, j] = np.linalg.solve(N_mat, section_g[:, j])
                    else:
                        nodal_elem[:, j] = np.linalg.lstsq(
                            N_mat, section_g[:, j], rcond=None
                        )[0]
                for node_idx, node_id in enumerate(node_ids):
                    if node_id >= self.n_nodes:
                        continue
                    if nodal_weight[node_id] == 0:
                        nodal_elem_mean[node_id] = elem_mean
                    nodal_section_forces[node_id] += nodal_elem[node_idx]
                    nodal_weight[node_id] += 1.0
            else:
                # Fallback: distribute element mean to nodes
                for node_id in node_ids:
                    if node_id >= self.n_nodes:
                        continue
                    if nodal_weight[node_id] == 0:
                        nodal_elem_mean[node_id] = elem_mean
                    nodal_section_forces[node_id] += elem_mean
                    nodal_weight[node_id] += 1.0

        # Normalize: divide by weight where weight > 0
        nonzero = nodal_weight > 0
        nodal_section_forces[nonzero] /= nodal_weight[nonzero, np.newaxis]

        # Boundary rule: nodes with exactly one contributing element use element mean
        boundary = nodal_weight == 1.0
        valid_mean = np.isfinite(nodal_elem_mean).any(axis=1)
        use_mean = boundary & valid_mean
        nodal_section_forces[use_mean] = nodal_elem_mean[use_mean]

        return nodal_section_forces
