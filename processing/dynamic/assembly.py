# processing/dynamic/assembly.py
"""
Assemble global K, M, and optionally C for dynamic analysis.
Copied and extended from modal assembly; no imports from processing.static.
"""

import logging
import numpy as np
from pathlib import Path
from scipy.sparse import coo_matrix, csr_matrix
from typing import List, Optional, Tuple, Union

_logger = logging.getLogger(__name__)


def _compute_local_global_dof_map(elements: List, total_dof: int) -> List[np.ndarray]:
    """Compute global DOF indices for each element."""
    local_global_dof_map = []
    for elem in elements:
        element_id = int(elem.element_id)
        dof = elem.assemble_global_dof_indices()
        validated_dof = np.asarray(dof, dtype=np.int32).ravel()
        if validated_dof.size == 0:
            raise ValueError(f"Element {element_id}: empty DOF mapping")
        if validated_dof.min() < 0:
            raise ValueError(f"Element {element_id}: negative DOF index")
        if validated_dof.max() >= total_dof:
            raise ValueError(f"Element {element_id}: DOF index >= total_dof {total_dof}")
        if len(np.unique(validated_dof)) != len(validated_dof):
            raise ValueError(f"Element {element_id}: duplicate DOF indices")
        local_global_dof_map.append(validated_dof)
    return local_global_dof_map


def _scatter_element_matrix(
    mat: Union[np.ndarray, coo_matrix], dof_map: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (rows, cols, data) in global indexing for one element matrix."""
    if dof_map.size == 0:
        return np.array([]), np.array([]), np.array([])
    if isinstance(mat, np.ndarray):
        coo = coo_matrix(mat)
    else:
        coo = mat.tocoo() if not isinstance(mat, coo_matrix) else mat
    if coo.shape != (dof_map.size, dof_map.size):
        raise ValueError(f"Matrix shape {coo.shape} does not match DOF size {dof_map.size}")
    rows = dof_map[coo.row]
    cols = dof_map[coo.col]
    return rows.astype(np.int32), cols.astype(np.int32), coo.data.astype(np.float64)


def assemble_global_system(
    elements: List,
    element_stiffness_matrices: Optional[List] = None,
    element_mass_matrices: Optional[List] = None,
    element_damping_matrices: Optional[List] = None,
    total_dof: Optional[int] = None,
    job_results_dir: Optional[str] = None,
) -> Tuple[csr_matrix, csr_matrix, Optional[csr_matrix], List[np.ndarray]]:
    """
    Assemble global K, M, and optionally C for dynamic analysis.

    Returns
    -------
    K_global : csr_matrix
    M_global : csr_matrix
    C_global : csr_matrix or None
    local_global_dof_map : list of ndarray
    """
    if not elements:
        raise ValueError("No elements provided")
    if total_dof is None or total_dof <= 0:
        raise ValueError("total_dof must be a positive integer")
    n_el = len(elements)
    if element_stiffness_matrices is not None and len(element_stiffness_matrices) != n_el:
        raise ValueError("element_stiffness_matrices length does not match elements")
    if element_mass_matrices is not None and len(element_mass_matrices) != n_el:
        raise ValueError("element_mass_matrices length does not match elements")
    if element_damping_matrices is not None and len(element_damping_matrices) != n_el:
        raise ValueError("element_damping_matrices length does not match elements")

    job_path = Path(job_results_dir) if job_results_dir else None
    if job_path:
        logs_dir = job_path.parent / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)

    local_global_dof_map = _compute_local_global_dof_map(elements, total_dof)

    all_k_rows, all_k_cols, all_k_data = [], [], []
    all_m_rows, all_m_cols, all_m_data = [], [], []
    all_c_rows, all_c_cols, all_c_data = [], [], []

    if element_stiffness_matrices is not None:
        for Ke, dof in zip(element_stiffness_matrices, local_global_dof_map):
            r, c, d = _scatter_element_matrix(Ke, dof)
            all_k_rows.append(r)
            all_k_cols.append(c)
            all_k_data.append(d)
    if element_mass_matrices is not None:
        for Me, dof in zip(element_mass_matrices, local_global_dof_map):
            r, c, d = _scatter_element_matrix(Me, dof)
            all_m_rows.append(r)
            all_m_cols.append(c)
            all_m_data.append(d)
    if element_damping_matrices is not None:
        for Ce, dof in zip(element_damping_matrices, local_global_dof_map):
            r, c, d = _scatter_element_matrix(Ce, dof)
            all_c_rows.append(r)
            all_c_cols.append(c)
            all_c_data.append(d)

    if all_k_data:
        K_global = coo_matrix(
            (
                np.concatenate(all_k_data),
                (np.concatenate(all_k_rows), np.concatenate(all_k_cols)),
            ),
            shape=(total_dof, total_dof),
            dtype=np.float64,
        ).tocsr()
    else:
        K_global = csr_matrix((total_dof, total_dof), dtype=np.float64)

    if all_m_data:
        M_global = coo_matrix(
            (
                np.concatenate(all_m_data),
                (np.concatenate(all_m_rows), np.concatenate(all_m_cols)),
            ),
            shape=(total_dof, total_dof),
            dtype=np.float64,
        ).tocsr()
    else:
        M_global = csr_matrix((total_dof, total_dof), dtype=np.float64)

    if all_c_data:
        C_global = coo_matrix(
            (
                np.concatenate(all_c_data),
                (np.concatenate(all_c_rows), np.concatenate(all_c_cols)),
            ),
            shape=(total_dof, total_dof),
            dtype=np.float64,
        ).tocsr()
    else:
        C_global = None

    _logger.info(
        "Dynamic assembly: K nnz=%s, M nnz=%s, C nnz=%s",
        K_global.nnz,
        M_global.nnz,
        C_global.nnz if C_global is not None else 0,
    )
    return K_global, M_global, C_global, local_global_dof_map
