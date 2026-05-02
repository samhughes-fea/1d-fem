"""Build formulation cache at converged displacements for nonlinear post-processing."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, List

import numpy as np

from processing.static.results.containers.formulation_results import FormulationResultSet


def build_converged_formulation_cache(
    *,
    elements: List[Any],
    reference_cache: FormulationResultSet,
    U_global: np.ndarray,
) -> FormulationResultSet:
    """
    Update cached element stiffness matrices to the converged tangent where available.

    External load contributions in ``force_objects`` are unchanged (same as linear assembly).
    """
    U = np.asarray(U_global, dtype=np.float64).reshape(-1)
    new_element_objects = []
    for i, elem in enumerate(elements):
        ref_eo = reference_cache.element_objects[i]
        dof_idx = elem.assemble_global_dof_indices()
        U_e = U[dof_idx].astype(np.float64)
        if hasattr(elem, "tangent_stiffness_matrix"):
            K_e = np.asarray(elem.tangent_stiffness_matrix(U_e), dtype=np.float64)
            new_element_objects.append(replace(ref_eo, K_e=K_e))
        else:
            new_element_objects.append(ref_eo)

    return FormulationResultSet(
        element_objects=new_element_objects,
        force_objects=list(reference_cache.force_objects),
    )
