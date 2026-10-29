"""
Mesh-wide nodal DOF layout.

Every node carries the six standard beam DOFs ``UX, UY, UZ, RX, RY, RZ`` followed by
the same number of *extra* DOFs: warping amplitudes and, later, cross-section mode
amplitudes. The number of extra DOFs is the largest any element row asks for:

* ``[warping] = 1`` (Vlasov ``chi`` or Benscoter ``psi``) asks for one;
* an optional ``element_dictionary["extra_dofs"]`` integer per row asks for that many
  (cross-section mode elements set it).

Global numbering is node-major with a fixed stride
``mesh_dof_per_node = 6 + mesh_extra_dofs``:
``global_dof = node_id * mesh_dof_per_node + local_index``.

Element vectors keep their own order: the six standard DOFs of node 1, those of node 2,
then the element's extra DOFs of node 1, then those of node 2 (see
``Element1DBase.assemble_global_dof_indices``). An element with fewer extra DOFs than the
mesh uses the first ones of each node; the rest are not touched by it.

DOF names in ``prescribed_displacement.txt``: the six standard names, ``Q1 … Qn`` for the
extra DOFs, and ``W`` / ``CHI`` as aliases of ``Q1``.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

from pre_processing.element_library.beam_warping import element_warping_stiffness_on, mesh_uses_warping_dof

STANDARD_DOFS_PER_NODE = 6
STANDARD_DOF_NAMES = ("UX", "UY", "UZ", "RX", "RY", "RZ")
WARPING_ALIASES = ("W", "CHI")


def row_extra_dofs(element_dictionary: Dict[str, Any], idx: int) -> int:
    """Extra nodal DOFs element row ``idx`` asks for (0 for a plain six-DOF element)."""
    n = 0
    extra = element_dictionary.get("extra_dofs")
    if extra is not None:
        n = max(n, int(np.asarray(extra).ravel()[idx]))
    types = element_dictionary.get("types")
    type_str = "" if types is None else str(np.asarray(types).ravel()[idx])
    if element_warping_stiffness_on(element_dictionary, idx, type_str):
        n = max(n, 1)
    return n


def mesh_extra_dofs(element_dictionary: Dict[str, Any]) -> int:
    """Extra DOFs every node carries: the largest any row asks for."""
    n = 0
    extra = element_dictionary.get("extra_dofs")
    if extra is not None and np.asarray(extra).size:
        n = int(np.max(np.asarray(extra)))
    if mesh_uses_warping_dof(element_dictionary):
        n = max(n, 1)
    return n


def mesh_dof_per_node(element_dictionary: Dict[str, Any] | None) -> int:
    """Nodal stride of the global numbering (6 for a plain mesh, 7 with warping, …)."""
    if element_dictionary is None:
        return STANDARD_DOFS_PER_NODE
    return STANDARD_DOFS_PER_NODE + mesh_extra_dofs(element_dictionary)


def dof_names(dof_per_node: int) -> List[str]:
    """Canonical nodal DOF names for a mesh stride (``W`` for a single extra DOF)."""
    _check(dof_per_node)
    n_extra = dof_per_node - STANDARD_DOFS_PER_NODE
    if n_extra == 1:
        return list(STANDARD_DOF_NAMES) + ["W"]
    return list(STANDARD_DOF_NAMES) + [f"Q{k}" for k in range(1, n_extra + 1)]


def dof_name_map(dof_per_node: int) -> Dict[str, int]:
    """Every accepted DOF name (upper case) to its local index for a mesh stride."""
    _check(dof_per_node)
    names = {n: i for i, n in enumerate(STANDARD_DOF_NAMES)}
    n_extra = dof_per_node - STANDARD_DOFS_PER_NODE
    for k in range(1, n_extra + 1):
        names[f"Q{k}"] = STANDARD_DOFS_PER_NODE + k - 1
    if n_extra >= 1:
        for alias in WARPING_ALIASES:
            names[alias] = STANDARD_DOFS_PER_NODE
    return names


def _check(dof_per_node: int) -> None:
    if int(dof_per_node) != dof_per_node or dof_per_node < STANDARD_DOFS_PER_NODE:
        raise ValueError(f"dof_per_node must be an integer >= {STANDARD_DOFS_PER_NODE}, got {dof_per_node}")
