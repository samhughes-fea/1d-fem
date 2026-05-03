# processing/boundary_supports/__init__.py
"""
Shared penalty-BC fixed-DOF resolution for spectral (eigen/buckling) and transient runners.

Grid files in this project do not carry a separate ``supports`` array; supports are expressed
through ``prescribed_displacement_dict`` (zeros = fixed) and optional ``fixed_node_id`` in
``[Eigen]`` / ``[Transient]`` (see ``simulation_settings_parser``).

Global DOF numbering follows ``Element1DBase.assemble_global_dof_indices``:
``global_dof = node_id * dof_per_node + local_dof_index``.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

import numpy as np


def global_dofs_for_clamped_node(node_id: int, dof_per_node: int) -> np.ndarray:
    """All global DOF indices belonging to one clamped node."""
    base = int(node_id) * int(dof_per_node)
    return np.arange(base, base + int(dof_per_node), dtype=np.int32)


def _validate_fixed_node_id(node_id: int, grid_node_ids: Optional[Sequence[Any]]) -> None:
    if grid_node_ids is None:
        return
    ids = np.asarray(grid_node_ids, dtype=np.int64).ravel()
    if int(node_id) not in set(int(x) for x in ids.tolist()):
        raise ValueError(f"fixed_node_id={node_id} is not present in mesh node_ids")


def resolve_penalty_fixed_dofs(
    *,
    total_dof: int,
    dof_per_node: int,
    prescribed_displacement_dict: Optional[Mapping[str, Any]],
    section_settings: Optional[Mapping[str, Any]],
    grid_node_ids: Optional[Sequence[Any]] = None,
) -> Optional[np.ndarray]:
    """
    Return ``fixed_dofs`` to pass into penalty BC helpers alongside ``prescribed_displacements``.

    * If ``prescribed_displacement_dict`` is **absent**: constrain either all DOFs at
      ``section_settings['fixed_node_id']`` when set, else the legacy anchor
      ``0 .. min(5, total_dof-1)`` (cantilever-style root when node numbering starts at 0).
    * If **present**: prescribed rows with zero value are handled inside
      ``apply_boundary_conditions``; this returns **only** extra DOFs from ``fixed_node_id``
      when set, otherwise ``None``.
    """
    nid = None
    if section_settings:
        raw = section_settings.get("fixed_node_id")
        if raw is not None:
            nid = int(raw)
            _validate_fixed_node_id(nid, grid_node_ids)

    if nid is not None:
        return global_dofs_for_clamped_node(nid, dof_per_node)

    if prescribed_displacement_dict is None:
        n_anchor = min(6, int(total_dof))
        return np.arange(n_anchor, dtype=np.int32) if n_anchor > 0 else None

    return None
