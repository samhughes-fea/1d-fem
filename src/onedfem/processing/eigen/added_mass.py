# processing/eigen/added_mass.py
"""Assemble the labelled strip-heave added-mass matrix onto global DOFs.

The structural mass is unchanged unless ``added_mass_from_settings`` returns a
model. ``M_total = M + M_a``. This is not a measured frequency in water.

Public API
----------
assemble_added_mass, add_labelled_added_mass, write_added_mass_gauss.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from onedfem.pre_processing.element_library.added_mass import (
    AddedMass,
    AddedMassObject,
    added_mass_from_settings,
    element_added_mass,
)

__all__ = [
    "assemble_added_mass",
    "add_labelled_added_mass",
    "write_added_mass_gauss",
    "mesh_arrays",
]


def mesh_arrays(settings: Dict[str, Any]):
    """Node coordinates and element connectivity from a runner ``settings`` dict."""
    grid = settings.get("grid_dictionary") or {}
    mesh = settings.get("mesh_dictionary") or {}
    coords = grid.get("coordinates")
    if coords is None:
        coords = mesh.get("coordinates")
    conn = (settings.get("element_dictionary") or {}).get("connectivity")
    return coords, conn


def _dof_maps(elements) -> List[np.ndarray]:
    return [np.asarray(e.assemble_global_dof_indices(), dtype=np.int64).ravel() for e in elements]


def assemble_added_mass(
    elements: Sequence,
    coordinates: npt.ArrayLike,
    connectivity: Sequence,
    model: AddedMass,
    total_dof: int,
) -> Tuple[sp.csr_matrix, List[AddedMassObject]]:
    """Global ``M_a`` (sparse) and the per-element added-mass objects."""
    xyz = np.asarray(coordinates, dtype=np.float64)
    rows, cols, data = [], [], []
    objs: List[AddedMassObject] = []
    for elem, dofs, (n0, n1) in zip(elements, _dof_maps(elements), connectivity):
        o = element_added_mass(elem, xyz[int(n0)], xyz[int(n1)], model, n_dof=dofs.size)
        objs.append(o)
        r, c = np.meshgrid(dofs, dofs, indexing="ij")
        rows.append(r.ravel())
        cols.append(c.ravel())
        data.append(o.M_a.ravel())
    if not data:
        return sp.csr_matrix((total_dof, total_dof)), objs
    M_a = sp.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(total_dof, total_dof),
    ).tocsr()
    return M_a, objs


def write_added_mass_gauss(objs: Sequence[AddedMassObject], out_dir: str) -> str:
    """Write the labelled model and ``added_mass_gauss.csv`` under ``out_dir``."""
    os.makedirs(out_dir, exist_ok=True)
    if objs:
        model = objs[0].model
        with open(os.path.join(out_dir, "added_mass_model.txt"), "w", encoding="utf-8") as handle:
            handle.write(f"model = {model.model}\n")
            handle.write(f"formula = {model.formula}\n")
            handle.write(f"source = {model.source}\n")
            handle.write(f"rho_fluid = {model.rho_fluid}\n")
            handle.write(f"chord = {model.chord}\n")
            handle.write(f"heave_local_axis = {model.heave_local_axis}\n")
            handle.write(
                "note = labelled potential-flow model; not a measured frequency in water\n"
            )
    path = os.path.join(out_dir, "added_mass_gauss.csv")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("element_id,xi,weight,jacobian,mass_per_length,e_x,e_y,e_z\n")
        for obj in objs:
            for g in obj.gauss_data:
                handle.write(
                    ",".join(
                        str(v)
                        for v in (
                            obj.element_id,
                            g.xi,
                            g.weight,
                            g.jacobian,
                            g.mass_per_length,
                            *g.direction,
                        )
                    )
                    + "\n"
                )
    return path


def add_labelled_added_mass(
    M: sp.spmatrix,
    *,
    elements: Sequence,
    coordinates: npt.ArrayLike,
    connectivity: Sequence,
    simulation_settings: Dict[str, Any],
    total_dof: int,
    results_dir: Optional[str] = None,
) -> Tuple[sp.csr_matrix, Optional[List[AddedMassObject]]]:
    """Add strip-heave ``M_a`` when the labelled model is on; otherwise return ``M``.

    When on and ``results_dir`` is set, writes ``added_mass/added_mass_gauss.csv``
    under that directory.
    """
    model = added_mass_from_settings(simulation_settings)
    M_csr = sp.csr_matrix(M)
    if model is None:
        return M_csr, None
    M_a, objs = assemble_added_mass(elements, coordinates, connectivity, model, total_dof)
    if results_dir:
        write_added_mass_gauss(objs, os.path.join(results_dir, "added_mass"))
    return (M_csr + M_a).tocsr(), objs
