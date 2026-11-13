"""A beam job from kernel stations and a reduced section.

The stations are span, twist and offset. Each element gets a reduced beam
section (area, centroidal moments, shear centre, warping constant) already
computed by the section model. This writer does not ask the geometry for a
cut, a wall or a shell graph.

Section axes follow the blade job: beam ``x`` is span, section ``y`` is chord
and section ``z`` is thickness. ``ixx`` (``∫ y² dA`` in the blade frame) is
written as ``I_y``; ``iyy`` (``∫ x² dA``) is written as ``I_z``. The shear
centre ``(x, y)`` in the blade frame is ``(y_sc, z_sc)``.
"""

from __future__ import annotations

import os
import shutil
from collections.abc import Callable, Mapping, Sequence
from typing import Optional, Union

import numpy as np

from onedfem.pre_processing.mesh_library.blade_beam_job import (
    NONLINEAR_SETTINGS,
    STANDARD_DOFS,
)

_SETTINGS = "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n"
_TYPES = (
    "LinearEulerBernoulliBeamElement3D",
    "LinearTimoshenkoBeamElement3D",
    "CorotationalEulerBernoulliBeamElement3D",
    "CorotationalTimoshenkoBeamElement3D",
)
_KAPPA = 5.0 / 6.0


def reduced_blade_beam_job(
    job_dir: str,
    stations,
    sections: Union[Sequence[object], Callable],
    n_elements: Optional[int] = None,
    *,
    z_nodes: Optional[Sequence[float]] = None,
    element_type: str = "LinearTimoshenkoBeamElement3D",
    pretwist: bool = True,
    young: float = 7.0e10,
    shear_modulus: float = 2.7e10,
    nu: float = 0.3,
    rho: float = 0.0,
    loads: Optional[Mapping[float, Sequence[float]]] = None,
    settings: Optional[str] = None,
):
    """Write a cantilever whose sections are reduced beam properties.

    ``stations`` is a sequence, or an object with ``.stations``, of spanwise
    stations (``z``, ``twist`` in degrees, ``offset`` ``(x, y)``). ``sections``
    is one reduced section per element, or a callable of the element mid-span
    ``z`` returning one. Twist and offset are interpolated on the station
    table.
    """
    if element_type not in _TYPES:
        raise ValueError(f"element_type must be one of {', '.join(_TYPES)}; got {element_type!r}")
    table = _table(stations)
    z = _nodes(table, n_elements, z_nodes)
    n = z.size - 1
    twist = np.array([_interp(table, "twist", float(zi)) for zi in z])
    offset = np.array([_offset(table, float(zi)) for zi in z])
    nodes = np.column_stack([z, offset])
    length = np.linalg.norm(np.diff(nodes, axis=0), axis=1)
    k_x0 = -np.radians(np.diff(twist)) / length if pretwist else np.zeros(n)
    z_cut = 0.5 * (z[:-1] + z[1:])
    props = _sections(sections, z_cut)

    if os.path.isdir(job_dir):
        shutil.rmtree(job_dir)
    os.makedirs(job_dir)

    def write(name, text):
        with open(os.path.join(job_dir, name), "w", encoding="utf-8") as fh:
            fh.write(text)

    write("grid.txt", "[Grid]\n[node_id] [x] [y] [z]\n"
          + "".join(f"{i} {x!r} {y!r} {zz!r}\n" for i, (x, y, zz) in enumerate(nodes.tolist())))
    write("element.txt",
          "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
          "[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order]\n"
          + "".join(f"{e} {e} {e + 1} {element_type} 2 2 2 2 2 2 2\n" for e in range(n)))
    if np.any(k_x0 != 0.0):
        write("precurvature.txt", "[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n"
              + "".join(f"{e} {float(k)!r} 0.0 0.0\n" for e, k in enumerate(k_x0)))
    write("material.txt", "[Material]\n[element_id] [E] [G] [nu] [rho]\n"
          + "".join(f"{e} {young!r} {shear_modulus!r} {nu!r} {rho!r}\n" for e in range(n)))
    write("section.txt",
          "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha] [y_sc] [z_sc] [Gamma]\n"
          + "".join(_section_row(e, p) for e, p in enumerate(props)))
    write("prescribed_displacement.txt",
          "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
          + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(STANDARD_DOFS)))
    rows = []
    for zl, force in (loads or {}).items():
        i = int(np.argmin(np.abs(z - float(zl))))
        if abs(z[i] - float(zl)) > 1e-9 * max(1.0, float(np.ptp(z))):
            raise ValueError(f"no node at z = {zl}")
        rows.append(" ".join(repr(float(v)) for v in nodes[i]) + " "
                    + " ".join(repr(float(v)) for v in np.asarray(force, float).reshape(6)) + "\n")
    write("point_load.txt", "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n" + "".join(rows))
    if settings is None:
        settings = NONLINEAR_SETTINGS if element_type.startswith("Corotational") else _SETTINGS
    write("simulation_settings.txt", settings)
    return {"z": z, "z_cut": z_cut, "nodes": nodes, "twist": twist, "k_x0": k_x0}


def _table(stations):
    rows = list(stations.stations if hasattr(stations, "stations") else stations)
    if len(rows) < 2:
        raise ValueError("need at least two stations")
    rows.sort(key=lambda s: float(s.z))
    if any(float(b.z) <= float(a.z) for a, b in zip(rows, rows[1:])):
        raise ValueError("station z must increase")
    return rows


def _nodes(table, n_elements, z_nodes) -> np.ndarray:
    z0, z1 = float(table[0].z), float(table[-1].z)
    if z_nodes is not None:
        z = np.asarray(z_nodes, dtype=np.float64)
        if z.ndim != 1 or z.size < 2 or np.any(np.diff(z) <= 0.0):
            raise ValueError("z_nodes must be at least two increasing span positions")
        if z[0] < z0 - 1e-9 or z[-1] > z1 + 1e-9:
            raise ValueError(f"z_nodes must lie between the first and last stations ({z0} to {z1})")
        return z
    if n_elements is None or int(n_elements) < 1:
        raise ValueError("give n_elements (>= 1) or z_nodes")
    return np.linspace(z0, z1, int(n_elements) + 1)


def _attr(station, name, default):
    if isinstance(station, Mapping):
        return station.get(name, default)
    return getattr(station, name, default)


def _interp(table, name, zq: float) -> float:
    zs = np.array([float(s.z) for s in table])
    values = np.array([float(_attr(s, name, 0.0)) for s in table])
    return float(np.interp(zq, zs, values))


def _offset(table, zq: float):
    zs = np.array([float(s.z) for s in table])
    ox, oy = [], []
    for station in table:
        off = _attr(station, "offset", (0.0, 0.0))
        ox.append(float(off[0]))
        oy.append(float(off[1]))
    return float(np.interp(zq, zs, ox)), float(np.interp(zq, zs, oy))


def _field(section, name, default=None):
    if isinstance(section, Mapping):
        if name in section and section[name] is not None:
            return section[name]
    elif hasattr(section, name):
        value = getattr(section, name)
        if value is not None:
            return value
    if default is not None:
        return default
    raise TypeError(f"reduced section needs {name}")


def _sections(sections, z_cut: np.ndarray):
    if callable(sections):
        return [sections(float(z)) for z in z_cut]
    props = list(sections)
    if len(props) != z_cut.size:
        raise ValueError(f"expected {z_cut.size} reduced sections, got {len(props)}")
    return props


def _section_row(eid: int, section) -> str:
    centre = _field(section, "shear_centre", _field(section, "centroid", (0.0, 0.0)))
    y_sc, z_sc = float(centre[0]), float(centre[1])
    kappa = float(_field(section, "kappa", _KAPPA))
    return (
        f"{eid} {_field(section, 'area')!r} 0.0 {_field(section, 'ixx')!r} {_field(section, 'iyy')!r} "
        f"{float(_field(section, 'torsion_constant', 0.0))!r} {kappa!r} 0.0 {y_sc!r} {z_sc!r} "
        f"{float(_field(section, 'warping_constant', 0.0))!r}\n"
    )
