# pre_processing/mesh_library/blade_beam_job.py
"""A 1d-fem beam job for a blade or wing lofted by implicit-geometry-kernel.

``blade_beam_job(job_dir, geometry, laminates, n_elements)`` writes the whole job for a
``BladeGeometry``: the reference axis, the elements and their twist rate, one section
per element cut at its mid-span, the root clamp and any nodal loads.

Axes
----
The geometry's blade frame has ``z`` along the span, ``x`` along the chord (towards the
trailing edge) and ``y`` through the thickness (suction side positive). The job's
global axes are a cyclic relabelling, which keeps them right handed::

    X = z (span),   Y = x (chord),   Z = y (flap)

so a straight blade lies along global ``X`` and its elements' own axes (``x`` along
the beam, ``y`` = ``Z cross x``, ``z`` up; ``element_frame.straight_frame``) are the
section axes ``(y_s, z_s) = (x, y)`` the section builders use. ``blade_to_job`` maps
blade-frame vectors to job axes.

Reference axis
--------------
Nodes sit on the pitch axis: node ``i`` at ``(z_i, offset_x(z_i), offset_y(z_i))``,
the station offsets (sweep, prebend) interpolated as the geometry interpolates them.
A prebent or swept blade is a chain of straight elements along that polyline. The
sections are measured from the pitch-axis point of their cut, so they ride on the
axis.

Sections and twist
------------------
Each element takes the section cut at its mid-span, ``geometry.section_at(z_mid)``,
as the geometry builds it: already twisted, so a layout's flapwise webs (``"-F"``)
stay normal to the rotor plane. That is why the ``[twist]`` column of ``element.txt``
is not used: a twisted ``-F`` cut is not a turned copy of the untwisted one. The twist
gradient becomes the element's pretwist rate, ``k_x0`` in ``precurvature.txt``. The
geometry's twist is nose up, a turn of the section by ``-t`` about the span, so

    k_x0 = -(t(z_1) - t(z_0)) / L          (radians per metre)

over an element of length ``L`` between nodes at ``z_0`` and ``z_1``. A pretwisted
element is the twisting arc of ``framed_element``: its section turns at ``k_x0`` along
it, matching the cut at its middle, and the helical fibres add their extension-twist
coupling (``SectionAnalysis.with_pretwist``).

Cuts lie in planes of constant ``z``, not normal to a prebent or swept axis; the
difference is second order in the axis slope.

``kind="thin_walled"`` writes ``laminate.txt`` and ``section_geometry.txt``
(``thin_walled_from_cut``); ``kind="2d"`` writes ``section_mesh_2d.txt``
(``section_mesh_from_cut``, a conforming layered mesh). Both write
``element_section.txt``, section ``e + 1`` for element ``e``.

Loads are nodal point loads (``loads``) and line loads per length of the axis
(``line_loads``, written to ``distributed_load.txt`` against span position). The
classical ``material.txt`` and ``section.txt`` rows the parsers need are placeholders:
the sections replace them.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Sequence, Union

import numpy as np

__all__ = ["BladeBeamJob", "blade_beam_job", "blade_to_job", "job_to_blade"]

STANDARD_DOFS = ("UX", "UY", "UZ", "RX", "RY", "RZ")
FRAMED_TYPES = ("LinearEulerBernoulliBeamElement3D", "LinearTimoshenkoBeamElement3D", "LinearBenscoterBeamElement3D")
_SETTINGS = "[Simulation]\n[Type]\nstatic\n\n[Solver]\ntype = direct\ntolerance = 1e-12\nmax_iterations = 500\n"


def blade_to_job(v: Sequence[float]) -> np.ndarray:
    """Blade-frame ``(x, y, z)`` components (or a force-moment 6-vector) in job axes ``(X, Y, Z)``."""
    v = np.asarray(v, dtype=np.float64)
    if v.shape == (3,):
        return v[[2, 0, 1]]
    if v.shape == (6,):
        return np.concatenate([v[[2, 0, 1]], v[[5, 3, 4]]])
    raise ValueError("expected 3 or 6 components")


def job_to_blade(v: Sequence[float]) -> np.ndarray:
    """Job-axis ``(X, Y, Z)`` components (or a 6-vector) in the blade frame ``(x, y, z)``."""
    v = np.asarray(v, dtype=np.float64)
    if v.shape == (3,):
        return v[[1, 2, 0]]
    if v.shape == (6,):
        return np.concatenate([v[[1, 2, 0]], v[[4, 5, 3]]])
    raise ValueError("expected 3 or 6 components")


@dataclass
class BladeBeamJob:
    """What ``blade_beam_job`` wrote.

    Attributes
    ----------
    job_dir : str
    kind : str
        ``"thin_walled"`` or ``"2d"``.
    z : ndarray (n + 1,)
        Span position of each node.
    nodes : ndarray (n + 1, 3)
        Node coordinates in job axes.
    z_cut : ndarray (n,)
        Span position of each element's section cut.
    twist : ndarray (n + 1,)
        Geometry twist at each node, degrees (nose up).
    pretwist : ndarray (n,)
        ``k_x0`` of each element, rad/m.
    sections : dict
        Section id to its ``GeometrySection`` (thin-walled) or ``SectionMesh`` (2D).
    element_type : str
    """

    job_dir: str
    kind: str
    z: np.ndarray
    nodes: np.ndarray
    z_cut: np.ndarray
    twist: np.ndarray
    pretwist: np.ndarray
    sections: Dict[int, object]
    element_type: str
    loads: Dict[int, np.ndarray] = field(default_factory=dict)

    @property
    def n_elements(self) -> int:
        return int(self.z_cut.size)

    def node_at(self, z: float) -> int:
        """Index of the node at span position ``z``."""
        i = int(np.argmin(np.abs(self.z - z)))
        if abs(self.z[i] - z) > 1e-9 * max(1.0, float(np.ptp(self.z))):
            raise ValueError(f"no node at z = {z}; nodes are at {self.z.tolist()}")
        return i


def _span_nodes(geometry, n_elements: Optional[int], z_nodes: Optional[Sequence[float]]) -> np.ndarray:
    z0, z1 = geometry.stations[0].z, geometry.stations[-1].z
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


def _laminate_key(lam):
    return (tuple((p.material, float(p.angle_deg), float(p.thickness)) for p in lam.plies), float(lam.shear_correction))


def blade_beam_job(
    job_dir: str,
    geometry,
    laminates: Mapping[str, object],
    n_elements: Optional[int] = None,
    kind: str = "thin_walled",
    *,
    z_nodes: Optional[Sequence[float]] = None,
    library: Optional[Mapping[int, object]] = None,
    element_type: str = "LinearTimoshenkoBeamElement3D",
    pretwist: bool = True,
    n_sub: int = 1,
    mesh_options: Optional[Mapping[str, object]] = None,
    loads: Optional[Mapping[float, Sequence[float]]] = None,
    line_loads: Optional[Mapping[float, Sequence[float]]] = None,
    thickness_tol: float = 0.01,
    settings: str = _SETTINGS,
) -> BladeBeamJob:
    """Write a cantilever beam job for ``geometry`` into ``job_dir`` (replaced).

    Parameters
    ----------
    job_dir : str
    geometry : implicit_geometry_kernel.BladeGeometry
    laminates : mapping
        Wall label, region or kind to a ``WallLaminate`` (or a laminate id in
        ``library``), as ``thin_walled_from_cut`` takes it.
    n_elements : int, optional
        Equal elements from the first station to the last. Or give ``z_nodes``.
    kind : {"thin_walled", "2d"}
    z_nodes : sequence of float, optional
        Span positions of the nodes.
    library : mapping int -> WallLaminate, optional
    element_type : str
        A beam type that takes rotated and curved elements: ``LinearTimoshenkoBeamElement3D``
        (default), ``LinearEulerBernoulliBeamElement3D`` or ``LinearBenscoterBeamElement3D``
        (thin-walled sections; the root warping is held).
    pretwist : bool
        Write the twist gradient as each element's pretwist rate (default). Off, each
        element keeps its mid-span section along its length.
    n_sub : int
        Thin-walled wall subdivision (``[Section_Mesh]``).
    mesh_options : mapping, optional
        Keyword arguments of ``section_mesh_from_cut`` for ``kind="2d"`` (``element_size``,
        ``rows_per_ply``).
    loads : mapping z -> 6 floats, optional
        Nodal loads ``[F_X, F_Y, F_Z, M_X, M_Y, M_Z]`` in job axes at nodes ``z``
        (``blade_to_job`` turns blade-frame loads).
    line_loads : mapping z -> 6 floats, optional
        Distributed loads per length of the axis, ``[q_X, q_Y, q_Z, m_X, m_Y, m_Z]`` in job
        axes, sampled at span positions ``z`` (linear between them; the samples must cover
        the blade).
    thickness_tol : float
    settings : str
        ``simulation_settings.txt`` text (default linear static, direct solver).

    Returns
    -------
    BladeBeamJob
    """
    if kind not in ("thin_walled", "2d"):
        raise ValueError(f"kind must be 'thin_walled' or '2d', got {kind!r}")
    if element_type not in FRAMED_TYPES:
        raise ValueError(f"element_type must be one of {', '.join(FRAMED_TYPES)} (types that take a blade's "
                         f"rotated and pretwisted elements); got {element_type!r}")
    benscoter = "Benscoter" in element_type
    if benscoter and kind != "thin_walled":
        raise ValueError("LinearBenscoterBeamElement3D needs thin-walled sections")
    z = _span_nodes(geometry, n_elements, z_nodes)
    n = z.size - 1
    stations = [geometry.station_at(float(zi)) for zi in z]
    twist = np.array([s[1] for s in stations])
    nodes = np.array([[zi, s[3][0], s[3][1]] for zi, s in zip(z, stations)], dtype=np.float64)
    L = np.linalg.norm(np.diff(nodes, axis=0), axis=1)
    k_x0 = -np.radians(np.diff(twist)) / L if pretwist else np.zeros(n)
    z_cut = 0.5 * (z[:-1] + z[1:])
    cuts = [geometry.section_at(float(zc)) for zc in z_cut]

    if os.path.isdir(job_dir):
        shutil.rmtree(job_dir)
    os.makedirs(job_dir)

    def write(name, text):
        with open(os.path.join(job_dir, name), "w", encoding="utf-8") as fh:
            fh.write(text)

    sections: Dict[int, object] = {}
    if kind == "thin_walled":
        from pre_processing.parsing.laminate_parser import format_laminate_file
        from pre_processing.parsing.section_geometry_parser import format_section_geometry
        from pre_processing.section_integration.thin_walled.from_geometry import thin_walled_from_cut
        from pre_processing.section_integration.thin_walled.thin_walled_geometry import ThinWalledSection, Wall

        lam_ids: Dict[object, int] = {}
        lam_by_id: Dict[int, object] = {}
        renumbered: Dict[int, ThinWalledSection] = {}
        for e, cut in enumerate(cuts):
            sid = e + 1
            g = thin_walled_from_cut(cut, laminates, library, thickness_tol=thickness_tol, section_id=sid)
            local = {}
            for lid, lam in g.laminates.items():
                key = _laminate_key(lam)
                if key not in lam_ids:
                    lam_ids[key] = len(lam_ids) + 1
                    lam_by_id[lam_ids[key]] = lam
                local[lid] = lam_ids[key]
            sec = g.section
            walls = [Wall(w.wall_id, w.start, w.end, local[w.laminate_id]) for w in sec.walls]
            renumbered[sid] = ThinWalledSection(dict(sec.points), walls, sorted(set(local.values())), section_id=sid)
            sections[sid] = g
        write("laminate.txt", format_laminate_file(lam_by_id))
        write("section_geometry.txt", format_section_geometry(renumbered, {sid: n_sub for sid in renumbered}))
    else:
        from pre_processing.section_integration.solid.from_geometry import section_mesh_from_cut, write_section_mesh_2d
        from pre_processing.section_integration.solid.material import SolidMaterial

        shared: Dict[object, SolidMaterial] = {}
        for value in list(laminates.values()) + list((library or {}).values()):
            for ply in getattr(value, "plies", ()):
                shared.setdefault(ply.material, SolidMaterial.from_ply(ply.material))
        opts = dict(mesh_options or {})
        for e, cut in enumerate(cuts):
            sections[e + 1] = section_mesh_from_cut(cut, laminates, library, thickness_tol=thickness_tol,
                                                    materials=shared, section_id=e + 1, **opts)
        write_section_mesh_2d(os.path.join(job_dir, "section_mesh_2d.txt"), sections)

    write("grid.txt", "[Grid]\n[node_id] [x] [y] [z]\n"
          + "".join(f"{i} {x!r} {y!r} {zz!r}\n" for i, (x, y, zz) in enumerate(nodes.tolist())))
    warp_col, warp = (" [warping]", " 1") if benscoter else ("", "")
    write("element.txt",
          "[Element]\n[element_id] [node1] [node2] [element_type] [axial_order] [bending_y_order] "
          f"[bending_z_order] [shear_y_order] [shear_z_order] [torsion_order] [load_order]{warp_col}\n"
          + "".join(f"{e} {e} {e + 1} {element_type} 2 2 2 2 2 2 2{warp}\n" for e in range(n)))
    if np.any(k_x0 != 0.0):
        write("precurvature.txt", "[Precurvature]\n[element_id] [k_x0] [k_y0] [k_z0]\n"
              + "".join(f"{e} {float(k)!r} 0.0 0.0\n" for e, k in enumerate(k_x0)))
    write("material.txt", "[Material]\n[element_id] [E] [G] [nu] [rho]\n" + "".join(f"{e} 7e10 2.7e10 0.3 0\n" for e in range(n)))
    write("section.txt", "[Section]\n[element_id] [A] [I_x] [I_y] [I_z] [J_t] [kappa] [alpha]\n"
          + "".join(f"{e} 1e-3 0.0 1e-8 1e-8 1e-9 0.8333 0.0\n" for e in range(n)))
    write("element_section.txt", "[Element_Section]\n[element_id] [section_id]\n" + "".join(f"{e} {e + 1}\n" for e in range(n)))
    dofs = list(STANDARD_DOFS) + (["Q1"] if benscoter else [])
    write("prescribed_displacement.txt", "[Prescribed Displacement]\n[id] [node_id] [dof] [value] [type] [comment]\n"
          + "".join(f"{i} 0 {d} 0.0 displacement\n" for i, d in enumerate(dofs)))
    job = BladeBeamJob(job_dir, kind, z, nodes, z_cut, twist, k_x0, sections, element_type)
    rows = []
    for zl, f in (loads or {}).items():
        i = job.node_at(float(zl))
        f = np.asarray(f, dtype=np.float64).reshape(6)
        job.loads[i] = job.loads.get(i, np.zeros(6)) + f
    for i, f in sorted(job.loads.items()):
        rows.append(" ".join(repr(float(v)) for v in nodes[i]) + " " + " ".join(repr(float(v)) for v in f) + "\n")
    if line_loads:
        zs = sorted(float(k) for k in line_loads)
        if zs[0] > z[0] + 1e-9 * max(1.0, float(np.ptp(z))) or zs[-1] < z[-1] - 1e-9 * max(1.0, float(np.ptp(z))):
            raise ValueError(f"line_loads must cover the blade, z = {z[0]} to {z[-1]}; they cover {zs[0]} to {zs[-1]}")
        write("distributed_load.txt", "[Distributed load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n"
              + "".join(f"{zl!r} 0.0 0.0 " + " ".join(repr(float(v)) for v in np.asarray(line_loads[zl], float).reshape(6))
                        + "\n" for zl in sorted(line_loads, key=float)))
    write("point_load.txt", "[Point load]\n[x] [y] [z] [F_x] [F_y] [F_z] [M_x] [M_y] [M_z]\n" + "".join(rows))
    write("simulation_settings.txt", settings)
    return job
