# pre_processing/parsing/section_geometry_parser.py
"""Parse ``section_geometry.txt``: mid-line points and walls of thin-walled sections.

File layout (comments start with ``#``)::

    [Section_Point]
    [section_id] [point_id] [y] [z]
    1  1  0.0  0.0
    1  2  0.2  0.0

    [Section_Wall]
    [section_id] [wall_id] [point_start] [point_end] [laminate_id]
    1  1  1  2  1

An optional block sets the number of wall elements per wall used by the section
analysis and the cross-section modes (default 1; more resolves local plate modes
and shear lag)::

    [Section_Mesh]
    [section_id] [n_sub]
    1  4

One file may hold several sections. See
docs/conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md.
"""

from __future__ import annotations

from typing import Dict, Iterable, Optional, Tuple

from pre_processing.parsing.laminate_parser import (
    to_float,
    to_int,
    preprocess_lines,
    split_blocks,
    tokens_match,
)
from pre_processing.section_integration.thin_walled.thin_walled_geometry import (
    ThinWalledSection,
    Wall,
)

POINT_HEADER = "[section_point]"
WALL_HEADER = "[section_wall]"
POINT_SUBHEADER = ["[section_id]", "[point_id]", "[y]", "[z]"]
WALL_SUBHEADER = ["[section_id]", "[wall_id]", "[point_start]", "[point_end]", "[laminate_id]"]
MESH_HEADER = "[section_mesh]"
MESH_SUBHEADER = ["[section_id]", "[n_sub]"]


def parse_section_mesh(filepath: str) -> Dict[int, int]:
    """Return ``{section_id: n_sub}`` from the optional ``[Section_Mesh]`` block (empty without it)."""
    blocks = split_blocks(preprocess_lines(filepath), [POINT_HEADER, WALL_HEADER, MESH_HEADER])
    if MESH_HEADER not in blocks:
        return {}
    sub, rows = blocks[MESH_HEADER]
    if not tokens_match(sub, MESH_SUBHEADER):
        raise ValueError(f"[Section_Mesh] sub-header must be: {' '.join(MESH_SUBHEADER)}")
    out: Dict[int, int] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 2:
            raise ValueError(f"Malformed [Section_Mesh] row (expected 2 columns): {ln!r}")
        sid, n = to_int(parts[0], "section_id", ln), to_int(parts[1], "n_sub", ln)
        if n < 1:
            raise ValueError(f"[Section_Mesh] n_sub must be >= 1: {ln!r}")
        out[sid] = n
    return out


def parse_section_geometry_file(
    filepath: str, laminate_ids: Optional[Iterable[int]] = None
) -> Dict[int, ThinWalledSection]:
    """Read ``section_geometry.txt`` and return ``{section_id: ThinWalledSection}``.

    Parameters
    ----------
    filepath : str
    laminate_ids : iterable of int, optional
        Known laminate ids (from ``laminate.txt``); walls naming others are rejected.
    """
    blocks = split_blocks(preprocess_lines(filepath), [POINT_HEADER, WALL_HEADER, MESH_HEADER])
    for h in (POINT_HEADER, WALL_HEADER):
        if h not in blocks:
            raise ValueError(f"Missing {h} block in {filepath}")

    sub, rows = blocks[POINT_HEADER]
    if not tokens_match(sub, POINT_SUBHEADER):
        raise ValueError(f"[Section_Point] sub-header must be: {' '.join(POINT_SUBHEADER)}")
    points: Dict[int, Dict[int, Tuple[float, float]]] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 4:
            raise ValueError(f"Malformed [Section_Point] row (expected 4 columns): {ln!r}")
        sid, pid = to_int(parts[0], "section_id", ln), to_int(parts[1], "point_id", ln)
        y, z = to_float(parts[2], "y", ln), to_float(parts[3], "z", ln)
        sec = points.setdefault(sid, {})
        if pid in sec:
            raise ValueError(f"Section {sid} has duplicate point_id {pid}")
        sec[pid] = (y, z)

    sub, rows = blocks[WALL_HEADER]
    if not tokens_match(sub, WALL_SUBHEADER):
        raise ValueError(f"[Section_Wall] sub-header must be: {' '.join(WALL_SUBHEADER)}")
    walls: Dict[int, list] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 5:
            raise ValueError(f"Malformed [Section_Wall] row (expected 5 columns): {ln!r}")
        sid, wid, a, b, lid = (to_int(t, "integer field", ln) for t in parts)
        walls.setdefault(sid, []).append(Wall(wall_id=wid, start=a, end=b, laminate_id=lid))

    missing_walls = sorted(set(points) - set(walls))
    if missing_walls:
        raise ValueError(f"Sections {missing_walls} have points but no walls")
    missing_points = sorted(set(walls) - set(points))
    if missing_points:
        raise ValueError(f"Sections {missing_points} have walls but no points")

    lam = None if laminate_ids is None else list(laminate_ids)
    return {
        sid: ThinWalledSection(points[sid], walls[sid], laminate_ids=lam, section_id=sid)
        for sid in sorted(points)
    }


def format_section_geometry(sections: Dict[int, ThinWalledSection], n_sub: Optional[Dict[int, int]] = None) -> str:
    """Text of a ``section_geometry.txt`` for ``{section_id: ThinWalledSection}``.

    ``n_sub`` (``{section_id: n_sub}``) adds a ``[Section_Mesh]`` block.
    """
    pts = ["[Section_Point]", " ".join(POINT_SUBHEADER)]
    walls = ["[Section_Wall]", " ".join(WALL_SUBHEADER)]
    for sid, sec in sections.items():
        pts += [f"{sid} {pid} {float(y)!r} {float(z)!r}" for pid, (y, z) in sec.points.items()]
        walls += [f"{sid} {w.wall_id} {w.start} {w.end} {w.laminate_id}" for w in sec.walls]
    out = pts + [""] + walls
    if n_sub:
        out += ["", "[Section_Mesh]", " ".join(MESH_SUBHEADER)] + [f"{sid} {int(n)}" for sid, n in n_sub.items()]
    return "\n".join(out) + "\n"
