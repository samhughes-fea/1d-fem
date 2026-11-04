# pre_processing/parsing/ply_strength_parser.py
"""Parse ``ply_strength.txt``: ply strengths per material id, for the failure criteria.

Format::

    [Ply_Strength]
    [material_id] [Xt] [Xc] [Yt] [Yc] [S12] [S23] [p_TL_t] [p_TL_c] [p_TT]
    1 1500e6 900e6 27e6 200e6 80e6 - 0.35 0.30 0.25

The first six columns are required. Optional columns, in any order after them, are any of
``[S23] [Zt] [Zc] [F_star] [p_TL_t] [p_TL_c] [p_TT]`` (see ``PlyStrength``); ``-`` keeps
the default. Material ids are those of the job's material tables: ``[Ply_Material]`` in
``laminate.txt`` and ``[Solid_Material]`` in ``section_mesh_2d.txt``.

Public API
----------
PLY_STRENGTH_FILE, parse_ply_strength_file, format_ply_strength_file.
"""

from __future__ import annotations

import os
from typing import Dict, Mapping

from pre_processing.section_integration.ply_failure import PlyStrength

PLY_STRENGTH_FILE = "ply_strength.txt"
_REQUIRED = ("material_id", "Xt", "Xc", "Yt", "Yc", "S12")
_OPTIONAL = ("S23", "Zt", "Zc", "F_star", "p_TL_t", "p_TL_c", "p_TT")


def parse_ply_strength_file(path: str) -> Dict[int, PlyStrength]:
    """``{material_id: PlyStrength}`` from ``path``; an empty dict when the file is absent."""
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        lines = [ln.split("#")[0].strip() for ln in f]
    lines = [ln for ln in lines if ln]
    if not lines or lines[0].lower() != "[ply_strength]":
        raise ValueError(f"{path}: the first line must be [Ply_Strength]")
    if len(lines) < 2 or not lines[1].startswith("["):
        raise ValueError(f"{path}: the second line must name the columns, e.g. [material_id] [Xt] [Xc] ...")
    cols = [c.strip("[]") for c in lines[1].split()]
    if tuple(cols[:6]) != _REQUIRED:
        raise ValueError(f"{path}: the first columns must be {' '.join(f'[{c}]' for c in _REQUIRED)}")
    unknown = [c for c in cols[6:] if c not in _OPTIONAL]
    if unknown:
        raise ValueError(f"{path}: unknown columns {unknown}; optional columns are {list(_OPTIONAL)}")
    out: Dict[int, PlyStrength] = {}
    for ln in lines[2:]:
        tok = ln.split()
        if len(tok) != len(cols):
            raise ValueError(f"{path}: expected {len(cols)} values, got {len(tok)} in '{ln}'")
        mid = int(tok[0])
        if mid in out:
            raise ValueError(f"{path}: material {mid} is given twice")
        kw = {c: float(t) for c, t in zip(cols[1:], tok[1:]) if t != "-"}
        out[mid] = PlyStrength(**kw)
    return out


def format_ply_strength_file(strengths: Mapping[int, PlyStrength]) -> str:
    """Text of ``ply_strength.txt`` with every column."""
    cols = _REQUIRED + _OPTIONAL
    rows = ["[Ply_Strength]", " ".join(f"[{c}]" for c in cols)]
    for mid, st in sorted(strengths.items()):
        vals = [getattr(st, c) for c in cols[1:]]
        rows.append(" ".join([str(int(mid))] + ["-" if v is None else repr(float(v)) for v in vals]))
    return "\n".join(rows) + "\n"
