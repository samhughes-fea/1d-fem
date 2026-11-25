# pre_processing/parsing/section_taper_parser.py
"""Parse ``section_taper.txt``: the taper coupling of each element's section (the Resal effect).

File layout (comments start with ``#``)::

    [Section_Taper]
    [element_id] [C_eps_y] [C_eps_z] [C_ky_y] [C_ky_z] [C_kz_y] [C_kz_z] [Q_yy] [Q_yz] [Q_zz]
    0  1.2e7  0.0  0.0  3.4e6  ...
    [Section_Taper_Laminate]
    [element_id] [laminate_id] [C_eps_y] [C_eps_z] [C_ky_y] [C_ky_z] [C_kz_y] [C_kz_z] [Q_yy] [Q_yz] [Q_zz]
    0  2  0.4  0.0  ...
    [Section_Taper_2D]
    [element_id] [section_start] [section_end] [distance]
    1  11  12  0.02

``C`` couples the right-handed section strains ``(eps_x, kappa_y, kappa_z)`` with the shear
strains ``(gamma_xy, gamma_xz)`` and ``Q`` adds to the shear block
(``thin_walled.taper``). Only elements with a section from ``element_section.txt`` take a row.

The optional ``[Section_Taper_Laminate]`` block splits each element's coupling by laminate,
per unit axial wall stiffness (``SectionTaper.parts``): layup sensitivities of a tapered
beam need it.

The optional ``[Section_Taper_2D]`` block gives an element the taper coupling of a section
meshed in 2D instead: ``section_start`` and ``section_end`` are two 2D sections of the job
(``[Section_Mesh_2D]`` in ``section_geometry.txt`` or ``section_mesh_2d.txt``) of the same
mesh topology, ``distance`` apart along the beam, ``section_end`` further along ``x``; the
coupling is the provider's ``solid.taper.solid_section_taper`` of the section midway between
them. It is differentiated through the two meshes, so layup sensitivities of a tapered 2D
section need it. An element takes a row in ``[Section_Taper]`` or here, not both, and a
file may hold either block alone.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

from onedfem.pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match
from onedfem.sections import section_import

SECTION_TAPER_HEADER = "[section_taper]"
SECTION_TAPER_SUBHEADER = ["[element_id]", "[c_eps_y]", "[c_eps_z]", "[c_ky_y]", "[c_ky_z]", "[c_kz_y]",
                           "[c_kz_z]", "[q_yy]", "[q_yz]", "[q_zz]"]
SECTION_TAPER_LAMINATE_HEADER = "[section_taper_laminate]"
SECTION_TAPER_LAMINATE_SUBHEADER = SECTION_TAPER_SUBHEADER[:1] + ["[laminate_id]"] + SECTION_TAPER_SUBHEADER[1:]
SECTION_TAPER_2D_HEADER = "[section_taper_2d]"
SECTION_TAPER_2D_SUBHEADER = ["[element_id]", "[section_start]", "[section_end]", "[distance]"]
_HEADERS = [SECTION_TAPER_HEADER, SECTION_TAPER_LAMINATE_HEADER, SECTION_TAPER_2D_HEADER]


def parse_section_taper_file(filepath: str) -> Dict[int, "SectionTaper"]:
    """Return ``{element_id: SectionTaper}`` of the ``[Section_Taper]`` block (empty without it).

    The ``[Section_Taper_2D]`` rows are ``parse_section_taper_2d``'s.
    """
    SectionTaper = section_import("thin_walled.taper", "SectionTaper")

    blocks = split_blocks(preprocess_lines(filepath), _HEADERS)
    if SECTION_TAPER_HEADER not in blocks and SECTION_TAPER_2D_HEADER not in blocks:
        raise ValueError(f"Missing [Section_Taper] (or [Section_Taper_2D]) block in {filepath}")
    if SECTION_TAPER_HEADER not in blocks:
        if SECTION_TAPER_LAMINATE_HEADER in blocks:
            raise ValueError("[Section_Taper_Laminate] splits [Section_Taper] rows; the file has no [Section_Taper] block")
        return {}
    sub, rows = blocks[SECTION_TAPER_HEADER]
    if not tokens_match(sub, SECTION_TAPER_SUBHEADER):
        raise ValueError(f"[Section_Taper] sub-header must be: {' '.join(SECTION_TAPER_SUBHEADER)}")
    out: Dict[int, SectionTaper] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 10:
            raise ValueError(f"Malformed [Section_Taper] row (expected 10 columns): {ln!r}")
        eid = to_int(parts[0], "element_id", ln)
        if eid in out:
            raise ValueError(f"Duplicate element_id in [Section_Taper]: {eid}")
        try:
            vals = [float(v) for v in parts[1:]]
        except ValueError as exc:
            raise ValueError(f"Non-numeric value in [Section_Taper] row: {ln!r}") from exc
        out[eid] = SectionTaper.from_row(vals)
    if SECTION_TAPER_LAMINATE_HEADER in blocks:
        import dataclasses

        import numpy as np

        sub, rows = blocks[SECTION_TAPER_LAMINATE_HEADER]
        if not tokens_match(sub, SECTION_TAPER_LAMINATE_SUBHEADER):
            raise ValueError(f"[Section_Taper_Laminate] sub-header must be: {' '.join(SECTION_TAPER_LAMINATE_SUBHEADER)}")
        parts: Dict[int, Dict[int, np.ndarray]] = {}
        for ln in rows:
            tok = ln.split()
            if len(tok) != 11:
                raise ValueError(f"Malformed [Section_Taper_Laminate] row (expected 11 columns): {ln!r}")
            eid, lid = to_int(tok[0], "element_id", ln), to_int(tok[1], "laminate_id", ln)
            if eid not in out:
                raise ValueError(f"[Section_Taper_Laminate] names element {eid}, which has no [Section_Taper] row")
            parts.setdefault(eid, {})[lid] = np.array([float(v) for v in tok[2:]])
        out = {e: dataclasses.replace(t, parts=parts.get(e)) for e, t in out.items()}
    return out


def parse_section_taper_2d(filepath: str) -> Dict[int, Tuple[int, int, float]]:
    """Return ``{element_id: (section_start, section_end, distance)}`` of the ``[Section_Taper_2D]``
    block (empty without it)."""
    blocks = split_blocks(preprocess_lines(filepath), _HEADERS)
    if SECTION_TAPER_2D_HEADER not in blocks:
        return {}
    sub, rows = blocks[SECTION_TAPER_2D_HEADER]
    if not tokens_match(sub, SECTION_TAPER_2D_SUBHEADER):
        raise ValueError(f"[Section_Taper_2D] sub-header must be: {' '.join(SECTION_TAPER_2D_SUBHEADER)}")
    out: Dict[int, Tuple[int, int, float]] = {}
    for ln in rows:
        tok = ln.split()
        if len(tok) != 4:
            raise ValueError(f"Malformed [Section_Taper_2D] row (expected 4 columns): {ln!r}")
        eid = to_int(tok[0], "element_id", ln)
        if eid in out:
            raise ValueError(f"Duplicate element_id in [Section_Taper_2D]: {eid}")
        start, end = to_int(tok[1], "section_start", ln), to_int(tok[2], "section_end", ln)
        try:
            distance = float(tok[3])
        except ValueError as exc:
            raise ValueError(f"Non-numeric distance in [Section_Taper_2D] row: {ln!r}") from exc
        if not distance > 0.0:
            raise ValueError(f"[Section_Taper_2D] distance must be positive: {ln!r}")
        if start == end:
            raise ValueError(f"[Section_Taper_2D] needs two different sections: {ln!r}")
        out[eid] = (start, end, distance)
    return out


def format_section_taper_file(tapers: Dict[int, "SectionTaper"],
                              tapers_2d: Optional[Dict[int, Tuple[int, int, float]]] = None) -> str:
    """Text of ``section_taper.txt`` for ``{element_id: SectionTaper}`` and, optionally, the
    ``[Section_Taper_2D]`` rows ``{element_id: (section_start, section_end, distance)}``."""
    lines = []
    if tapers:
        lines += ["[Section_Taper]", " ".join(SECTION_TAPER_SUBHEADER).replace("[c_", "[C_").replace("[q_", "[Q_")]
    for eid in sorted(tapers):
        lines.append(f"{int(eid)} " + " ".join(repr(float(v)) for v in tapers[eid].as_row()))
    if any(t.parts for t in tapers.values()):
        lines += ["[Section_Taper_Laminate]",
                  " ".join(SECTION_TAPER_LAMINATE_SUBHEADER).replace("[c_", "[C_").replace("[q_", "[Q_")]
        for eid in sorted(tapers):
            for lid, row in sorted((tapers[eid].parts or {}).items()):
                lines.append(f"{int(eid)} {int(lid)} " + " ".join(repr(float(v)) for v in row))
    if tapers_2d:
        lines += ["[Section_Taper_2D]", " ".join(SECTION_TAPER_2D_SUBHEADER)]
        lines += [f"{int(e)} {int(a)} {int(b)} {float(d)!r}" for e, (a, b, d) in sorted(tapers_2d.items())]
    return "\n".join(lines) + "\n"
