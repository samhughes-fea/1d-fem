# pre_processing/parsing/section_taper_parser.py
"""Parse ``section_taper.txt``: the taper coupling of each element's section (the Resal effect).

File layout (comments start with ``#``)::

    [Section_Taper]
    [element_id] [C_eps_y] [C_eps_z] [C_ky_y] [C_ky_z] [C_kz_y] [C_kz_z] [Q_yy] [Q_yz] [Q_zz]
    0  1.2e7  0.0  0.0  3.4e6  ...
    [Section_Taper_Laminate]
    [element_id] [laminate_id] [C_eps_y] [C_eps_z] [C_ky_y] [C_ky_z] [C_kz_y] [C_kz_z] [Q_yy] [Q_yz] [Q_zz]
    0  2  0.4  0.0  ...

``C`` couples the right-handed section strains ``(eps_x, kappa_y, kappa_z)`` with the shear
strains ``(gamma_xy, gamma_xz)`` and ``Q`` adds to the shear block
(``thin_walled.taper``). Only elements with a section from ``element_section.txt`` take a row.

The optional ``[Section_Taper_Laminate]`` block splits each element's coupling by laminate,
per unit axial wall stiffness (``SectionTaper.parts``): layup sensitivities of a tapered
beam need it.
"""

from __future__ import annotations

from typing import Dict

from onedfem.pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match
from onedfem.sections import section_import

SECTION_TAPER_HEADER = "[section_taper]"
SECTION_TAPER_SUBHEADER = ["[element_id]", "[c_eps_y]", "[c_eps_z]", "[c_ky_y]", "[c_ky_z]", "[c_kz_y]",
                           "[c_kz_z]", "[q_yy]", "[q_yz]", "[q_zz]"]
SECTION_TAPER_LAMINATE_HEADER = "[section_taper_laminate]"
SECTION_TAPER_LAMINATE_SUBHEADER = SECTION_TAPER_SUBHEADER[:1] + ["[laminate_id]"] + SECTION_TAPER_SUBHEADER[1:]


def parse_section_taper_file(filepath: str) -> Dict[int, "SectionTaper"]:
    """Return ``{element_id: SectionTaper}``."""
    SectionTaper = section_import("thin_walled.taper", "SectionTaper")

    blocks = split_blocks(preprocess_lines(filepath), [SECTION_TAPER_HEADER, SECTION_TAPER_LAMINATE_HEADER])
    if SECTION_TAPER_HEADER not in blocks:
        raise ValueError(f"Missing [Section_Taper] block in {filepath}")
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


def format_section_taper_file(tapers: Dict[int, "SectionTaper"]) -> str:
    """Text of ``section_taper.txt`` for ``{element_id: SectionTaper}``."""
    lines = ["[Section_Taper]", " ".join(SECTION_TAPER_SUBHEADER).replace("[c_", "[C_").replace("[q_", "[Q_")]
    for eid in sorted(tapers):
        lines.append(f"{int(eid)} " + " ".join(repr(float(v)) for v in tapers[eid].as_row()))
    if any(t.parts for t in tapers.values()):
        lines += ["[Section_Taper_Laminate]",
                  " ".join(SECTION_TAPER_LAMINATE_SUBHEADER).replace("[c_", "[C_").replace("[q_", "[Q_")]
        for eid in sorted(tapers):
            for lid, row in sorted((tapers[eid].parts or {}).items()):
                lines.append(f"{int(eid)} {int(lid)} " + " ".join(repr(float(v)) for v in row))
    return "\n".join(lines) + "\n"
