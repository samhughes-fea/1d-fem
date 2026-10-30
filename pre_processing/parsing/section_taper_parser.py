# pre_processing/parsing/section_taper_parser.py
"""Parse ``section_taper.txt``: the taper coupling of each element's section (the Resal effect).

File layout (comments start with ``#``)::

    [Section_Taper]
    [element_id] [C_eps_y] [C_eps_z] [C_ky_y] [C_ky_z] [C_kz_y] [C_kz_z] [Q_yy] [Q_yz] [Q_zz]
    0  1.2e7  0.0  0.0  3.4e6  ...

``C`` couples the right-handed section strains ``(eps_x, kappa_y, kappa_z)`` with the shear
strains ``(gamma_xy, gamma_xz)`` and ``Q`` adds to the shear block
(``thin_walled.taper``). Only elements with a section from ``element_section.txt`` take a row.
"""

from __future__ import annotations

from typing import Dict

from pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match

SECTION_TAPER_HEADER = "[section_taper]"
SECTION_TAPER_SUBHEADER = ["[element_id]", "[c_eps_y]", "[c_eps_z]", "[c_ky_y]", "[c_ky_z]", "[c_kz_y]",
                           "[c_kz_z]", "[q_yy]", "[q_yz]", "[q_zz]"]


def parse_section_taper_file(filepath: str) -> Dict[int, "SectionTaper"]:
    """Return ``{element_id: SectionTaper}``."""
    from pre_processing.section_integration.thin_walled.taper import SectionTaper

    blocks = split_blocks(preprocess_lines(filepath), [SECTION_TAPER_HEADER])
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
    return out


def format_section_taper_file(tapers: Dict[int, "SectionTaper"]) -> str:
    """Text of ``section_taper.txt`` for ``{element_id: SectionTaper}``."""
    lines = ["[Section_Taper]", " ".join(SECTION_TAPER_SUBHEADER).replace("[c_", "[C_").replace("[q_", "[Q_")]
    for eid in sorted(tapers):
        lines.append(f"{int(eid)} " + " ".join(repr(float(v)) for v in tapers[eid].as_row()))
    return "\n".join(lines) + "\n"
