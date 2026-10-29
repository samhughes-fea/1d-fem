# pre_processing/parsing/element_section_parser.py
"""Parse ``element_section.txt``: which thin-walled section each element uses.

File layout (comments start with ``#``)::

    [Element_Section]
    [element_id] [section_id]
    1  1
    2  1

Elements not listed keep the classical ``section.txt`` properties. See
docs/conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md.
"""

from __future__ import annotations

from typing import Dict

from pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match

ELEMENT_SECTION_HEADER = "[element_section]"
ELEMENT_SECTION_SUBHEADER = ["[element_id]", "[section_id]"]


def parse_element_section_file(filepath: str) -> Dict[int, int]:
    """Return ``{element_id: section_id}``."""
    blocks = split_blocks(preprocess_lines(filepath), [ELEMENT_SECTION_HEADER])
    if ELEMENT_SECTION_HEADER not in blocks:
        raise ValueError(f"Missing [Element_Section] block in {filepath}")
    sub, rows = blocks[ELEMENT_SECTION_HEADER]
    if not tokens_match(sub, ELEMENT_SECTION_SUBHEADER):
        raise ValueError(
            f"[Element_Section] sub-header must be: {' '.join(ELEMENT_SECTION_SUBHEADER)}"
        )
    out: Dict[int, int] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 2:
            raise ValueError(f"Malformed [Element_Section] row (expected 2 columns): {ln!r}")
        eid, sid = to_int(parts[0], "element_id", ln), to_int(parts[1], "section_id", ln)
        if eid in out:
            raise ValueError(f"Duplicate element_id in [Element_Section]: {eid}")
        out[eid] = sid
    return out
