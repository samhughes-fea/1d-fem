# pre_processing/parsing/element_section_parser.py
"""Parse ``element_section.txt``: which thin-walled section each element uses.

File layout (comments start with ``#``)::

    [Element_Section]
    [element_id] [section_id]
    1  1
    2  1

An optional third column ``[modes]`` gives the number of cross-section
distortional modes the element carries (``LinearGBTBeamElement3D``)::

    [element_id] [section_id] [modes]
    1  1  2

Elements not listed keep the classical ``section.txt`` properties. See
docs/conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md.
"""

from __future__ import annotations

from typing import Dict, Tuple

from pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match

ELEMENT_SECTION_HEADER = "[element_section]"
ELEMENT_SECTION_SUBHEADER = ["[element_id]", "[section_id]"]
ELEMENT_SECTION_MODES_SUBHEADER = ELEMENT_SECTION_SUBHEADER + ["[modes]"]


def parse_element_section_file(filepath: str) -> Dict[int, int]:
    """Return ``{element_id: section_id}``."""
    return {e: sid for e, (sid, _) in parse_element_section_modes(filepath).items()}


def parse_element_section_modes(filepath: str) -> Dict[int, Tuple[int, int]]:
    """Return ``{element_id: (section_id, modes)}``; ``modes`` is 0 without the ``[modes]`` column."""
    blocks = split_blocks(preprocess_lines(filepath), [ELEMENT_SECTION_HEADER])
    if ELEMENT_SECTION_HEADER not in blocks:
        raise ValueError(f"Missing [Element_Section] block in {filepath}")
    sub, rows = blocks[ELEMENT_SECTION_HEADER]
    if tokens_match(sub, ELEMENT_SECTION_SUBHEADER):
        n_cols = 2
    elif tokens_match(sub, ELEMENT_SECTION_MODES_SUBHEADER):
        n_cols = 3
    else:
        raise ValueError(
            f"[Element_Section] sub-header must be: {' '.join(ELEMENT_SECTION_SUBHEADER)} "
            f"(optionally followed by [modes])"
        )
    out: Dict[int, Tuple[int, int]] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != n_cols:
            raise ValueError(f"Malformed [Element_Section] row (expected {n_cols} columns): {ln!r}")
        eid, sid = to_int(parts[0], "element_id", ln), to_int(parts[1], "section_id", ln)
        modes = to_int(parts[2], "modes", ln) if n_cols == 3 else 0
        if modes < 0:
            raise ValueError(f"[modes] must be >= 0: {ln!r}")
        if eid in out:
            raise ValueError(f"Duplicate element_id in [Element_Section]: {eid}")
        out[eid] = (sid, modes)
    return out
