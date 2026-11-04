# pre_processing/parsing/element_section_parser.py
"""Parse ``element_section.txt``: which thin-walled section each element uses.

File layout (comments start with ``#``)::

    [Element_Section]
    [element_id] [section_id]
    1  1
    2  1

An optional third column ``[modes]`` gives the number of cross-section
distortional modes the element carries (``LinearGBTBeamElement3D``), and an
optional fourth ``[shear_modes]`` the number of shear-lag warping modes::

    [element_id] [section_id] [modes] [shear_modes]
    1  1  2  1

Elements not listed keep the classical ``section.txt`` properties. See
docs/conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md.
"""

from __future__ import annotations

from typing import Dict, Tuple

from pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match

ELEMENT_SECTION_HEADER = "[element_section]"
ELEMENT_SECTION_SUBHEADER = ["[element_id]", "[section_id]"]
ELEMENT_SECTION_MODES_SUBHEADER = ELEMENT_SECTION_SUBHEADER + ["[modes]"]
ELEMENT_SECTION_SHEAR_SUBHEADER = ELEMENT_SECTION_MODES_SUBHEADER + ["[shear_modes]"]


def parse_element_section_file(filepath: str) -> Dict[int, int]:
    """Return ``{element_id: section_id}``."""
    return {e: v[0] for e, v in parse_element_section_modes(filepath).items()}


def parse_element_section_modes(filepath: str) -> Dict[int, Tuple[int, int, int]]:
    """Return ``{element_id: (section_id, modes, shear_modes)}`` (0 for absent columns)."""
    blocks = split_blocks(preprocess_lines(filepath), [ELEMENT_SECTION_HEADER])
    if ELEMENT_SECTION_HEADER not in blocks:
        raise ValueError(f"Missing [Element_Section] block in {filepath}")
    sub, rows = blocks[ELEMENT_SECTION_HEADER]
    if tokens_match(sub, ELEMENT_SECTION_SUBHEADER):
        n_cols = 2
    elif tokens_match(sub, ELEMENT_SECTION_MODES_SUBHEADER):
        n_cols = 3
    elif tokens_match(sub, ELEMENT_SECTION_SHEAR_SUBHEADER):
        n_cols = 4
    else:
        raise ValueError(
            f"[Element_Section] sub-header must be: {' '.join(ELEMENT_SECTION_SUBHEADER)} "
            f"(optionally followed by [modes] and [shear_modes])"
        )
    out: Dict[int, Tuple[int, int, int]] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != n_cols:
            raise ValueError(f"Malformed [Element_Section] row (expected {n_cols} columns): {ln!r}")
        eid, sid = to_int(parts[0], "element_id", ln), to_int(parts[1], "section_id", ln)
        modes = to_int(parts[2], "modes", ln) if n_cols >= 3 else 0
        shear = to_int(parts[3], "shear_modes", ln) if n_cols == 4 else 0
        if modes < 0 or shear < 0:
            raise ValueError(f"[modes] and [shear_modes] must be >= 0: {ln!r}")
        if eid in out:
            raise ValueError(f"Duplicate element_id in [Element_Section]: {eid}")
        out[eid] = (sid, modes, shear)
    return out
