# pre_processing/parsing/wall_load_parser.py
"""Parse ``wall_load.txt``: uniform surface loads on the walls of thin-walled sections.

File layout (comments start with ``#``)::

    [Wall_Load]
    [element_id] [wall_id] [p_x] [p_s] [p_n]
    0  1  0.0  0.0  -2000.0

``p_x``, ``p_s`` and ``p_n`` are force per area along the beam axis, the wall
tangent (start to end point) and the wall normal ``n = (-t_z, t_y)``. Several
rows for the same element and wall add up. Only elements with a thin-walled
section (``element_section.txt``) may carry wall loads.
"""

from __future__ import annotations

from typing import Dict, Tuple

from pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_float, to_int, tokens_match

WALL_LOAD_HEADER = "[wall_load]"
WALL_LOAD_SUBHEADER = ["[element_id]", "[wall_id]", "[p_x]", "[p_s]", "[p_n]"]


def parse_wall_load_file(filepath: str) -> Dict[int, Dict[int, Tuple[float, float, float]]]:
    """Return ``{element_id: {wall_id: (p_x, p_s, p_n)}}``."""
    blocks = split_blocks(preprocess_lines(filepath), [WALL_LOAD_HEADER])
    if WALL_LOAD_HEADER not in blocks:
        raise ValueError(f"Missing [Wall_Load] block in {filepath}")
    sub, rows = blocks[WALL_LOAD_HEADER]
    if not tokens_match(sub, WALL_LOAD_SUBHEADER):
        raise ValueError(f"[Wall_Load] sub-header must be: {' '.join(WALL_LOAD_SUBHEADER)}")
    out: Dict[int, Dict[int, Tuple[float, float, float]]] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 5:
            raise ValueError(f"Malformed [Wall_Load] row (expected 5 columns): {ln!r}")
        eid, wid = to_int(parts[0], "element_id", ln), to_int(parts[1], "wall_id", ln)
        p = tuple(to_float(parts[k], name, ln) for k, name in ((2, "p_x"), (3, "p_s"), (4, "p_n")))
        old = out.setdefault(eid, {}).get(wid, (0.0, 0.0, 0.0))
        out[eid][wid] = tuple(a + b for a, b in zip(old, p))
    return out
