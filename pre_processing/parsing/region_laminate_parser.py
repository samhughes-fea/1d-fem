# pre_processing/parsing/region_laminate_parser.py
"""Parse ``region_laminate.txt``: which ``laminate.txt`` laminate each geometry region uses.

Sections cut from an implicit-geometry-kernel blade name their walls by region label
(``skin_suction_surface:bay_0``, ``spar_cap_pressure_surface``, ``shear_web_1``). This
file binds labels to laminate ids. A wall looks up its full label first, then its
region (the label without the ``:part`` suffix), then its kind (``skin``, ``cap``,
``web``), so one line per kind is enough for a uniform blade::

    [Region_Laminate]
    [region] [laminate_id]
    skin                          1
    cap                           2
    web                           3
    skin_suction_surface:bay_2    4      # a core panel aft of the rear web

See ``pre_processing/section_integration/thin_walled/from_geometry.py`` and
docs/conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md.
"""

from __future__ import annotations

from typing import Dict, Iterable, Optional

from pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match

REGION_LAMINATE_HEADER = "[region_laminate]"
REGION_LAMINATE_SUBHEADER = ["[region]", "[laminate_id]"]


def parse_region_laminate_file(filepath: str, laminate_ids: Optional[Iterable[int]] = None) -> Dict[str, int]:
    """Read ``region_laminate.txt`` and return ``{label: laminate_id}``.

    With ``laminate_ids`` every id must be one of them.
    """
    blocks = split_blocks(preprocess_lines(filepath), [REGION_LAMINATE_HEADER])
    if REGION_LAMINATE_HEADER not in blocks:
        raise ValueError(f"Missing [Region_Laminate] block in {filepath}")
    sub, rows = blocks[REGION_LAMINATE_HEADER]
    if not tokens_match(sub, REGION_LAMINATE_SUBHEADER):
        raise ValueError(f"[Region_Laminate] sub-header must be: {' '.join(REGION_LAMINATE_SUBHEADER)}")
    known = None if laminate_ids is None else {int(i) for i in laminate_ids}
    out: Dict[str, int] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 2:
            raise ValueError(f"Malformed [Region_Laminate] row (expected 2 columns): {ln!r}")
        label, lid = parts[0], to_int(parts[1], "laminate_id", ln)
        if label in out:
            raise ValueError(f"Duplicate region {label!r} in {filepath}")
        if known is not None and lid not in known:
            raise ValueError(f"region {label!r} uses unknown laminate_id {lid}")
        out[label] = lid
    if not out:
        raise ValueError("[Region_Laminate] block has no rows")
    return out


def format_region_laminate(mapping: Dict[str, int]) -> str:
    """Text of a ``region_laminate.txt`` for ``{label: laminate_id}``."""
    width = max(len(k) for k in mapping) + 2
    rows = "".join(f"{label:<{width}}{lid}\n" for label, lid in mapping.items())
    return "[Region_Laminate]\n[region] [laminate_id]\n" + rows
