# pre_processing/parsing/section_forces_parser.py
"""Parse ``section_forces.txt``: section forces from outside the job, for ply checks.

For load cases whose section forces come from elsewhere (an aeroelastic tool, a load
envelope from a certification run, hand calculations). With ``ply_strength.txt`` the job
checks the plies of every listed element's section under each set of forces and writes the
envelope over the load cases (``post_processing/failure/section_force_checks.py``). The
beam solve plays no part in it.

File layout (comments start with ``#``)::

    [Section_Forces]
    [load_case] [element_id] [x] [frame] [f1] [f2] [f3] [f4] [f5] [f6]
    dlc1.1_12ms  0  0.25  engine   1.2e5 3.4e5 -2.0e4 1.1e4 6.0e3 8.0e2
    dlc1.1_12ms  1  0.75  beamdyn  ...

- ``load_case`` is any name without spaces;
- ``element_id`` picks the section (the element's, from ``element_section.txt``);
- ``x`` is where along the beam the forces act, a label for the output only;
- ``frame`` is the frame the six forces are written in, as in ``section_matrix.txt``:

  - ``engine``: 1d-fem's ``[N, M_y, M_z, V_y, V_z, T]``, beam axis x, right-handed;
  - ``public``: anisotropic-section-model's ``[F_x, F_y, N, M_x, M_y, M_z]``, beam axis z;
  - ``beamdyn``: the same order in BeamDyn's mirrored blade frame.

The forces are about the beam axis through the element nodes (the section origin).

Public API
----------
SECTION_FORCES_FILE, FRAMES, parse_section_forces_file, format_section_forces_file, to_engine_forces.
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

from onedfem.pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match

SECTION_FORCES_FILE = "section_forces.txt"
SECTION_FORCES_HEADER = "[section_forces]"
SECTION_FORCES_SUBHEADER = ["[load_case]", "[element_id]", "[x]", "[frame]", "[f1]", "[f2]", "[f3]", "[f4]", "[f5]",
                            "[f6]"]
FRAMES = ("engine", "public", "beamdyn")

#: one row: (load_case, element_id, x, frame, forces (6,))
SectionForceRow = Tuple[str, int, float, str, np.ndarray]


def parse_section_forces_file(filepath: str) -> List[SectionForceRow]:
    """The rows of ``section_forces.txt``, in file order."""
    blocks = split_blocks(preprocess_lines(filepath), [SECTION_FORCES_HEADER])
    if SECTION_FORCES_HEADER not in blocks:
        raise ValueError(f"Missing [Section_Forces] block in {filepath}")
    sub, rows = blocks[SECTION_FORCES_HEADER]
    if not tokens_match(sub, SECTION_FORCES_SUBHEADER):
        raise ValueError(f"[Section_Forces] sub-header must be: {' '.join(SECTION_FORCES_SUBHEADER)}")
    out: List[SectionForceRow] = []
    seen = set()
    for ln in rows:
        tok = ln.split()
        if len(tok) != 10:
            raise ValueError(f"Malformed [Section_Forces] row (expected 10 columns): {ln!r}")
        case, eid, frame = tok[0], to_int(tok[1], "element_id", ln), tok[3].lower()
        if frame not in FRAMES:
            raise ValueError(f"[Section_Forces] frame must be one of {FRAMES}: {ln!r}")
        try:
            x = float(tok[2])
            f = np.array([float(v) for v in tok[4:]])
        except ValueError as exc:
            raise ValueError(f"Non-numeric value in [Section_Forces] row: {ln!r}") from exc
        if (case, eid, x) in seen:
            raise ValueError(f"[Section_Forces] load case {case!r}, element {eid} at x = {x} is given twice")
        seen.add((case, eid, x))
        out.append((case, eid, x, frame, f))
    if not out:
        raise ValueError(f"[Section_Forces] in {filepath} has no rows")
    return out


def to_engine_forces(forces, frame: str) -> np.ndarray:
    """Six section forces in ``frame`` to 1d-fem's right-handed ``[N, M_y, M_z, V_y, V_z, T]``."""
    f = np.asarray(forces, dtype=np.float64)
    if frame == "engine":
        return f.copy()
    from onedfem.sections import section_module

    axes = section_module("axes")
    if frame == "beamdyn":
        f = axes.MIRROR @ f          # MIRROR is its own inverse
    elif frame != "public":
        raise ValueError(f"frame must be one of {FRAMES}, got {frame!r}")
    return axes.to_engine(f)


def format_section_forces_file(rows) -> str:
    """The text of ``section_forces.txt`` for rows ``(load_case, element_id, x, frame, forces)``."""
    lines = ["[Section_Forces]", " ".join(SECTION_FORCES_SUBHEADER)]
    for case, eid, x, frame, f in rows:
        lines.append(f"{case} {eid} {x:.17g} {frame} " + " ".join(f"{v:.17g}" for v in np.asarray(f, float)))
    return "\n".join(lines) + "\n"
