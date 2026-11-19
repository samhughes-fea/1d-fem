# pre_processing/parsing/section_matrix_parser.py
"""Parse ``section_matrix.txt``: sections given only as their 6 x 6 stiffness and mass tables.

For section properties from elsewhere (a BeamDyn blade file, windIO 2.0 elastic properties,
VABS or BECAS output). The section provider builds each one
(``SectionAnalysis.from_matrices``); ``element_section.txt`` assigns them to elements like
thin-walled and 2D sections, and their ids must differ from those.

File layout (comments start with ``#``), six rows of each matrix per section::

    [Section_Matrix]
    [section_id] [frame] [matrix] [row] [c1] [c2] [c3] [c4] [c5] [c6]
    1 beamdyn stiffness 1  1.2e10 0 0 0 0 3.1e8
    ...
    1 beamdyn stiffness 6  ...
    1 beamdyn mass 1  678.9 0 0 0 0 -12.3
    ...

``frame`` is the frame the tables are written in, the same on every row of a section:

- ``engine``: 1d-fem's beam axes (beam axis x, right-handed strains
  ``[eps, kappa_y, kappa_z, gamma_xy, gamma_xz, phi]``, mass on ``[u_x, u_y, u_z, theta_x,
  theta_y, theta_z]``);
- ``public``: anisotropic-section-model's public frame (beam axis z, BeamDyn's order);
- ``beamdyn``: BeamDyn's own mirrored blade frame, as in a BeamDyn blade file.

Both matrices are about the beam axis through the element nodes. The mass rows may be left
out (static analyses only). The mass must be a rigid-section mass (mass, its centre and the
second moments). A table section carries no second-order, warping, third-order or GBT data,
and the tables are taken as complete: pretwisted and curved elements keep them as they are.

Public API
----------
SECTION_MATRIX_FILE, parse_section_matrix_file, format_section_matrix_file, section_matrix_analyses.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from onedfem.pre_processing.parsing.laminate_parser import preprocess_lines, split_blocks, to_int, tokens_match

SECTION_MATRIX_FILE = "section_matrix.txt"
SECTION_MATRIX_HEADER = "[section_matrix]"
SECTION_MATRIX_SUBHEADER = ["[section_id]", "[frame]", "[matrix]", "[row]", "[c1]", "[c2]", "[c3]", "[c4]", "[c5]",
                            "[c6]"]
FRAMES = ("engine", "public", "beamdyn")
MATRICES = ("stiffness", "mass")


def parse_section_matrix_file(filepath: str) -> Dict[int, Tuple[str, np.ndarray, Optional[np.ndarray]]]:
    """Return ``{section_id: (frame, stiffness 6 x 6, mass 6 x 6 or None)}``."""
    blocks = split_blocks(preprocess_lines(filepath), [SECTION_MATRIX_HEADER])
    if SECTION_MATRIX_HEADER not in blocks:
        raise ValueError(f"Missing [Section_Matrix] block in {filepath}")
    sub, rows = blocks[SECTION_MATRIX_HEADER]
    if not tokens_match(sub, SECTION_MATRIX_SUBHEADER):
        raise ValueError(f"[Section_Matrix] sub-header must be: {' '.join(SECTION_MATRIX_SUBHEADER)}")
    frames: Dict[int, str] = {}
    data: Dict[Tuple[int, str], Dict[int, np.ndarray]] = {}
    for ln in rows:
        tok = ln.split()
        if len(tok) != 10:
            raise ValueError(f"Malformed [Section_Matrix] row (expected 10 columns): {ln!r}")
        sid = to_int(tok[0], "section_id", ln)
        frame, matrix = tok[1].lower(), tok[2].lower()
        if frame not in FRAMES:
            raise ValueError(f"[Section_Matrix] frame must be one of {FRAMES}: {ln!r}")
        if matrix not in MATRICES:
            raise ValueError(f"[Section_Matrix] matrix must be one of {MATRICES}: {ln!r}")
        if frames.setdefault(sid, frame) != frame:
            raise ValueError(f"[Section_Matrix] section {sid} is written in two frames ({frames[sid]}, {frame})")
        r = to_int(tok[3], "row", ln)
        if not 1 <= r <= 6:
            raise ValueError(f"[Section_Matrix] row must be 1 to 6: {ln!r}")
        try:
            values = np.array([float(v) for v in tok[4:]])
        except ValueError as exc:
            raise ValueError(f"Non-numeric value in [Section_Matrix] row: {ln!r}") from exc
        rows_of = data.setdefault((sid, matrix), {})
        if r in rows_of:
            raise ValueError(f"[Section_Matrix] section {sid} {matrix} row {r} is given twice")
        rows_of[r] = values
    out: Dict[int, Tuple[str, np.ndarray, Optional[np.ndarray]]] = {}
    for sid, frame in frames.items():
        mats = []
        for matrix in MATRICES:
            rows_of = data.get((sid, matrix))
            if rows_of is None:
                if matrix == "stiffness":
                    raise ValueError(f"[Section_Matrix] section {sid} has no stiffness rows")
                mats.append(None)
                continue
            if sorted(rows_of) != [1, 2, 3, 4, 5, 6]:
                raise ValueError(f"[Section_Matrix] section {sid} {matrix} needs rows 1 to 6, has {sorted(rows_of)}")
            mats.append(np.array([rows_of[r] for r in range(1, 7)]))
        out[sid] = (frame, mats[0], mats[1])
    return out


def format_section_matrix_file(sections: Dict[int, Tuple[str, np.ndarray, Optional[np.ndarray]]]) -> str:
    """The text of ``section_matrix.txt`` for ``{section_id: (frame, stiffness, mass or None)}``."""
    lines = ["[Section_Matrix]", " ".join(SECTION_MATRIX_SUBHEADER)]
    for sid, (frame, K, M) in sorted(sections.items()):
        for name, A in (("stiffness", K), ("mass", M)):
            if A is None:
                continue
            A = np.asarray(A, dtype=float)
            for r in range(6):
                lines.append(f"{sid} {frame} {name} {r + 1} " + " ".join(f"{v:.17g}" for v in A[r]))
    return "\n".join(lines) + "\n"


def section_matrix_analyses(filepath: str) -> Dict[int, object]:
    """``{section_id: SectionAnalysis}`` for every section of ``section_matrix.txt``, built by the
    section provider (``SectionAnalysis.from_matrices``)."""
    from onedfem.sections import section_import

    SectionAnalysis = section_import("section_analysis", "SectionAnalysis")
    out = {}
    for sid, (frame, K, M) in parse_section_matrix_file(filepath).items():
        try:
            out[sid] = SectionAnalysis.from_matrices(K, M, frame=frame, section_id=sid)
        except ValueError as exc:
            raise ValueError(f"section_matrix.txt, section {sid}: {exc}") from exc
    return out
