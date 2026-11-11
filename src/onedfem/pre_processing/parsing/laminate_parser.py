# pre_processing/parsing/laminate_parser.py
"""Parse ``laminate.txt``: ply materials and wall laminates for thin-walled sections.

File layout (comments start with ``#``)::

    [Ply_Material]
    [material_id] [E1] [E2] [G12] [nu12] [G13] [G23]          # optional 8th column [rho]
    1  140e9 10e9 5e9 0.3 5e9 3.5e9

    [Laminate]
    [laminate_id] [ply] [material_id] [angle_deg] [thickness]
    1  1  1   0.0  0.125e-3
    1  2  1  90.0  0.125e-3

    [Laminate_Offset]                                         # optional
    [laminate_id] [offset]
    1  -0.004

Plies are numbered 1..n within each laminate, bottom ply first. ``[Laminate_Offset]``
gives the height of a laminate's mid-plane above its wall's reference line (along the
wall normal); laminates not listed have none. See
docs/conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from onedfem.sections import section_import

PLY_MATERIAL_HEADER = "[ply_material]"
LAMINATE_HEADER = "[laminate]"
LAMINATE_OFFSET_HEADER = "[laminate_offset]"

PLY_MATERIAL_SUBHEADER = ["[material_id]", "[E1]", "[E2]", "[G12]", "[nu12]", "[G13]", "[G23]"]
PLY_MATERIAL_SUBHEADER_RHO = PLY_MATERIAL_SUBHEADER + ["[rho]"]
LAMINATE_SUBHEADER = ["[laminate_id]", "[ply]", "[material_id]", "[angle_deg]", "[thickness]"]
LAMINATE_OFFSET_SUBHEADER = ["[laminate_id]", "[offset]"]


def preprocess_lines(filepath: str) -> List[str]:
    """Stripped non-empty lines with ``#`` comments (whole-line or trailing) removed."""
    out: List[str] = []
    with open(filepath, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.split("#", 1)[0].strip()
            if ln:
                out.append(ln)
    return out


def split_blocks(lines: Sequence[str], headers: Sequence[str]) -> Dict[str, Tuple[str, List[str]]]:
    """Split lines into ``{header: (subheader_line, data_lines)}`` for the known block headers.

    A block header is a line equal (case-insensitive) to one of ``headers``. Each
    header may appear once and must be followed by a sub-header line.
    """
    known = {h.lower() for h in headers}
    blocks: Dict[str, Tuple[str, List[str]]] = {}
    current = None
    for ln in lines:
        key = ln.lower()
        if key in known:
            if key in blocks:
                raise ValueError(f"Duplicate {ln} block")
            blocks[key] = ("", [])
            current = key
            continue
        if current is None:
            raise ValueError(f"Line outside any known block {sorted(known)}: {ln!r}")
        sub, rows = blocks[current]
        if not sub:
            blocks[current] = (ln, rows)
        else:
            rows.append(ln)
    for key, (sub, _) in blocks.items():
        if not sub:
            raise ValueError(f"Block {key} has no sub-header")
    return blocks


def tokens_match(line: str, expected: Sequence[str]) -> bool:
    return [t.lower() for t in line.split()] == [e.lower() for e in expected]


def to_float(tok: str, what: str, line: str) -> float:
    try:
        return float(tok)
    except ValueError as exc:
        raise TypeError(f"Bad {what} in line {line!r}: {exc}") from exc


def to_int(tok: str, what: str, line: str) -> int:
    try:
        return int(tok)
    except ValueError as exc:
        raise TypeError(f"Bad {what} in line {line!r}: {exc}") from exc


def parse_ply_materials(sub: str, rows: Sequence[str]) -> Dict[int, PlyMaterial]:
    PlyMaterial = section_import("thin_walled.laminate", "PlyMaterial")
    if tokens_match(sub, PLY_MATERIAL_SUBHEADER_RHO):
        ncol = 8
    elif tokens_match(sub, PLY_MATERIAL_SUBHEADER):
        ncol = 7
    else:
        raise ValueError(
            "[Ply_Material] sub-header must be (case-insensitive): "
            f"{' '.join(PLY_MATERIAL_SUBHEADER)} [optional [rho]]"
        )
    mats: Dict[int, PlyMaterial] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != ncol:
            raise ValueError(f"Malformed [Ply_Material] row (expected {ncol} columns): {ln!r}")
        mid = to_int(parts[0], "material_id", ln)
        if mid in mats:
            raise ValueError(f"Duplicate material_id: {mid}")
        E1, E2, G12, nu12, G13, G23 = (to_float(t, "material value", ln) for t in parts[1:7])
        rho = to_float(parts[7], "rho", ln) if ncol == 8 else 0.0
        mats[mid] = PlyMaterial(E1=E1, E2=E2, G12=G12, nu12=nu12, G13=G13, G23=G23, rho=rho)
    if not mats:
        raise ValueError("[Ply_Material] block has no rows")
    return mats


def parse_laminate_offsets(sub: str, rows: Sequence[str]) -> Dict[int, float]:
    """``{laminate_id: offset}`` of a ``[Laminate_Offset]`` block."""
    if not tokens_match(sub, LAMINATE_OFFSET_SUBHEADER):
        raise ValueError(
            f"[Laminate_Offset] sub-header must be (case-insensitive): {' '.join(LAMINATE_OFFSET_SUBHEADER)}"
        )
    out: Dict[int, float] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 2:
            raise ValueError(f"Malformed [Laminate_Offset] row (expected 2 columns): {ln!r}")
        lid = to_int(parts[0], "laminate_id", ln)
        if lid in out:
            raise ValueError(f"Duplicate [Laminate_Offset] laminate_id: {lid}")
        out[lid] = to_float(parts[1], "offset", ln)
    return out


def parse_laminates(sub: str, rows: Sequence[str], materials: Dict[int, PlyMaterial],
                    offsets: Dict[int, float] | None = None) -> Dict[int, WallLaminate]:
    Ply, WallLaminate = section_import("thin_walled.laminate", "Ply", "WallLaminate")
    if not tokens_match(sub, LAMINATE_SUBHEADER):
        raise ValueError(
            f"[Laminate] sub-header must be (case-insensitive): {' '.join(LAMINATE_SUBHEADER)}"
        )
    by_lam: Dict[int, Dict[int, Ply]] = {}
    for ln in rows:
        parts = ln.split()
        if len(parts) != 5:
            raise ValueError(f"Malformed [Laminate] row (expected 5 columns): {ln!r}")
        lid = to_int(parts[0], "laminate_id", ln)
        k = to_int(parts[1], "ply", ln)
        mid = to_int(parts[2], "material_id", ln)
        angle = to_float(parts[3], "angle_deg", ln)
        t = to_float(parts[4], "thickness", ln)
        if mid not in materials:
            raise ValueError(f"Laminate {lid} ply {k} uses unknown material_id {mid}")
        plies = by_lam.setdefault(lid, {})
        if k in plies:
            raise ValueError(f"Laminate {lid} has duplicate ply {k}")
        plies[k] = Ply(materials[mid], angle, t)
    if not by_lam:
        raise ValueError("[Laminate] block has no rows")
    out: Dict[int, WallLaminate] = {}
    for lid, plies in by_lam.items():
        order = sorted(plies)
        if order != list(range(1, len(order) + 1)):
            raise ValueError(f"Laminate {lid} plies must be numbered 1..{len(order)}, got {order}")
        offset = (offsets or {}).get(lid, 0.0)
        out[lid] = WallLaminate([plies[k] for k in order], offset=offset) if offset else WallLaminate(
            [plies[k] for k in order])
    unknown = sorted(set(offsets or {}) - set(out))
    if unknown:
        raise ValueError(f"[Laminate_Offset] names unknown laminate_id {unknown}")
    return out


def parse_laminate_file(filepath: str) -> Dict[int, WallLaminate]:
    """Read ``laminate.txt`` and return ``{laminate_id: WallLaminate}``."""
    blocks = split_blocks(preprocess_lines(filepath), [PLY_MATERIAL_HEADER, LAMINATE_HEADER, LAMINATE_OFFSET_HEADER])
    for h in (PLY_MATERIAL_HEADER, LAMINATE_HEADER):
        if h not in blocks:
            raise ValueError(f"Missing {h} block in {filepath}")
    materials = parse_ply_materials(*blocks[PLY_MATERIAL_HEADER])
    offsets = parse_laminate_offsets(*blocks[LAMINATE_OFFSET_HEADER]) if LAMINATE_OFFSET_HEADER in blocks else None
    return parse_laminates(*blocks[LAMINATE_HEADER], materials, offsets)


def parse_ply_material_file(filepath: str) -> Dict[int, PlyMaterial]:
    """Read only the ``[Ply_Material]`` block of ``laminate.txt``: ``{material_id: PlyMaterial}``."""
    blocks = split_blocks(preprocess_lines(filepath), [PLY_MATERIAL_HEADER, LAMINATE_HEADER, LAMINATE_OFFSET_HEADER])
    if PLY_MATERIAL_HEADER not in blocks:
        raise ValueError(f"Missing {PLY_MATERIAL_HEADER} block in {filepath}")
    return parse_ply_materials(*blocks[PLY_MATERIAL_HEADER])


def format_laminate_file(laminates: Dict[int, WallLaminate]) -> str:
    """Text of a ``laminate.txt`` for ``{laminate_id: WallLaminate}``.

    Ply materials are numbered from 1 in order of first use (equal materials share an
    id). The file has no shear correction column, so every laminate must use the
    default 5/6. Laminates with an ``offset`` are listed in a ``[Laminate_Offset]``
    block.
    """
    DEFAULT_SHEAR_CORRECTION = section_import("thin_walled.laminate", "DEFAULT_SHEAR_CORRECTION")

    mat_id: Dict[PlyMaterial, int] = {}
    rows = []
    for lid, lam in laminates.items():
        if abs(lam.shear_correction - DEFAULT_SHEAR_CORRECTION) > 1e-12:
            raise ValueError(f"laminate {lid}: laminate.txt has no shear correction column (needs 5/6)")
        for k, ply in enumerate(lam.plies, start=1):
            mid = mat_id.setdefault(ply.material, len(mat_id) + 1)
            rows.append(f"{lid} {k} {mid} {float(ply.angle_deg)!r} {float(ply.thickness)!r}\n")
    mats = "".join(
        f"{mid} " + " ".join(repr(float(v)) for v in (m.E1, m.E2, m.G12, m.nu12, m.G13, m.G23, m.rho)) + "\n"
        for m, mid in mat_id.items())
    offsets = "".join(f"{lid} {float(lam.offset)!r}\n" for lid, lam in laminates.items()
                      if float(getattr(lam, "offset", 0.0)) != 0.0)
    tail = ("\n[Laminate_Offset]\n" + " ".join(LAMINATE_OFFSET_SUBHEADER) + "\n" + offsets) if offsets else ""
    return ("[Ply_Material]\n" + " ".join(PLY_MATERIAL_SUBHEADER_RHO) + "\n" + mats
            + "\n[Laminate]\n" + " ".join(LAMINATE_SUBHEADER) + "\n" + "".join(rows) + tail)
