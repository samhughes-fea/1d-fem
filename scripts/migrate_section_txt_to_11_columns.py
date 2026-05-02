"""
One-shot migration: jobs/**/section.txt from 6-column to 11-column subheader.

Appends: kappa=5/6, alpha=I_z/A, y_sc=0, z_sc=0, Gamma=0.

Usage (repo root):
    python scripts/migrate_section_txt_to_11_columns.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

SUBHEADER_6_TOK = [
    "[element_id]",
    "[a]",
    "[i_x]",
    "[i_y]",
    "[i_z]",
    "[j_t]",
]

SUBHEADER_11_LINE = (
    "[element_id]  [A]  [I_x]  [I_y]  [I_z]  [J_t]  "
    "[kappa]  [alpha]  [y_sc]  [z_sc]  [Gamma]"
)

KAPPA_VALUE = 5.0 / 6.0


def _norm_header_tokens(line: str) -> list[str]:
    return [t.lower() for t in line.split()]


def migrate_text(text: str) -> str | None:
    """
    Return migrated file text, or None if already 11-column.

    Raises ValueError on unsupported or malformed content.
    """
    lines = text.splitlines()
    try:
        sec_i = next(i for i, ln in enumerate(lines) if ln.strip().lower() == "[section]")
    except StopIteration:
        raise ValueError("no [Section] header") from None

    if sec_i + 1 >= len(lines):
        raise ValueError("[Section] not followed by subheader")

    sub_line = lines[sec_i + 1]
    sub_tok = _norm_header_tokens(sub_line)
    sub11_tok = [
        "[element_id]",
        "[a]",
        "[i_x]",
        "[i_y]",
        "[i_z]",
        "[j_t]",
        "[kappa]",
        "[alpha]",
        "[y_sc]",
        "[z_sc]",
        "[gamma]",
    ]
    if sub_tok == sub11_tok:
        return None
    if sub_tok != SUBHEADER_6_TOK:
        raise ValueError(
            f"expected 6-column subheader, got tokens {sub_line.split()!r}"
        )

    head = lines[: sec_i + 1]
    out_lines = list(head)
    out_lines.append(SUBHEADER_11_LINE)

    for ln in lines[sec_i + 2 :]:
        st = ln.strip()
        if st == "":
            out_lines.append(ln.rstrip("\r\n"))
            continue
        if ln.lstrip().startswith("#"):
            out_lines.append(ln.rstrip("\r\n"))
            continue
        parts = ln.split()
        if len(parts) != 6:
            raise ValueError(f"data row expected 6 columns, got {len(parts)}: {ln!r}")
        eid, a_s, ix_s, iy_s, iz_s, jt_s = parts
        a_f = float(a_s)
        iz_f = float(iz_s)
        alpha = (iz_f / a_f) if a_f > 0.0 else 0.0
        row = (
            f"{eid}  {a_s}  {ix_s}  {iy_s}  {iz_s}  {jt_s}  "
            f"{KAPPA_VALUE:.17g}  {alpha:.17g}  0.0  0.0  0.0"
        )
        out_lines.append(row)

    return "\n".join(out_lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=PROJECT_ROOT,
        help="Repository root (default: parent of scripts/)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions only; do not write files",
    )
    args = parser.parse_args()
    root: Path = args.root.resolve()

    paths = sorted(root.glob("jobs/**/section.txt"))
    migrated = 0
    skipped = 0
    errors: list[str] = []

    for path in paths:
        rel = path.relative_to(root)
        try:
            text = path.read_text(encoding="utf-8")
            new_text = migrate_text(text)
            if new_text is None:
                print(f"skip (already 11 col): {rel}")
                skipped += 1
                continue
            if args.dry_run:
                print(f"would migrate: {rel}")
                migrated += 1
            else:
                path.write_text(new_text, encoding="utf-8", newline="\n")
                print(f"migrated: {rel}")
                migrated += 1
        except Exception as e:
            errors.append(f"{rel}: {e}")

    print(f"Done. migrated={migrated} skipped={skipped} errors={len(errors)}")
    for msg in errors:
        print(msg, file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
