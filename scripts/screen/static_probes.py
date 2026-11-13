"""Probes 2 to 7 of the 1d-fem screen: scans of the library source.

2. ``names``: element type names and name fragments used outside the element library
   (R1, G2). A capability decided by a name list instead of a class attribute.
3. ``private``: ``obj._attr`` reads where ``obj`` is not ``self``/``cls`` and ``_attr`` is
   not defined in the same module (C1).
4. ``domain``: blade, airfoil, wing, spar, chord, rotor, turbine words in the core (G1).
5. ``outputs``: every result file the library writes, with its columns where a module
   declares them, and whether those say where (element, xi, wall, s, ply, z) and what kind
   of value (Gauss, nodal, projected) each row is (P1, P3).
6. ``plots``: smoothing or interpolation calls and figures saved outside ``figure_style``
   in the library, examples and validation (P2).
7. ``coverage``: for each registered element, its theory pages and the tests that name it
   (R2).

Usage::

    python scripts/screen/static_probes.py --out screen_out

Writes ``static_probes.json`` and one CSV per probe.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "onedfem"
ELEMENT_LIBRARY = SRC / "pre_processing" / "element_library"
ADAPTERS = {"pre_processing/mesh_library/blade_beam_job.py",
            "pre_processing/section_integration/thin_walled/from_geometry.py",
            "pre_processing/section_integration/solid/from_geometry.py",
            "interfaces/__init__.py",
            "interfaces/blade/__init__.py",
            "interfaces/blade/blade_beam_job.py",
            "interfaces/blade/thin_walled_sections.py",
            "interfaces/blade/solid_sections.py"}
TYPE_NAME = re.compile(r"^(Linear|Nonlinear|Corotational|UpdatedLagrangian|GEBT)\w*Element3D$")
FRAGMENTS = ("GBT", "Benscoter", "Reddy", "Levinson", "Timoshenko", "EulerBernoulli", "Corotational", "Nonlinear",
             "Truss", "Bar")
DOMAIN = re.compile(r"\b(blade|airfoil|aerofoil|wing|spar|chord|rotor|turbine)s?\b", re.I)
SMOOTHING = re.compile(r"\b(interp1d|make_interp_spline|UnivariateSpline|CubicSpline|splrep|splev|savgol_filter|"
                       r"gaussian_filter1d|np\.interp|numpy\.interp|PchipInterpolator|Akima1DInterpolator)\b")
WHERE = ("element", "xi", "wall", "s", "ply", "z", "n", "x", "node", "station", "gauss")
KIND = ("gauss", "nodal", "projected", "kind", "point_kind")


def _py(base: Path):
    return sorted(p for p in base.rglob("*.py") if "__pycache__" not in p.parts)


def _rel(p: Path) -> str:
    return str(p.relative_to(SRC)) if SRC in p.parents else str(p.relative_to(ROOT))


def _parse(p: Path):
    return ast.parse(p.read_text(encoding="utf-8-sig"), filename=str(p))


def probe_names():
    rows = []
    for p in _py(SRC):
        in_library = ELEMENT_LIBRARY in p.parents
        tree = _parse(p)
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                v = node.value
                if TYPE_NAME.match(v) and not in_library:
                    rows.append({"file": _rel(p), "line": node.lineno, "kind": "type name", "text": v})
            if isinstance(node, ast.Compare) and any(isinstance(o, (ast.In, ast.NotIn)) for o in node.ops):
                left = node.left
                if isinstance(left, ast.Constant) and isinstance(left.value, str) and left.value in FRAGMENTS:
                    rows.append({"file": _rel(p), "line": node.lineno, "kind": "name fragment test",
                                 "text": ast.unparse(node)[:120]})
            if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "isinstance" and len(node.args) == 2:
                t = ast.unparse(node.args[1])
                if "Element" in t and not in_library:
                    rows.append({"file": _rel(p), "line": node.lineno, "kind": "isinstance on element class",
                                 "text": ast.unparse(node)[:120]})
    return rows


def probe_private():
    rows = []
    for p in _py(SRC):
        tree = _parse(p)
        defined = set()
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                defined.add(node.name)
            elif isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id in ("self", "cls"):
                defined.add(node.attr)
            elif isinstance(node, (ast.Name,)):
                defined.add(node.id)
            elif isinstance(node, ast.alias):
                defined.add(node.asname or node.name.split(".")[-1])
        for node in ast.walk(tree):
            if not isinstance(node, ast.Attribute):
                continue
            a = node.attr
            if not a.startswith("_") or a.startswith("__"):
                continue
            base = node.value
            if isinstance(base, ast.Name) and base.id in ("self", "cls"):
                continue
            if isinstance(base, ast.Call) and getattr(base.func, "id", None) == "super":
                continue
            if a in defined:
                continue
            rows.append({"file": _rel(p), "line": node.lineno, "attribute": a, "text": ast.unparse(node)[:100]})
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("onedfem"):
                for al in node.names:
                    if al.name.startswith("_") and not al.name.startswith("__"):
                        rows.append({"file": _rel(p), "line": node.lineno, "attribute": al.name,
                                     "text": f"from {node.module} import {al.name}"})
    return rows


def probe_domain():
    rows = []
    for p in _py(SRC):
        rel = _rel(p)
        if rel in ADAPTERS:
            continue
        src = p.read_text(encoding="utf-8-sig")
        tree = _parse(p)
        doc_lines = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and hasattr(node, "end_lineno"):
                doc_lines.update(range(node.lineno, node.end_lineno + 1))
        for i, line in enumerate(src.splitlines(), 1):
            for m in DOMAIN.finditer(line):
                where = "comment or docstring" if (i in doc_lines or line.strip().startswith("#")) else "code"
                rows.append({"file": rel, "line": i, "word": m.group(0).lower(), "where": where,
                             "text": line.strip()[:120]})
    return rows


def probe_outputs():
    rows = []
    pat = re.compile(r"""["']([A-Za-z0-9_]+\.(csv|txt|json))["']""")
    for p in _py(SRC):
        rel = _rel(p)
        if "pre_processing/parsing" in rel:
            continue
        src = p.read_text(encoding="utf-8-sig")
        if not re.search(r"\b(open\([^)]*['\"]w|to_csv|savetxt|csv\.writer|DictWriter|write_text)", src):
            continue
        files = sorted({m.group(1) for m in pat.finditer(src)})
        cols = []
        for m in re.finditer(r"([A-Z_]*COLUMNS)\s*=\s*(\([^)]*\)|\[[^\]]*\])", src, re.S):
            try:
                cols.extend(ast.literal_eval(m.group(2)))
            except (ValueError, SyntaxError):
                pass
        low = [str(c).lower() for c in cols]
        rows.append({"file": rel, "writes": " ".join(files), "n_columns": len(cols),
                     "says_where": any(w in c for c in low for w in WHERE) if cols else None,
                     "says_kind": any(k in c for c in low for k in KIND) if cols else None,
                     "columns": " ".join(map(str, cols))[:300]})
    return rows


def probe_plots():
    rows = []
    for base in (SRC, ROOT / "examples", ROOT / "validation"):
        if not base.exists():
            continue
        for p in _py(base):
            src = p.read_text(encoding="utf-8-sig")
            uses_style = "figure_style" in src or "resolution_plotting_utils" in src
            for i, line in enumerate(src.splitlines(), 1):
                if SMOOTHING.search(line):
                    rows.append({"file": _rel(p), "line": i, "issue": "smoothing or interpolation call",
                                 "text": line.strip()[:120]})
                if ".savefig(" in line and not uses_style:
                    rows.append({"file": _rel(p), "line": i, "issue": "figure saved without figure_style",
                                 "text": line.strip()[:120]})
    return rows


def probe_coverage():
    sys.path.insert(0, str(ROOT / "src"))
    from onedfem.pre_processing.element_library.element_factory import ElementFactory

    docs = {p: p.read_text(encoding="utf-8") for p in (ROOT / "docs" / "element_library").rglob("*.md")}
    tests = {p: p.read_text(encoding="utf-8") for p in (ROOT / "tests").glob("test_*.py")}
    rows = []
    for name in sorted(ElementFactory.ELEMENT_CLASS_MAP):
        d = [p.name for p, t in docs.items() if name in t]
        t = [p.name for p, s in tests.items() if name in s]
        rows.append({"element": name, "theory_pages": " ".join(d), "n_theory_pages": len(d),
                     "tests": " ".join(t), "n_tests": len(t)})
    return rows


PROBES = {"names": probe_names, "private": probe_private, "domain": probe_domain, "outputs": probe_outputs,
          "plots": probe_plots, "coverage": probe_coverage}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "screen_out"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name, fn in PROBES.items():
        rows = fn()
        summary[name] = len(rows)
        if rows:
            with open(out / f"probe_{name}.csv", "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
    (out / "static_probes.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
