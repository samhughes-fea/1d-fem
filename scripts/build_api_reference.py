"""Build the API reference in ``docs/api/`` from the modules' own docstrings.

Run from the repository root::

    python scripts/build_api_reference.py           # write docs/api/
    python scripts/build_api_reference.py --check   # exit 1 if docs/api/ is out of date

A module is documented when it declares ``__all__`` or has a ``Public API`` section in its
docstring. Each page lists, per module, the module summary (the first paragraph of its
docstring) and each public name with its signature and the first paragraph of its docstring.
The element classes are listed from ``ElementFactory.ELEMENT_CLASS_MAP``.

Modules that cannot be imported (an optional dependency missing) are listed by name with the
import error, so a page never depends on which optional packages are installed.
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import os
import re
import sys
import warnings
from typing import Dict, List, Optional, Tuple

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SRC = os.path.join(REPO, "src")
OUT = os.path.join(REPO, "docs", "api")

#: page file -> (title, package roots, intro); a root "+pkg" lists every module under pkg
PAGES: Dict[str, Tuple[str, Tuple[str, ...], str]] = {
    "sections.md": ("Section analysis", ("onedfem.pre_processing.section_integration",),
                    "Thin-walled and 2D (solid) cross-section analysis: stiffness, mass, warping, "
                    "GBT modes, stress recovery, ply failure and layup sensitivities."),
    "elements.md": ("Elements", ("onedfem.pre_processing.element_library",),
                    "Beam, bar and truss elements and their shared utilities."),
    "job_input.md": ("Job input", ("+onedfem.pre_processing.parsing", "onedfem.pre_processing.mesh_library",
                                   "onedfem.pre_processing.load_library", "onedfem.pre_processing.containers"),
                     "Parsers for the job directory's text files, meshes and loads."),
    "solvers.md": ("Solvers", ("onedfem.processing", "onedfem.simulation_runner"),
                   "Static, nonlinear, eigen, buckling, harmonic, transient and spectral solvers and "
                   "their runners."),
    "post_processing.md": ("Post-processing", ("onedfem.post_processing.failure", "onedfem.post_processing.sensitivity",
                                               "onedfem.post_processing.solid", "onedfem.post_processing.graphical_visualisers"),
                           "Result files written after a solve: stresses, failure margins, "
                           "sensitivities and field sampling for plots."),
    "workflow.md": ("Running jobs", ("+onedfem.workflow_orchestrator.run_job",),
                    "Running a job directory end to end."),
}
SKIP_DIRS = {"__pycache__", "results", "generated", "tests"}


def _modules(root: str) -> List[str]:
    base = os.path.join(SRC, *root.split("."))
    if os.path.isfile(base + ".py"):
        return [root]
    out = []
    for dirpath, dirnames, filenames in os.walk(base):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
        rel = os.path.relpath(dirpath, SRC).replace(os.sep, ".")
        for f in sorted(filenames):
            if f.endswith(".py"):
                out.append(rel if f == "__init__.py" else f"{rel}.{f[:-3]}")
    return sorted(set(out))


def _path(mod_name: str) -> str:
    base = os.path.join(SRC, *mod_name.split("."))
    return os.path.join(base, "__init__.py") if os.path.isdir(base) else base + ".py"


def _declares_api(path: str) -> bool:
    try:
        with open(path, encoding="utf-8") as f:
            src = f.read()
    except OSError:
        return False
    return bool(re.search(r"^__all__\s*=", src, re.M)) or "Public API" in src


def _first_paragraph(doc: Optional[str]) -> str:
    if not doc:
        return ""
    doc = inspect.cleandoc(doc)
    para = doc.split("\n\n")[0]
    return " ".join(line.strip() for line in para.splitlines())


def _signature(obj) -> str:
    try:
        sig = str(inspect.signature(obj))
    except (TypeError, ValueError):
        return ""
    sig = re.sub(r" at 0x[0-9a-f]+", "", sig)
    sig = re.sub(r"(: |-> )'([^']*)'", r"\1\2", sig)     # string annotations
    return re.sub(r"<([A-Za-z_.]+) object>", r"\1()", sig)


def _public_names(mod) -> List[str]:
    names = getattr(mod, "__all__", None)
    if names is not None:
        return [n for n in names if hasattr(mod, n)]
    out = []
    for n, v in vars(mod).items():
        if n.startswith("_"):
            continue
        if (inspect.isclass(v) or inspect.isfunction(v)) and getattr(v, "__module__", None) == mod.__name__:
            out.append(n)
    return out


def _describe(mod_name: str) -> str:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mod = importlib.import_module(mod_name)
    except Exception as e:  # noqa: BLE001 - optional dependencies
        return f"### `{mod_name}`\n\nNot imported here ({type(e).__name__}: {str(e).splitlines()[0][:120]}).\n"
    path = os.path.relpath(getattr(mod, "__file__", "") or "", REPO).replace(os.sep, "/")
    lines = [f"### `{mod_name}`", "", f"[{path}](../../{path})", ""]
    summary = _first_paragraph(mod.__doc__)
    if summary:
        lines += [summary, ""]
    for n in _public_names(mod):
        obj = getattr(mod, n)
        if inspect.isclass(obj):
            lines.append(f"- **class `{n}{_signature(obj)}`**")
        elif callable(obj):
            lines.append(f"- **`{n}{_signature(obj)}`**")
        else:
            lines.append(f"- **`{n}`** (constant)")
        d = _first_paragraph(obj.__doc__) if (inspect.isclass(obj) or inspect.isroutine(obj)) else ""
        if d:
            lines.append(f"  {d}")
    lines.append("")
    return "\n".join(lines)


def _element_table() -> str:
    from onedfem.pre_processing.element_library.element_factory import ElementFactory

    rows = ["| Element type | Module | Summary |", "|---|---|---|"]
    for etype, mod_name in sorted(ElementFactory.ELEMENT_CLASS_MAP.items()):
        try:
            cls = getattr(importlib.import_module(mod_name), etype)
            s = _first_paragraph(cls.__doc__).replace("|", "\\|")
        except Exception as e:  # noqa: BLE001
            s = f"not imported here ({type(e).__name__})"
        rows.append(f"| `{etype}` | `{mod_name}` | {s} |")
    return "\n".join(rows) + "\n"


def build() -> Dict[str, str]:
    for p in (REPO, SRC):
        if p not in sys.path:
            sys.path.insert(0, p)
    import logging

    logging.disable(logging.CRITICAL)
    pages: Dict[str, str] = {}
    index = ["# API reference", "",
             "Generated from the modules' docstrings by `python scripts/build_api_reference.py`; "
             "do not edit by hand. A module is listed when it declares `__all__` or a `Public API` "
             "section. The theory behind each part is in [element_library](../element_library/) and the "
             "job file formats in [conventions](../conventions/).", "",
             "| Page | Covers |", "|---|---|"]
    for fname, (title, roots, intro) in PAGES.items():
        body = [f"# {title}", "", intro, ""]
        if fname == "elements.md":
            body += ["## Element types", "", "The `element_type` names a job's `element.txt` accepts.", "",
                     _element_table(), "## Modules", ""]
        mods = []
        for r in roots:                         # "+root": every module under it, not only declared APIs
            every = r.startswith("+")
            mods += [m for m in _modules(r.lstrip("+")) if every or _declares_api(_path(m))]
        for m in mods:
            body.append(_describe(m))
        pages[fname] = "\n".join(body).rstrip() + "\n"
        index.append(f"| [{title}]({fname}) | {intro} |")
    pages["README.md"] = "\n".join(index) + "\n"
    logging.disable(logging.NOTSET)
    return pages


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="exit 1 if docs/api/ is out of date")
    args = ap.parse_args(argv)
    pages = build()
    if args.check:
        stale = []
        for f, text in pages.items():
            p = os.path.join(OUT, f)
            if not os.path.exists(p) or open(p, encoding="utf-8").read() != text:
                stale.append(f)
        if stale:
            print("docs/api is out of date: " + ", ".join(stale) + ". Run python scripts/build_api_reference.py")
            return 1
        return 0
    os.makedirs(OUT, exist_ok=True)
    for f, text in pages.items():
        with open(os.path.join(OUT, f), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    print(f"wrote {len(pages)} pages to {os.path.relpath(OUT, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
