# interfaces/blade/solid_sections.py
"""Conforming 2D section meshes of an implicit-geometry-kernel blade section.

``section_mesh_from_cut`` meshes a ``SectionCut`` and ``write_section_mesh_2d`` writes
``section_mesh_2d.txt``. Both are implemented by the section provider
(``solid.from_geometry``). This module is the 1d-fem import path and re-exports them.
"""

from __future__ import annotations

from typing import Any

__all__ = ["section_mesh_from_cut", "write_section_mesh_2d"]

_LOADED = False


def _load() -> None:
    global _LOADED
    if _LOADED:
        return
    from onedfem.sections import section_import

    mesh = section_import("solid.from_geometry", "section_mesh_from_cut", "write_section_mesh_2d")
    names = globals()
    names["section_mesh_from_cut"], names["write_section_mesh_2d"] = mesh
    _LOADED = True


def __getattr__(name: str) -> Any:
    if name in __all__:
        _load()
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list:
    return sorted(set(globals()) | set(__all__))
