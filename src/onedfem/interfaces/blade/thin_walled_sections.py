# interfaces/blade/thin_walled_sections.py
"""Thin-walled sections cut from an implicit-geometry-kernel blade.

``thin_walled_from_cut`` turns a ``SectionCut`` (``BladeGeometry.section_at(z)``) into a
``ThinWalledSection`` and its laminates. ``section_taper_from_geometry`` is the taper
(Resal) coupling at a span station. Both are implemented by the section provider
(``thin_walled.from_geometry`` and ``thin_walled.taper``). This module is the 1d-fem
import path and re-exports them, so the core does not keep a second copy.
"""

from __future__ import annotations

from typing import Any

__all__ = ["GeometrySection", "thin_walled_from_cut", "laminate_lookup", "section_taper_from_geometry"]

_LOADED = False


def _load() -> None:
    global _LOADED
    if _LOADED:
        return
    from onedfem.sections import section_import

    geometry = section_import(
        "thin_walled.from_geometry", "GeometrySection", "thin_walled_from_cut", "laminate_lookup")
    taper = section_import("thin_walled.taper", "section_taper_from_geometry")
    names = globals()
    names["GeometrySection"], names["thin_walled_from_cut"], names["laminate_lookup"] = geometry
    names["section_taper_from_geometry"] = taper
    _LOADED = True


def __getattr__(name: str) -> Any:
    if name in __all__:
        _load()
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list:
    return sorted(set(globals()) | set(__all__))
