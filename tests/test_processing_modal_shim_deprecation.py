"""``processing.modal.assembly`` / ``boundary_conditions`` / ``buckling`` shims removed — import canonical modules."""

from __future__ import annotations

import importlib
import sys

import pytest


@pytest.mark.parametrize(
    "mod",
    (
        "processing.modal.assembly",
        "processing.modal.boundary_conditions",
        "processing.modal.buckling",
    ),
)
def test_processing_modal_submodule_shims_removed(mod: str) -> None:
    sys.modules.pop(mod, None)
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module(mod)


def test_processing_modal_package_doc_stub_importable() -> None:
    import processing.modal  # noqa: F401 — package placeholder remains
