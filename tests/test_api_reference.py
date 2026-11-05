"""The API reference generator runs and documents every module it lists."""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import build_api_reference as api  # noqa: E402


def test_every_listed_module_imports_and_every_page_has_content():
    pages = api.build()
    assert set(pages) == set(api.PAGES) | {"README.md"}
    for name, text in pages.items():
        if name == "README.md":
            continue
        assert "### `" in text, name
        assert "Not imported here" not in text or "implicit_geometry_kernel" in text, name
    assert "`LinearTimoshenkoBeamElement3D`" in pages["elements.md"]
    assert "process_job(" in pages["workflow.md"]
    assert "set_section_provider(" in pages["sections.md"]


def test_committed_reference_lists_the_same_modules():
    """Signatures may print differently between Python versions, so only the module list is compared."""
    for name, text in api.build().items():
        path = os.path.join(api.OUT, name)
        assert os.path.exists(path), f"docs/api/{name} missing: run python scripts/build_api_reference.py"
        with open(path, encoding="utf-8") as f:
            committed = f.read()
        heads = [ln for ln in text.splitlines() if ln.startswith("### `")]
        assert heads == [ln for ln in committed.splitlines() if ln.startswith("### `")], (
            f"docs/api/{name} is out of date: run python scripts/build_api_reference.py")
