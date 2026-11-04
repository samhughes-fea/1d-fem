"""v0.4.0 release gate checks (shim removal) — see docs/conventions/API_STANDARDS.md."""

from __future__ import annotations

import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _python_grep_legacy_modal_imports() -> list[str]:
    bad: list[str] = []
    for path in PROJECT_ROOT.rglob("*.py"):
        rel = path.relative_to(PROJECT_ROOT)
        s = str(rel).replace("\\", "/")
        if s.startswith("tests/") or s.startswith("docs/"):
            continue
        if "processing/modal/" in s:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for pat in (
            "from onedfem.processing.modal.assembly",
            "import onedfem.processing.modal.assembly",
            "from onedfem.processing.modal.boundary_conditions",
            "import onedfem.processing.modal.boundary_conditions",
            "from onedfem.processing.modal.buckling",
            "import onedfem.processing.modal.buckling",
        ):
            if pat in text:
                bad.append(f"{rel}: {pat}")
    return bad


def test_changelog_documents_v04_shim_removal_window() -> None:
    """Gate checklist: CHANGELOG tracks shim deprecation / removal (see API_STANDARDS.md)."""
    cl = (PROJECT_ROOT / "docs" / "CHANGELOG.md").read_text(encoding="utf-8")
    assert "## [Unreleased]" in cl
    assert re.search(
        r"v0\.4|Removed|processing\.modal",
        cl,
        re.IGNORECASE,
    ), "CHANGELOG should reference v0.4.0 gate and/or onedfem.processing.modal shim (API_STANDARDS.md)"


def test_no_library_imports_processing_modal_shim_paths() -> None:
    """Library code (excluding ``processing/modal`` shims and tests) must not import legacy paths."""
    bad = _python_grep_legacy_modal_imports()
    assert not bad, "Found legacy onedfem.processing.modal submodule imports:\n" + "\n".join(bad)


def test_simulation_runner_modal_shim_package_removed() -> None:
    """``onedfem.simulation_runner.modal`` re-export shims were removed; use ``onedfem.simulation_runner.spectral``."""
    legacy = PROJECT_ROOT / "src" / "onedfem" / "simulation_runner" / "modal"
    assert not legacy.is_dir(), (
        f"expected {legacy} to be removed; import from onedfem.simulation_runner.spectral instead"
    )

