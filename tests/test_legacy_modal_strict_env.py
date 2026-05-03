"""Opt-in hard failure for legacy ``[Modal]`` / ``[Type] modal`` (``FEM_LEGACY_MODAL_ERROR``)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings


def test_legacy_modal_fixture_raises_when_strict_env(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("FEM_LEGACY_MODAL_ERROR", "1")
    fixture = (
        Path(__file__).resolve().parent.parent
        / "jobs"
        / "fixtures"
        / "simulation_settings_legacy_modal_vibration.txt"
    )
    assert fixture.is_file()
    with pytest.raises(ValueError, match="Legacy \\[Modal\\]"):
        parse_simulation_settings(str(fixture))


def test_legacy_modal_fixture_parses_when_strict_env_off(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("FEM_LEGACY_MODAL_ERROR", raising=False)
    fixture = (
        Path(__file__).resolve().parent.parent
        / "jobs"
        / "fixtures"
        / "simulation_settings_legacy_modal_vibration.txt"
    )
    s = parse_simulation_settings(str(fixture))
    assert s["type"] == "eigen"
