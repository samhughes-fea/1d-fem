"""The section-provider contract: versions, and when the provider is imported."""

from __future__ import annotations

import types
import warnings

import pytest

from onedfem import sections
from onedfem.sections import (
    SECTION_API,
    SectionProviderError,
    check_api_version,
    section_provider,
    set_section_provider,
)


def _module(name, version=None):
    m = types.ModuleType(name)
    if version is not None:
        m.SECTION_API_VERSION = version
    return m


def test_the_default_provider_meets_the_contract():
    assert check_api_version(section_provider())[0] == SECTION_API[0]


@pytest.mark.parametrize("version", [(SECTION_API[0] + 1, 0), (SECTION_API[0] - 1, 9)])
def test_another_major_version_is_refused(version):
    with pytest.raises(SectionProviderError, match="section contract"):
        check_api_version(_module("future_sections", version))


def test_an_older_minor_version_is_refused(monkeypatch):
    monkeypatch.setattr(sections, "SECTION_API", (1, 3))
    with pytest.raises(SectionProviderError, match="x >= 3"):
        check_api_version(_module("old_sections", (1, 2)))
    assert check_api_version(_module("new_sections", (1, 4))) == (1, 4)


def test_a_provider_without_a_version_is_taken_as_1_0_with_a_warning():
    m = _module("unversioned_sections")
    with pytest.warns(UserWarning, match="SECTION_API_VERSION"):
        assert check_api_version(m) == (1, 0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # once only
        check_api_version(m)


def test_a_provider_module_is_checked_when_set():
    set_section_provider(_module("future_sections", (99, 0)))
    try:
        with pytest.raises(SectionProviderError, match="section contract"):
            section_provider()
    finally:
        set_section_provider(None)


def test_attaching_no_section_files_needs_no_provider(tmp_path):
    """With no section files, attaching sections imports no provider. (The isotropic elements
    themselves take their ``D`` from it: ``test_isotropic_sections_from_provider.py``.)"""
    from onedfem.pre_processing.section_job import attach_section_analyses

    job = tmp_path / "job"
    job.mkdir()
    ed = {"ids": [0, 1]}
    set_section_provider("no_such_section_package")
    try:
        assert attach_section_analyses(str(job), ed) is ed
        (job / "element_section.txt").write_text("[Element_Section]\n[element_id] [section_id]\n0 1\n")
        with pytest.raises(SectionProviderError, match="no_such_section_package"):
            attach_section_analyses(str(job), ed)
    finally:
        set_section_provider(None)


def test_isotropic_elements_need_the_provider(tmp_path):
    """An isotropic section.txt element gets its section from the provider, so a missing one fails."""
    import numpy as np

    from onedfem.pre_processing.element_library.isotropic_section import isotropic_section

    set_section_provider("no_such_section_package")
    try:
        with pytest.raises(SectionProviderError, match="no_such_section_package"):
            isotropic_section(E=1.0, G=1.0, A=1.0, I_y=1.0, I_z=1.0, J_t=1.0)
    finally:
        set_section_provider(None)
    assert np.isfinite(isotropic_section(E=1.0, G=1.0, A=1.0, I_y=1.0, I_z=1.0, J_t=1.0).stiffness).all()
