"""1d-fem's section job input (``pre_processing/section_job.py``) and its section provider seam."""

import pytest

from onedfem.sections import (
    DEFAULT_SECTION_PROVIDER,
    SectionProviderError,
    check_section_provider,
    section_import,
    section_provider,
    set_section_provider,
)

E, NU = 70e9, 0.3
G = E / (2 * (1 + NU))


def test_default_provider_has_everything_1d_fem_uses():
    assert section_provider().__name__ == DEFAULT_SECTION_PROVIDER
    assert check_section_provider() == []


def test_a_missing_provider_is_reported():
    set_section_provider("no_such_section_package")
    try:
        with pytest.raises(SectionProviderError, match="no_such_section_package"):
            section_import("thin_walled.laminate", "WallLaminate")
    finally:
        set_section_provider(None)
    assert section_provider().__name__ == DEFAULT_SECTION_PROVIDER


def test_another_provider_plugs_in():
    import types

    import anisotropic_section_model

    wrapper = types.ModuleType("wrapped_sections")
    wrapper.__name__ = anisotropic_section_model.__name__  # same modules, different object
    set_section_provider(wrapper)
    try:
        assert section_provider() is wrapper
        assert section_import("ply_failure", "CRITERIA") == ("max_stress", "tsai_wu", "hashin", "puck")
    finally:
        set_section_provider(None)


def test_section_analyses_from_files(tmp_path):
    from onedfem.pre_processing.section_job import section_analyses_from_files

    (tmp_path / "laminate.txt").write_text(
        "[Ply_Material]\n[material_id] [E1] [E2] [G12] [nu12] [G13] [G23]\n"
        f"1 {E} {E} {G} {NU} {G} {G}\n"
        "[Laminate]\n[laminate_id] [ply] [material_id] [angle_deg] [thickness]\n1 1 1 0.0 0.002\n",
        encoding="utf-8",
    )
    (tmp_path / "section_geometry.txt").write_text(
        "[Section_Point]\n[section_id] [point_id] [y] [z]\n"
        "1 1 -0.2 -0.1\n1 2 0.2 -0.1\n1 3 0.2 0.1\n1 4 -0.2 0.1\n"
        "[Section_Wall]\n[section_id] [wall_id] [point_start] [point_end] [laminate_id]\n"
        "1 1 1 2 1\n1 2 2 3 1\n1 3 3 4 1\n1 4 4 1 1\n",
        encoding="utf-8",
    )
    out = section_analyses_from_files(str(tmp_path / "laminate.txt"), str(tmp_path / "section_geometry.txt"), hoop="free")
    K = out[1].stiffness
    assert K[0, 0] == pytest.approx(E * 1.2 * 0.002, rel=1e-10)
