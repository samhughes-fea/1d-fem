"""Smoke tests of examples/blade_section_2d and examples/blade_beam."""

import importlib.util
import os

import pytest

pytest.importorskip("implicit_geometry_kernel")

ROOT = os.path.join(os.path.dirname(__file__), "..")


def _skip_without_reference_blades():
    """The blade examples read the IEA 15 MW stations from the private reference-blades checkout."""
    spec = importlib.util.spec_from_file_location(
        "reference_blades_locator", os.path.join(ROOT, "examples", "iea15mw", "reference_blades.py"))
    blades = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(blades)
    try:
        blades.reference_blade_file("iea15/derived/sections/precomp/precomp_sections.json")
    except blades.ReferenceBladesUnavailable as exc:
        pytest.skip(str(exc), allow_module_level=True)


_skip_without_reference_blades()


def _load(rel, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, *rel.split("/")))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_blade_section_2d_example(tmp_path):
    mod = _load("examples/blade_section_2d/run.py", "blade_section_2d_run")
    s = mod.main(plot=False, out_dir=str(tmp_path))
    # caps: the 2D section holds more carbon than the thin-walled caps, so a little less stress
    assert -0.08 < s["cap_sigma_xx"]["difference"] < -0.02
    w = s["max_interface_shear"]
    assert w is not None and 1e6 < abs(w["tau_xn"]) < 5e6
    assert any(lab.startswith("spar_cap") for lab in w["labels"])
    assert s["panel"]["n_elements"] == 8          # four plies, two rows each
    assert os.path.exists(tmp_path / "summary.json")


def test_blade_beam_example_static(tmp_path):
    mod = _load("examples/blade_beam/run.py", "blade_beam_run")
    r = mod.main(n_elements=4, plot=False, work_dir=str(tmp_path), parts=("static",), out_dir=str(tmp_path))
    ref = mod.reference()["static"]["tip"]
    st = r["static"]
    # edgewise bending and Benscoter torsion within 3 % of the shell blade
    assert st[mod.TW]["edge"][0] == pytest.approx(ref["edge"][0], rel=0.03)
    assert st[mod.BENSCOTER]["torque"][5] == pytest.approx(ref["torque"][5], rel=0.03)
    # flap: the beam is softer on the tapered blade (the shell's inclined caps carry shear)
    assert 0.0 < st[mod.TW]["flap"][1] / ref["flap"][1] - 1.0 < 0.08
