"""Corotational Timoshenko beam against the CalculiX ``NLGEOM`` solid composite box.

The reference (``post_processing/validation_visualisers/nonlinear_solid_box/reference.json``)
is a ``C3D20R`` model of the CAS CFRP box with mitred corners, clamped at the root, under a
dead flapwise tip force up to 3 P (P gives 20 % of the span in linear theory).
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
for p in (REPO, os.path.join(REPO, "post_processing", "validation_visualisers", "iea15mw_sections")):
    if p not in sys.path:
        sys.path.insert(0, p)

REFERENCE = os.path.join(REPO, "post_processing", "validation_visualisers", "nonlinear_solid_box", "reference.json")


@pytest.fixture(scope="module")
def ref():
    with open(REFERENCE, encoding="utf-8") as f:
        return json.load(f)


def _solid_tip(ref, level):
    nl = ref["nlgeom"]
    return np.array(nl["tip"][nl["levels"].index(level)])


def test_linear_solid_matches_the_beam(ref):
    """The solid's linear tip deflection is the beam's: the section and the solid agree."""
    from post_processing.validation_visualisers.abaqus.shell_validation.cases import BOX_L

    assert ref["linear"]["tip"][0][2] == pytest.approx(0.2 * BOX_L, rel=2e-3)


def test_corotational_beam_follows_the_solid(tmp_path, ref):
    """At 2 P (35 % of the span) the tip deflection is within 1.3 %, the twist within 2 %,
    and the tip shortening within 1 % of the deflection."""
    from post_processing.validation_visualisers.nonlinear_solid_box.run import run_job

    beam = run_job(str(tmp_path), "thin_walled", 2.0)[-1]
    solid = _solid_tip(ref, 2.0)
    assert beam[2] == pytest.approx(solid[2], rel=0.015)
    assert beam[3] == pytest.approx(solid[3], rel=0.025)
    assert abs(beam[0] - solid[0]) < 0.01 * solid[2]


def test_solid_flattens_like_brazier(ref):
    """The solid's flanges close in under bending (Brazier): about quadratically in the load,
    3 % of the depth at a quarter span at 3 P, which is why the beam is 5 % stiff there."""
    nl = ref["nlgeom"]
    f = np.array(nl["flattening"]["0.25"])
    lv = np.array(nl["levels"])
    assert f[-1] == pytest.approx(0.0033, rel=0.1)
    assert np.all(np.diff(f) > 0)
    assert f[lv.tolist().index(2.0)] / f[lv.tolist().index(1.0)] > 3.0


@pytest.mark.slow
def test_gbt_beam_flattens_like_the_solid(tmp_path, ref):
    """The corotational GBT beam's Brazier term flattens the flanges: at P its quarter-span
    flattening is 86 % of the solid's, and its tip deflection is within 0.2 %."""
    from post_processing.validation_visualisers.nonlinear_solid_box.run import (GBT_MODES, N_EL, flattening,
                                                                                gbt_modes, run_gbt_job)

    modes = gbt_modes()
    U = run_gbt_job(str(tmp_path), 1.0, modes)
    k = ref["nlgeom"]["levels"].index(1.0)
    assert U[-1, 2] == pytest.approx(_solid_tip(ref, 1.0)[2], rel=2e-3)
    a = U[N_EL // 4, 6 + modes.n_sv_warping::2][:GBT_MODES]
    assert flattening(modes, a) / ref["nlgeom"]["flattening"]["0.25"][k] == pytest.approx(0.86, abs=0.03)
