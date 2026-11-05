"""The thin-walled section drawings take their values from the section's own solution."""

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from onedfem.post_processing.graphical_visualisers.thin_walled_drawing import (  # noqa: E402
    contour_field,
    draw_section,
    gauss_s,
    unrolled_field,
    wall_samples,
)
from anisotropic_section_model.thin_walled import (  # noqa: E402
    ThinWalledSection,
    ThinWalledSectionAnalysis,
    Wall,
    WallLaminate,
)
from anisotropic_section_model.thin_walled.laminate import Ply, PlyMaterial  # noqa: E402
from anisotropic_section_model.thin_walled.stress_recovery import SaintVenantRecovery  # noqa: E402


def _box():
    cfrp = PlyMaterial(135e9, 10e9, 5e9, 0.3, 5e9, 3.8e9, 1550.0)
    cap = WallLaminate([Ply(cfrp, a, 0.5e-3) for a in (45, -45, 0, 0, -45, 45)])
    web = WallLaminate([Ply(cfrp, a, 0.5e-3) for a in (45, -45)])
    pts = {1: (-0.3, -0.15), 2: (0.3, -0.15), 3: (0.3, 0.15), 4: (-0.3, 0.15)}
    walls = [Wall(1, 1, 2, 1), Wall(2, 2, 3, 2), Wall(3, 3, 4, 1), Wall(4, 4, 1, 2)]
    return ThinWalledSectionAnalysis(ThinWalledSection(pts, walls), {1: cap, 2: web}, n_sub=3), cfrp


def test_samples_and_gauss_points_lie_in_each_wall_element():
    an, _ = _box()
    L = an.section.length(an.section.wall(1))
    s = wall_samples(an, 1, per_element=5)
    assert len(s) == 3 and s[0][0] == 0.0 and s[-1][-1] == pytest.approx(L)
    g = gauss_s(an, 1)
    assert len(g) == 3 * 5 and np.all((g > 0) & (g < L))


def test_contour_field_is_the_recovered_field():
    import matplotlib.pyplot as plt

    an, cfrp = _box()
    rec = SaintVenantRecovery(an)
    F = [0.0, -1e5, 0.0, 0.0, 2e4, 5e3]
    fig, ax = plt.subplots()
    handles = draw_section(ax, an, [("CFRP", cfrp)])
    assert [h.get_label() for h in handles][0] == "CFRP"
    out = contour_field(ax, an, lambda w, s: rec.at(w, s, F).resultants[2], 1e-6)
    s, v = out[2]
    for si, vi in zip(s[::7], v[::7]):
        assert vi == pytest.approx(rec.at(2, float(si), F).resultants[2])
    starts = unrolled_field(ax, an, lambda w, s: rec.at(w, s, F).resultants[0], [1, 2, 3, 4])
    assert starts[-1] == pytest.approx(2 * 0.6 + 2 * 0.3)
    plt.close(fig)
