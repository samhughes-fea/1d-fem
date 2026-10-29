"""Thin-walled section geometry and its input parsers."""

import numpy as np
import pytest

from pre_processing.parsing.laminate_parser import parse_laminate_file
from pre_processing.parsing.section_geometry_parser import parse_section_geometry_file
from pre_processing.section_integration.thin_walled import ThinWalledSection, Wall, WallLaminate


def box(n_webs_inside=0):
    """Rectangle 0.4 x 0.2 with optional interior vertical webs."""
    xs = np.linspace(0.0, 0.4, n_webs_inside + 2)
    pts, walls = {}, []
    for i, x in enumerate(xs):
        pts[2 * i + 1] = (x, 0.0)
        pts[2 * i + 2] = (x, 0.2)
    wid = 1
    for i in range(len(xs) - 1):
        walls.append(Wall(wid, 2 * i + 1, 2 * i + 3, 1)); wid += 1
        walls.append(Wall(wid, 2 * i + 4, 2 * i + 2, 1)); wid += 1
    for i in range(len(xs)):
        walls.append(Wall(wid, 2 * i + 1, 2 * i + 2, 1)); wid += 1
    return ThinWalledSection(pts, walls)


def test_cell_counts():
    channel = ThinWalledSection(
        {1: (0.1, 0.0), 2: (0.0, 0.0), 3: (0.0, 0.2), 4: (0.1, 0.2)},
        [Wall(1, 1, 2, 1), Wall(2, 2, 3, 1), Wall(3, 3, 4, 1)],
    )
    assert channel.n_cells == 0 and not channel.is_closed
    assert box().n_cells == 1
    assert box(1).n_cells == 2
    assert box(2).n_cells == 3


def test_tangent_and_normal():
    sec = ThinWalledSection({1: (0.0, 0.0), 2: (3.0, 4.0)}, [Wall(7, 1, 2, 1)])
    w = sec.wall(7)
    assert sec.length(w) == pytest.approx(5.0)
    np.testing.assert_allclose(sec.tangent(w), [0.6, 0.8])
    np.testing.assert_allclose(sec.normal(w), [-0.8, 0.6])


@pytest.mark.parametrize(
    "points, walls, msg",
    [
        ({1: (0, 0), 2: (1, 0)}, [Wall(1, 1, 3, 1)], "unknown point"),
        ({1: (0, 0), 2: (0, 0)}, [Wall(1, 1, 2, 1)], "zero length"),
        ({1: (0, 0), 2: (1, 0)}, [Wall(1, 1, 1, 1)], "starts and ends"),
        ({1: (0, 0), 2: (1, 0)}, [Wall(1, 1, 2, 1), Wall(1, 2, 1, 1)], "duplicate wall_id"),
        ({1: (0, 0), 2: (1, 0)}, [Wall(1, 1, 2, 1), Wall(2, 2, 1, 1)], "more than one wall"),
        ({1: (0, 0), 2: (1, 0), 3: (2, 0)}, [Wall(1, 1, 2, 1)], "not on any wall"),
        ({1: (0, 0), 2: (1, 0)}, [], "at least one wall"),
    ],
)
def test_geometry_errors(points, walls, msg):
    with pytest.raises(ValueError, match=msg):
        ThinWalledSection(points, walls)


def test_unknown_laminate_rejected():
    with pytest.raises(ValueError, match="unknown laminate"):
        ThinWalledSection({1: (0, 0), 2: (1, 0)}, [Wall(1, 1, 2, 5)], laminate_ids=[1, 2])


LAMINATE_TXT = """# test laminates
[Ply_Material]
[material_id] [E1] [E2] [G12] [nu12] [G13] [G23] [rho]
1  140e9 10e9 5e9 0.3 5e9 3.5e9 1600
2  70e9  70e9 26.923076923e9 0.3 26.923076923e9 26.923076923e9 2700

[Laminate]
[laminate_id] [ply] [material_id] [angle_deg] [thickness]
1  2  1  90.0  0.5e-3   # top ply listed first on purpose
1  1  1   0.0  0.5e-3
2  1  2   0.0  2e-3
"""

GEOMETRY_TXT = """[Section_Point]
[section_id] [point_id] [y] [z]
1 1 0.0 0.0
1 2 0.4 0.0
1 3 0.4 0.2
1 4 0.0 0.2
2 1 0.0 0.0
2 2 0.0 0.1

[Section_Wall]
[section_id] [wall_id] [point_start] [point_end] [laminate_id]
1 1 1 2 1
1 2 2 3 2
1 3 3 4 1
1 4 4 1 2
2 1 1 2 2
"""


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_parse_laminate_file(tmp_path):
    lams = parse_laminate_file(write(tmp_path, "laminate.txt", LAMINATE_TXT))
    assert sorted(lams) == [1, 2]
    assert [p.angle_deg for p in lams[1].plies] == [0.0, 90.0]
    assert lams[1].thickness == pytest.approx(1e-3)
    assert lams[1].mass_per_area == pytest.approx(1.6)
    iso = WallLaminate.isotropic(70e9, 0.3, 2e-3)
    np.testing.assert_allclose(lams[2].ABD(), iso.ABD(), rtol=1e-8)


def test_parse_section_geometry_file(tmp_path):
    lams = parse_laminate_file(write(tmp_path, "laminate.txt", LAMINATE_TXT))
    secs = parse_section_geometry_file(write(tmp_path, "section_geometry.txt", GEOMETRY_TXT), lams.keys())
    assert sorted(secs) == [1, 2]
    assert secs[1].n_cells == 1
    assert secs[2].n_cells == 0
    assert secs[1].wall(2).laminate_id == 2


@pytest.mark.parametrize(
    "text, err",
    [
        (LAMINATE_TXT.replace("[G23] [rho]", "[G23] [density]"), ValueError),
        (LAMINATE_TXT.replace("1  1  1   0.0", "1  3  1   0.0"), ValueError),
        (LAMINATE_TXT.replace("2  1  2   0.0", "2  1  9   0.0"), ValueError),
        (LAMINATE_TXT.replace("140e9", "abc"), TypeError),
        (LAMINATE_TXT.split("[Laminate]")[0], ValueError),
        (LAMINATE_TXT.replace("2  70e9", "1  70e9"), ValueError),
    ],
)
def test_laminate_parser_errors(tmp_path, text, err):
    with pytest.raises(err):
        parse_laminate_file(write(tmp_path, "laminate.txt", text))


@pytest.mark.parametrize(
    "text",
    [
        GEOMETRY_TXT.replace("[point_end]", "[end]"),
        GEOMETRY_TXT.replace("1 4 0.0 0.2", "1 3 0.0 0.2"),
        GEOMETRY_TXT.replace("2 1 1 2 2", "3 1 1 2 2"),
        GEOMETRY_TXT.replace("1 4 4 1 2", "1 4 4 1 7"),
    ],
)
def test_geometry_parser_errors(tmp_path, text):
    with pytest.raises(ValueError):
        parse_section_geometry_file(write(tmp_path, "section_geometry.txt", text), [1, 2])
