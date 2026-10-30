"""Distributed and wall loads on rotated, curved and pretwisted beam elements.

References:

- A straight beam turned to any direction under a line load per length gives the
  unturned beam's response, turned: the load is spread along the element's own axis.
- A quarter-circle cantilever under a uniform out-of-plane line load: four curved
  elements against 64 straight ones.
- The same quarter circle under a wall load (a section load, which turns with the
  section): four curved elements against 64 straight ones.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.test_element_frames import BS, EB, TS, Mesh, _line
from tests.test_second_order import channel

L = 2.0


def _rotation_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _table(xs, q6):
    """``distributed_load.txt`` rows at global ``X`` positions ``xs`` with loads ``q6(i)``."""
    return np.array([[x, 0.0, 0.0, *q6(i)] for i, x in enumerate(xs)])


def _turn_nodal(U, R, dpn):
    """Nodal translations and rotations turned by ``R``."""
    V = U.reshape(-1, dpn).copy()
    V[:, 0:3] = V[:, 0:3] @ R.T
    V[:, 3:6] = V[:, 3:6] @ R.T
    return V


@pytest.mark.parametrize("etype", [EB, TS, BS])
def test_turned_beam_under_a_varying_line_load(etype):
    """Triangular load across the beam on a cantilever turned 30 degrees about Z."""
    n, a = 6, np.radians(30.0)
    R = _rotation_z(a)
    s = np.linspace(0.0, L, 13)
    if etype == EB:   # the reference kernel integrates line loads with two points: exact for uniform ones
        q_ref = lambda i: np.array([0.0, 300.0, 500.0, 20.0, 0.0, 0.0])
    else:
        q_ref = lambda i: np.array([0.0, 300.0 * s[i] / L, 500.0 * (1 - s[i] / L), 20.0, 0.0, 0.0])
    sa = channel() if etype == BS else None
    ref = Mesh(etype, _line([0, 0, 0], [1, 0, 0], L, n), sa=sa, dl=_table(s, q_ref))
    # turned: samples at X = s cos a carry the reference loads turned by R
    turned_q = lambda i: np.concatenate([R @ q_ref(i)[:3], R @ q_ref(i)[3:]])
    tur = Mesh(etype, _line([0, 0, 0], R[:, 0], L, n), sa=sa, dl=_table(s * np.cos(a), turned_q))
    fixed = list(range(ref.dpn))
    U0 = ref.solve(ref.F, fixed)
    U1 = tur.solve(tur.F, fixed)
    V0 = _turn_nodal(U0, R, ref.dpn)
    V1 = U1.reshape(-1, tur.dpn)
    assert np.abs(V0[:, :6]).max() > 0
    np.testing.assert_allclose(V1[:, :6], V0[:, :6], rtol=1e-9, atol=1e-12 * np.abs(V0).max())


def test_quarter_circle_under_a_uniform_line_load():
    """Four curved elements against 64 straight ones, out-of-plane uniform load."""
    Rc, q = 1.0, 100.0

    def mesh(n, curved):
        th = np.linspace(0.0, np.pi / 2, n + 1)
        pts = np.column_stack([Rc * np.sin(th), Rc * (1 - np.cos(th)), np.zeros(n + 1)])
        dl = _table(np.linspace(0.0, Rc, 5), lambda i: [0, 0, -q, 0, 0, 0])
        return Mesh(TS, pts, curvature=[0, 0, 1 / Rc] if curved else None, dl=dl)

    fine, arc = mesh(64, False), mesh(4, True)
    tip_f = fine.solve(fine.F, list(range(6)))[fine.dof(64, 2)]
    tip_a = arc.solve(arc.F, list(range(6)))[arc.dof(4, 2)]
    # total load on the arc: q pi R / 2 on both
    assert abs(fine.F[2::fine.dpn].sum()) == pytest.approx(q * np.pi * Rc / 2, rel=1e-3)
    assert abs(arc.F[2::arc.dpn].sum()) == pytest.approx(q * np.pi * Rc / 2, rel=1e-6)
    assert tip_a == pytest.approx(tip_f, rel=3e-3)


@pytest.mark.parametrize("etype", [TS, BS])
def test_turned_beam_under_a_wall_load(etype):
    """A section load on a turned element turns with it."""
    sa = channel()
    wall = np.array([0.0, 40.0, -120.0, 3.0, 0.0, 0.0, 0.0])
    n, R = 4, _rotation_z(np.radians(50.0))
    ref = Mesh(etype, _line([0, 0, 0], [1, 0, 0], L, n), sa=sa, wall=wall)
    tur = Mesh(etype, _line([0, 0, 0], R[:, 0], L, n), sa=sa, wall=wall)
    fixed = list(range(ref.dpn))
    V0 = _turn_nodal(ref.solve(ref.F, fixed), R, ref.dpn)
    V1 = tur.solve(tur.F, fixed).reshape(-1, tur.dpn)
    np.testing.assert_allclose(V1, V0, rtol=1e-9, atol=1e-12 * np.abs(V0).max())


def test_quarter_circle_under_a_wall_load():
    """A section load on curved elements: four curved elements against 64 straight ones."""
    sa = channel()
    wall = np.array([0.0, 30.0, -100.0, 2.0, 0.0, 0.0, 0.0])
    Rc = 1.0

    def mesh(n, curved):
        th = np.linspace(0.0, np.pi / 2, n + 1)
        pts = np.column_stack([Rc * np.sin(th), Rc * (1 - np.cos(th)), np.zeros(n + 1)])
        return Mesh(TS, pts, sa=sa, curvature=[0, 0, 1 / Rc] if curved else None, wall=wall)

    fine, arc = mesh(64, False), mesh(4, True)
    tip_f = fine.solve(fine.F, list(range(6))).reshape(-1, 6)[-1]
    tip_a = arc.solve(arc.F, list(range(6))).reshape(-1, 6)[-1]
    np.testing.assert_allclose(tip_a, tip_f, rtol=5e-3, atol=5e-3 * np.abs(tip_f).max())


def test_line_load_outside_the_samples_is_refused():
    n, R = 2, _rotation_z(np.radians(30.0))
    with pytest.raises(ValueError, match="distributed load samples cover"):
        Mesh(TS, _line([0, 0, 0], R[:, 0], L, n), dl=_table([0.0, 0.5], lambda i: [0, 0, 1.0, 0, 0, 0]))
