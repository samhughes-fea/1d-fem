"""Corotational beams that start rotated, twisted, kinked and pretwisted.

References:

- A corotational cantilever turned to any direction and twisted about its axis follows
  the unturned beam's path, turned, to large deflection.
- An initially curved cantilever (a quarter circle of radius ``R`` built from straight
  chords, stress-free) under a tip moment ``M`` keeps a uniform curvature ``1/R + M/EI``
  (pure bending of an inextensible rod): with ``M = EI / R`` it closes into a semicircle
  of radius ``R / 2``, with its tip at ``(0, R)`` and turned by ``pi``.
- A prebent cantilever's large-deflection path tends to the straight one as the prebend
  goes to zero, with a difference of the order of the prebend.
- At small load a corotational element is its linear core, here on rotated, kinked and
  pretwisted meshes.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from pre_processing.element_library.element_frame import straight_frame
from tests.test_element_frames import Mesh, _line
from tests.test_second_order import E, channel

CR_EB, CR_TS, CR_BS = (
    "CorotationalEulerBernoulliBeamElement3D",
    "CorotationalTimoshenkoBeamElement3D",
    "CorotationalBenscoterBeamElement3D",
)
LIN = {CR_EB: "LinearEulerBernoulliBeamElement3D", CR_TS: "LinearTimoshenkoBeamElement3D",
       CR_BS: "LinearBenscoterBeamElement3D"}
B, H = 0.05, 0.02            # rectangle: bending about local y is the soft direction
RECT = dict(A=B * H, I_x=B * H**3 / 12 + H * B**3 / 12, I_y=B * H**3 / 12, I_z=H * B**3 / 12,
            J_t=0.229 * B * H**3, kappa=5 / 6)


_JOB = tempfile.mkdtemp(prefix="cr_frames_")
for _sub in ("element_stiffness_matrices", "element_force_vectors", "logs"):
    os.makedirs(os.path.join(_JOB, _sub), exist_ok=True)


def _mesh(etype, coords, **kw):
    sa = channel() if "Benscoter" in etype else None
    mesh = Mesh(etype, coords, sa=sa, props=None if sa is not None else RECT, **kw)
    for el in mesh.elements:      # the corotational set-up builds its core lazily, and logs
        el.job_results_dir = _JOB
    return mesh


def _assemble(mesh, U):
    Fi = np.zeros(mesh.n)
    Kt = np.zeros((mesh.n, mesh.n))
    for el, g in zip(mesh.elements, mesh.maps):
        Fi[g] += el.internal_force_vector(U[g])
        Kt[np.ix_(g, g)] += el.tangent_stiffness_matrix(U[g])
    return Fi, Kt


def _newton(mesh, F, fixed, steps, tol=1e-8, max_it=40):
    free = np.setdiff1d(np.arange(mesh.n), fixed)
    U = np.zeros(mesh.n)
    for lam in np.linspace(1.0 / steps, 1.0, steps):
        for _ in range(max_it):
            Fi, Kt = _assemble(mesh, U)
            r = (lam * F - Fi)[free]
            if np.linalg.norm(r) <= tol * np.linalg.norm(F):
                break
            U[free] += np.linalg.solve(Kt[np.ix_(free, free)], r)
        else:
            raise AssertionError(f"Newton did not converge at load factor {lam}")
    return U


def _tip(mesh, node, f6):
    F = np.zeros(mesh.n)
    F[mesh.dof(node, 0):mesh.dof(node, 6)] = f6
    return F


@pytest.mark.parametrize("etype", [CR_EB, CR_TS, CR_BS])
def test_turned_corotational_cantilever_follows_the_turned_path(etype):
    """Tip load bending the tip by about L/3, on a cantilever turned off x and twisted."""
    L, n, twist = 1.0, 8, 35.0
    d = np.array([1.0, 2.0, -2.0]) / 3.0
    ref = _mesh(etype, _line([0, 0, 0], [1, 0, 0], L, n), twist=twist)
    tur = _mesh(etype, _line([0.3, -0.2, 0.1], d, L, n), twist=twist)
    R0 = straight_frame([0, 0, 0], [1, 0, 0], twist)
    R1 = straight_frame([0, 0, 0], d, twist)
    Q = R1 @ R0.T                          # takes the reference beam onto the turned one
    ref.elements[0]._cr_setup()
    EI = ref.elements[0]._K_core[2, 2] * L**3 / (12 * n**3)    # 12 EI / l^3 of one element
    P = R0 @ np.array([0.0, 0.0, -1.0]) * 1.2 * EI / L**2      # along the soft local z
    U0 = _newton(ref, _tip(ref, n, np.r_[P, 0, 0, 0]), list(range(ref.dpn)), 6)
    U1 = _newton(tur, _tip(tur, n, np.r_[Q @ P, 0, 0, 0]), list(range(tur.dpn)), 6)
    V0 = U0.reshape(-1, ref.dpn).copy()
    V0[:, 0:3] = V0[:, 0:3] @ Q.T
    V0[:, 3:6] = V0[:, 3:6] @ Q.T          # rotation vectors turn with the beam
    assert np.linalg.norm(V0[-1, :3]) > 0.04 * L          # the channel twists and bends less
    np.testing.assert_allclose(U1.reshape(-1, tur.dpn), V0, atol=1e-8 * L)


@pytest.mark.parametrize("etype", [CR_EB, CR_TS])
def test_curved_cantilever_closes_into_a_semicircle(etype):
    """Quarter circle of chords under the tip moment EI / R: a semicircle of radius R / 2."""
    Rc, n = 1.0, 48
    th = np.linspace(0.0, np.pi / 2, n + 1)
    pts = np.column_stack([Rc * np.sin(th), Rc * (1 - np.cos(th)), np.zeros(n + 1)])
    mesh = _mesh(etype, pts)
    EI = E * RECT["I_z"]                    # in-plane bending, about local (and global) z
    U = _newton(mesh, _tip(mesh, n, [0, 0, 0, 0, 0, EI / Rc]), list(range(6)), 8)
    tip = pts[-1] + U.reshape(-1, 6)[-1, :3]
    # inextensible exact answer (0, R); the chords add an O(1/n^2) error, and the beam stretches a little
    np.testing.assert_allclose(tip, [0.0, Rc, 0.0], atol=2e-3 * Rc)
    assert U.reshape(-1, 6)[-1, 5] == pytest.approx(np.pi / 2, rel=2e-3)   # tangent from +Y to -X


def test_prebent_path_tends_to_the_straight_path():
    """Tip deflection under a large flap load; prebend delta at the tip, parabolic."""
    L, n = 1.0, 16
    x = np.linspace(0.0, L, n + 1)
    straight = None
    gaps = []
    for delta in (0.0, 0.04, 0.02, 0.01):
        pts = np.column_stack([x, np.zeros_like(x), delta * (x / L) ** 2])
        mesh = _mesh(CR_TS, pts)
        EI = E * RECT["I_y"]
        U = _newton(mesh, _tip(mesh, n, [0, 0, -EI / L**2, 0, 0, 0]), list(range(6)), 4)
        tip = pts[-1] + U.reshape(-1, 6)[-1, :3]
        if straight is None:
            straight = tip
            # Bisshopp and Drucker: P L^2 / EI = 1 drops the tip by 0.30172 L and pulls it in by 0.05643 L
            assert straight[2] == pytest.approx(-0.30172 * L, rel=3e-3)
            assert L - straight[0] == pytest.approx(0.05643 * L, rel=1e-2)
            continue
        gaps.append(np.linalg.norm(tip - straight - [0, 0, delta]))
    assert gaps[0] > gaps[1] > gaps[2] > 0
    assert gaps[1] / gaps[2] == pytest.approx(2.0, rel=0.1)     # linear in the prebend


@pytest.mark.parametrize("etype", [CR_EB, CR_TS, CR_BS])
def test_small_load_on_a_kinked_twisted_mesh_is_the_linear_core(etype):
    """A bent, swept and twisted line of elements: the corotational tangent at rest is the linear stiffness."""
    s = np.linspace(0.0, 1.0, 7)
    pts = np.column_stack([s, 0.15 * s**2, -0.1 * s**3])
    twist = np.linspace(0.0, 40.0, 6)
    cr = _mesh(etype, pts, twist=twist)
    lin = _mesh(LIN[etype], pts, twist=twist)
    _, Kt = _assemble(cr, np.zeros(cr.n))
    np.testing.assert_allclose(Kt, lin.K, atol=1e-9 * np.abs(lin.K).max())


@pytest.mark.parametrize("etype", [CR_TS, CR_BS])
def test_pretwisted_corotational_element_near_rest(etype):
    """Pretwist rate through precurvature.txt: at small load, close to the linear pretwisted element."""
    L, n, k = 1.0, 8, np.radians(60.0)     # 60 degrees of pretwist over the beam
    pts = _line([0, 0, 0], [1, 0, 0], L, n)
    twist = np.degrees(k * np.linspace(0.0, L, n + 1)[:-1])
    cr = _mesh(etype, pts, twist=twist, curvature=[k, 0.0, 0.0])
    lin = _mesh(LIN[etype], pts, twist=twist, curvature=[k, 0.0, 0.0])
    F = _tip(cr, n, [200.0, 30.0, -40.0, 5.0, 0.0, 0.0])
    fixed = list(range(cr.dpn))
    _, Kt = _assemble(cr, np.zeros(cr.n))
    free = np.setdiff1d(np.arange(cr.n), fixed)
    U_cr = np.zeros(cr.n)
    U_cr[free] = np.linalg.solve(Kt[np.ix_(free, free)], F[free])
    U_lin = lin.solve(F, fixed)
    # the corotational element is one straight element with its section frame at mid-length; the
    # linear pretwisted element is condensed from 8 sub-elements, so it is the more refined of the two
    np.testing.assert_allclose(U_cr, U_lin, rtol=0, atol=5e-3 * np.abs(U_lin).max())


def test_curved_axis_in_precurvature_is_refused():
    pts = _line([0, 0, 0], [1, 0, 0], 1.0, 2)
    with pytest.raises(NotImplementedError, match="corotational elements are straight"):
        _mesh(CR_TS, pts, curvature=[0.0, 0.0, 0.5])
