# examples/thin_walled/gbt_buckling/finite_strip.py
"""Semi-analytical finite strip method: the reference for the GBT buckling example.

An independent model of the same column. The walls are split into flat strips; along the
column every strip displaces as one half-wave of a simply supported plate,

    u = U(s) cos(pi x / a),   v = V(s) sin(pi x / a),   w = W(s) sin(pi x / a),

with ``U``, ``V`` linear and ``W`` a cubic Hermite across the strip (Cheung 1976; the
method of CUFSM). Membrane and plate bending are both full Kirchhoff plate theory, with no
GBT hypothesis (the walls may shear and stretch transversely). The prestress is a uniform
longitudinal compression, which acts on the slopes of ``v`` and ``w``.

``critical_stress(walls, t, E, nu, a)`` is the lowest buckling stress for half-wavelength
``a``; a simply supported column of length ``L`` buckles at the lowest over ``a = L / m``.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np
import scipy.linalg as sla


def _hermite(x: float, le: float):
    H = np.array([1 - 3 * x**2 + 2 * x**3, le * (x - 2 * x**2 + x**3), 3 * x**2 - 2 * x**3, le * (-x**2 + x**3)])
    dH = np.array([-6 * x + 6 * x**2, le * (1 - 4 * x + 3 * x**2), 6 * x - 6 * x**2, le * (-2 * x + 3 * x**2)]) / le
    d2H = np.array([-6 + 12 * x, le * (-4 + 6 * x), 6 - 12 * x, le * (-2 + 6 * x)]) / le**2
    return H, dH, d2H


def critical_stress(walls: Sequence[Tuple[Tuple[float, float], Tuple[float, float]]], t: float, E: float,
                    nu: float, a: float, n_sub: int = 12) -> float:
    """Lowest compressive buckling stress of one half-wave ``a`` (isotropic walls of thickness ``t``)."""
    nodes = []

    def node(p):
        for i, q in enumerate(nodes):
            if np.allclose(q, p, atol=1e-12):
                return i
        nodes.append(np.asarray(p, float))
        return len(nodes) - 1

    strips = []
    for p0, p1 in walls:
        p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
        ids = [node(p0 + (p1 - p0) * k / n_sub) for k in range(n_sub + 1)]
        strips += list(zip(ids[:-1], ids[1:]))
    nd = 4 * len(nodes)                         # per node: U, Y, Z, beta (rotation about x)
    K = np.zeros((nd, nd))
    Kg = np.zeros((nd, nd))
    k = np.pi / a
    A = E * t / (1 - nu**2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
    D = A * t**2 / 12
    xg, wg = np.polynomial.legendre.leggauss(6)
    for i, j in strips:
        d = nodes[j] - nodes[i]
        le = float(np.linalg.norm(d))
        tt = d / le
        nn = np.array([-tt[1], tt[0]])           # dw/ds = beta for this normal
        dofs = np.r_[4 * i:4 * i + 4, 4 * j:4 * j + 4]
        Ke = np.zeros((8, 8))
        Ge = np.zeros((8, 8))
        for x, w in zip(xg, wg):
            xi, wt = 0.5 * (1 + x), 0.5 * w * le
            U = np.zeros(8)
            U[[0, 4]] = [1 - xi, xi]
            Us = np.zeros(8)
            Us[[0, 4]] = [-1 / le, 1 / le]
            V = np.zeros(8)
            V[1:3], V[5:7] = (1 - xi) * tt, xi * tt
            Vs = np.zeros(8)
            Vs[1:3], Vs[5:7] = -tt / le, tt / le
            H, dH, d2H = _hermite(xi, le)
            Wr = []
            for h in (H, dH, d2H):
                r = np.zeros(8)
                r[1:3], r[3], r[5:7], r[7] = h[0] * nn, h[1], h[2] * nn, h[3]
                Wr.append(r)
            W, Ws, Wss = Wr
            # sin terms: eps_x, eps_s, kappa_x, kappa_s; cos terms: gamma_xs, kappa_xs
            Em = np.vstack([-k * U, Vs])
            Eb = np.vstack([k**2 * W, -Wss])
            Gm, Gb = Us + k * V, -2 * k * Ws
            Ke += wt * (Em.T @ A[:2, :2] @ Em + A[2, 2] * np.outer(Gm, Gm)
                        + Eb.T @ D[:2, :2] @ Eb + D[2, 2] * np.outer(Gb, Gb))
            Ge += wt * t * k**2 * (np.outer(V, V) + np.outer(W, W))
        K[np.ix_(dofs, dofs)] += Ke
        Kg[np.ix_(dofs, dofs)] += Ge
    mu = sla.eigh(Kg, K, eigvals_only=True)      # mu = 1 / sigma; K is positive definite
    return float(1.0 / mu.max())


def column_stress(walls, t: float, E: float, nu: float, length: float, n_sub: int = 12,
                  shortest_wave: float = 0.02) -> Tuple[float, int]:
    """Lowest buckling stress of a simply supported column of ``length`` and its number of half-waves.

    Half-waves down to ``shortest_wave`` are tried.
    """
    max_waves = int(length / shortest_wave) + 1
    s = [critical_stress(walls, t, E, nu, length / m, n_sub) for m in range(1, max_waves + 1)]
    m = int(np.argmin(s))
    return s[m], m + 1
