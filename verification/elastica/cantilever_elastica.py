"""
Large-deflection (elastica) analytical solution for cantilever beams.

Uses an iterative fixed-point method: bending moment M(s) is computed from
the current deformed horizontal coordinates x(s), updated each iteration via
numerical integration of the elastica ODE (EI dθ/ds = M). Converges in
~15-25 iterations for δ/L ≤ 30 %.

Supported load_type values:
  "end"        — vertical point load at the free end (tip)
  "mid"        — vertical point load at s = L/2
  "quarter"    — vertical point load at s = L/4
  "udl"        — uniformly distributed vertical load, intensity w per unit length
  "triangular" — linearly increasing load q(s) = w_max · s/L (zero at root, max at tip)
  "parabolic"  — parabolically increasing load q(s) = w_max · (s/L)² (zero at root, max at tip)
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid


def _compute_moment(
    s: np.ndarray,
    x: np.ndarray,
    load_type: str,
    magnitude: float,
    L: float,
) -> np.ndarray:
    n = len(s)
    M = np.zeros(n)
    if load_type == "end":
        x_tip = x[-1]
        M = magnitude * (x_tip - x)
        M[-1] = 0.0
    elif load_type in ("mid", "quarter"):
        a = L / 2.0 if load_type == "mid" else L / 4.0
        x_a = float(np.interp(a, s, x))
        for i in range(n):
            M[i] = magnitude * (x_a - x[i]) if s[i] <= a else 0.0
    elif load_type == "udl":
        # M(si) = w * integral_{si}^{L} (x(s') - x(si)) ds'
        for i in range(n):
            mask = s >= s[i]
            if mask.sum() >= 2:
                M[i] = magnitude * np.trapezoid(x[mask] - x[i], s[mask])
    elif load_type == "triangular":
        # q(s) = w_max * s/L
        for i in range(n):
            mask = s >= s[i]
            if mask.sum() >= 2:
                q = magnitude * s[mask] / L
                M[i] = np.trapezoid(q * (x[mask] - x[i]), s[mask])
    elif load_type == "parabolic":
        # q(s) = w_max * (s/L)^2
        for i in range(n):
            mask = s >= s[i]
            if mask.sum() >= 2:
                q = magnitude * (s[mask] / L) ** 2
                M[i] = np.trapezoid(q * (x[mask] - x[i]), s[mask])
    else:
        raise ValueError(f"Unknown load_type: {load_type!r}")
    return M


def elastica_tip_deflection(
    L: float,
    E: float,
    I: float,
    magnitude: float,
    load_type: str,
    n_pts: int = 200,
    n_iter: int = 30,
) -> float:
    """Return the absolute tip deflection (m) for a cantilever in the large-deflection regime.

    Parameters
    ----------
    L : float
        Beam length (m).
    E : float
        Young's modulus (Pa).
    I : float
        Second moment of area about the bending axis (m⁴).
    magnitude : float
        Load magnitude: force P (N) for point-load types; intensity w or w_max (N/m)
        for distributed types. Always positive.
    load_type : str
        One of "end", "mid", "quarter", "udl", "triangular", "parabolic".
    n_pts : int
        Number of arc-length sample points (default 200).
    n_iter : int
        Maximum fixed-point iterations (default 30).

    Returns
    -------
    float
        Absolute tip deflection |y(L)| in metres.
    """
    EI = E * I
    s = np.linspace(0.0, L, n_pts)
    x = s.copy()  # initial guess: undeformed geometry

    for _ in range(n_iter):
        M = _compute_moment(s, x, load_type, magnitude, L)
        theta = cumulative_trapezoid(M / EI, s, initial=0.0)
        x_new = cumulative_trapezoid(np.cos(theta), s, initial=0.0)
        y = cumulative_trapezoid(np.sin(theta), s, initial=0.0)
        if np.max(np.abs(x_new - x)) < 1e-9 * L:
            break
        x = x_new

    return float(abs(y[-1]))
