"""Exact natural frequencies of a prismatic beam whose 6 x 6 section is fully coupled.

The reference for 1d-fem's eigen solver on anisotropic sections: extension-twist,
bending-twist and shear-bending couplings in the stiffness, and mass-centre offsets in the
mass, all at once. A linear Timoshenko beam (the linearised geometrically exact beam of
Reissner and Simo) is twelve first-order equations along the span in the state
``y = [u, theta, F, M]``::

    u'     = (C f)_F - e_x x theta        theta' = (C f)_M
    F'     = -omega^2 (M q)_u             M'     = -e_x x F - omega^2 (M q)_theta

with ``f = [F, M]``, ``q = [u, theta]`` and ``C`` the section compliance. A prismatic beam
integrates them exactly, ``y(L) = expm(A L) y(0)``; a natural frequency is one where some
``y(0)`` meeting the root condition also meets the tip condition, a zero of the 6 x 6 end
matrix. Moved here from anisotropic-section-model's modal work (its PR 29), where it was
the check of a spectral-element modal beam against BeamDyn.

Everything is in 1d-fem's right-handed beam axes: beam along x, strains
``[eps, kappa_y, kappa_z, gamma_xy, gamma_xz, phi]`` against ``[N, M_y, M_z, V_y, V_z, T]``,
mass on ``[u_x, u_y, u_z, theta_x, theta_y, theta_z]`` (``SectionAnalysis.mass_matrix``).
An element-signed stiffness (``SectionAnalysis.stiffness``) is right-handed after
``ELEMENT_SIGN @ K @ ELEMENT_SIGN``.

Public API
----------
prismatic_frequencies, ENDS.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

ENDS = ("clamped", "free", "pinned")

#: state [u, theta, F, M] from the strain-ordered stiffness: forces [N, V_y, V_z] and moments
#: [T, M_y, M_z] are strain entries [0, 3, 4] and [5, 1, 2]
_FORCE_ORDER = np.array([0, 3, 4, 5, 1, 2])


def _skew(v: np.ndarray) -> np.ndarray:
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def _state_matrix(stiffness: np.ndarray, mass: np.ndarray, omega: float) -> np.ndarray:
    """``A(omega)`` of ``y' = A y`` for the right-handed strain-ordered ``stiffness``."""
    K = 0.5 * (stiffness + stiffness.T)
    # compliance on [F_x, F_y, F_z, M_x, M_y, M_z] against [u' - .., theta']
    C = np.linalg.inv(K[np.ix_(_FORCE_ORDER, _FORCE_ORDER)])
    Mm = 0.5 * (mass + mass.T)
    E = _skew(np.array([1.0, 0.0, 0.0]))
    A = np.zeros((12, 12))
    A[0:6, 6:12] = C
    A[0:3, 3:6] -= E
    A[6:12, 0:6] = -omega ** 2 * Mm
    A[9:12, 6:9] -= E
    return A


def _end_rows(end: str) -> np.ndarray:
    """Indices of the state ``[u, theta, F, M]`` that vanish at an end."""
    if end == "clamped":
        return np.arange(6)
    if end == "free":
        return np.arange(6, 12)
    if end == "pinned":  # translations and twist held, bending moments free
        return np.array([0, 1, 2, 3, 10, 11])
    raise ValueError(f"an end is one of {ENDS}")


def prismatic_frequencies(stiffness, mass, length: float, f_max: float, ends: Sequence[str] = ("clamped", "free"),
                          n_grid: int = 2000) -> np.ndarray:
    """Exact natural frequencies (Hz) below ``f_max`` of a prismatic, fully coupled beam.

    ``stiffness`` is the right-handed 6 x 6 in strain order, ``mass`` the 6 x 6 per length,
    both about the beam axis. ``ends`` are the root and tip conditions: ``"clamped"``,
    ``"free"`` or ``"pinned"`` (translations and twist held). A repeated frequency (a
    symmetric section) is returned once per independent mode. Accurate while ``expm`` stays
    well conditioned: the first ten or so modes of a slender beam.
    """
    from scipy.linalg import expm
    from scipy.optimize import minimize_scalar

    S = np.asarray(stiffness, dtype=float)
    Mm = np.asarray(mass, dtype=float)
    root, tip = _end_rows(ends[0]), _end_rows(ends[1])
    free0 = np.setdiff1d(np.arange(12), root)
    # scale forces and moments so the end matrix is balanced
    EI = np.mean(np.abs(np.diag(S)[1:3]))
    fs = np.ones(12)
    fs[0:3] = 1.0 / length
    fs[6:9] = EI / length ** 2
    fs[9:12] = EI / length
    D_scale = np.diag(1.0 / fs)

    def sv(f):
        Phi = D_scale @ expm(_state_matrix(S, Mm, 2.0 * np.pi * f) * length) @ np.diag(fs)
        s = np.linalg.svd(Phi[np.ix_(tip, free0)], compute_uv=False)
        return s / s[0]

    grid = np.linspace(f_max / n_grid, 1.02 * f_max, n_grid)  # a little past f_max: a root at the end is a dip
    smin = np.array([sv(f)[-1] for f in grid])
    out = []
    for i in range(1, n_grid - 1):
        if smin[i] <= smin[i - 1] and smin[i] <= smin[i + 1]:
            r = minimize_scalar(lambda f: sv(f)[-1], bracket=(grid[i - 1], grid[i], grid[i + 1]), tol=1e-13)
            f = float(r.x)
            s = sv(f)
            # a root is a deep dip; a repeated root has as many singular values that small
            if s[-1] < 1e-2 * max(smin[i - 1], smin[i + 1]):
                out.extend([f] * int(np.sum(s <= 1e3 * s[-1] + 1e-300)))
    return np.array(sorted(f for f in out if f <= f_max))
