# pre_processing/element_library/isotropic_section.py
"""Isotropic section properties, from the section provider.

1D-FEM computes no section properties of its own: the isotropic beam constants of
``section.txt`` and ``material.txt`` (``E``, ``G``, ``rho``, ``A``, ``I_y``, ``I_z``, ``J_t``,
``kappa``, ...) go through the provider's ``SectionAnalysis.from_isotropic``, and every
element reads its ``D`` and its section mass from the result. This module is the one place
that call is made, so the bar, truss, Euler-Bernoulli, Timoshenko and corotational elements
and the element ``section_analysis`` all see the same numbers.

Two geometric torsion constants [m^4], not to be confused:

- ``J_t``, the Saint-Venant torsion constant, independent of the reference axis. It enters
  only the torsional stiffness ``D[5, 5] = G J_t`` [N m^2] (torque per twist rate).
- ``I_x``, the polar second moment of area ``I_p = int (y^2 + z^2) dA`` about the centroid
  (the beam reference axis of a ``section.txt`` element). It enters the twist rotary
  inertia ``M[3, 3] = rho I_p`` [kg m] about that axis, and the Wagner term ``I_p / A``.
  ``section.txt`` often leaves it ``0.0``; a value that is not positive counts as absent
  and the provider takes ``I_y + I_z``.

They are equal only for circular sections; for any other shape ``J_t < I_p``.

Public API
----------
isotropic_section, isotropic_section_of, isotropic_stiffness, isotropic_mass, nodal_mass_diagonal
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Optional

import numpy as np

from onedfem.sections import section_import

__all__ = ["isotropic_section", "isotropic_section_of", "isotropic_stiffness", "isotropic_mass", "nodal_mass_diagonal"]


@lru_cache(maxsize=4096)
def _cached(E, G, A, I_y, I_z, J_t, kappa, I_x, y_sc, z_sc, Gamma, rho, provider):
    SectionAnalysis = section_import("section_analysis", "SectionAnalysis")
    section = SectionAnalysis.from_isotropic(
        E=E, G=G, A=A, I_y=I_y, I_z=I_z, J_t=J_t, kappa=kappa, I_x=I_x,
        y_sc=y_sc, z_sc=z_sc, Gamma=Gamma, rho=rho,
    )
    # one object is shared by every element with these constants: make its arrays read-only
    for value in vars(section).values():
        if isinstance(value, np.ndarray):
            value.flags.writeable = False
    return section


def isotropic_section(
    *,
    E: float,
    G: float,
    A: float,
    I_y: float = 0.0,
    I_z: float = 0.0,
    J_t: float = 0.0,
    kappa: float = 5.0 / 6.0,
    I_x: Optional[float] = None,
    y_sc: float = 0.0,
    z_sc: float = 0.0,
    Gamma: float = 0.0,
    rho: float = 0.0,
):
    """The provider's ``SectionAnalysis.from_isotropic`` for these constants (cached, read-only).

    ``I_x`` that is ``None`` or not positive is passed as ``None`` (polar moment ``I_y + I_z``).
    """
    from onedfem.sections import section_provider

    ix = None if I_x is None or not float(I_x) > 0.0 else float(I_x)
    return _cached(float(E), float(G), float(A), float(I_y), float(I_z), float(J_t), float(kappa), ix,
                   float(y_sc), float(z_sc), float(Gamma), float(rho), section_provider().__name__)


def isotropic_section_of(element):
    """``isotropic_section`` from an element's ``section_array`` and ``material_array``.

    ``section_array``: ``A, I_x, I_y, I_z, J_t`` then optional ``kappa, alpha``,
    ``y_sc, z_sc`` and ``Gamma``; ``kappa`` defaults to 5/6. ``material_array``: ``E, G, nu, rho``.
    """
    sa = np.asarray(element.section_array, dtype=np.float64)
    E, G, _nu, rho = (float(v) for v in element.material_array)

    def at(i: int, default: float = 0.0) -> float:
        return float(sa[i]) if sa.size > i else default

    return isotropic_section(
        E=E, G=G, A=at(0), I_x=at(1), I_y=at(2), I_z=at(3), J_t=at(4),
        kappa=at(5) if sa.size >= 7 else 5.0 / 6.0,
        y_sc=at(7) if sa.size >= 9 else 0.0,
        z_sc=at(8) if sa.size >= 9 else 0.0,
        Gamma=at(9) if sa.size >= 10 else 0.0,
        rho=rho,
    )


def isotropic_stiffness(**constants) -> np.ndarray:
    """Section ``D`` (6, 6) on ``[eps, kappa_y, kappa_z, gamma_xy, gamma_xz, phi]`` (a copy)."""
    return np.array(isotropic_section(**constants).stiffness, dtype=np.float64)


def isotropic_mass(section) -> Dict[str, float]:
    """Mass per length of a section (inertia, about the section origin = beam reference axis).

    ``m = rho A``; ``I_p`` the polar mass moment ``rho int (y^2 + z^2) dA`` (twist inertia, not
    the torsion constant ``J_t``); ``I_yy``, ``I_zz`` the bending rotary inertias.
    """
    mass = dict(getattr(section, "mass", None) or {})
    return {k: float(mass.get(k, 0.0)) for k in ("m", "I_p", "I_yy", "I_zz")}


def nodal_mass_diagonal(section, n_dof: int = 12) -> np.ndarray:
    """Per-DOF mass densities ``mu`` of a two-node element from the section mass.

    Translations ``m``, twist ``I_p`` (polar mass moment, not ``J_t``), bending rotations ``I_yy`` (``theta_y``) and ``I_zz``
    (``theta_z``), node 1 then node 2 on DOFs ``0..5`` and ``6..11``; any further DOFs are 0.
    """
    ms = isotropic_mass(section)
    mu = np.zeros(int(n_dof), dtype=np.float64)
    for i in (0, 1, 2, 6, 7, 8):
        mu[i] = ms["m"]
    for i in (3, 9):
        mu[i] = ms["I_p"]
    for i in (4, 10):
        mu[i] = ms["I_yy"]
    for i in (5, 11):
        mu[i] = ms["I_zz"]
    return mu
