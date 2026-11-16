# pre_processing/element_library/added_mass.py
"""Added mass in water as a labelled strip-heave model.

This is a potential-flow model, not a measured frequency in water. No still-water
frequency, hydrodynamic coefficient, or empirical factor is implied or fitted.

Labelled model
    ``strip_heave``: two-dimensional unbounded-fluid added mass of a flat plate
    (thin-airfoil heave), applied as a consistent mass on the section-normal
    (heave) translation only. Edgewise (chordwise) added mass of a thin plate is
    left at zero. Axial, torsion, and extra DOFs (warping, GBT) carry none.

Formula
    Per unit span, ``m_a' = rho_fluid * pi * (chord / 2)**2``
    ``         = rho_fluid * pi * chord**2 / 4``.

Source
    Newman, J. N. (1977). *Marine Hydrodynamics*. MIT Press, Cambridge, MA.
    Two-dimensional added mass of a flat plate of chord ``c = 2b`` for motion
    normal to the plate is ``rho * pi * b**2`` per unit span (the circumscribed
    cylinder). Theodorsen, T. (1935). *General Theory of Aerodynamic Instability
    and the Mechanism of Flutter*, NACA Report 496: the apparent mass of a thin
    airfoil of semi-chord ``b`` is ``rho * pi * b**2`` per unit span.

Weak form
    ``M_a += N_h.T @ N_h * m_a' * w_g * detJ`` at Gauss points, with ``N_h`` the
    heave (local ``y`` or ``z``) translation of the element's own field
    interpolation (Hermite bending for Euler-Bernoulli, linear otherwise),
    rotated to global DOFs by the node frames. ``detJ = L/2``. Straight elements
    only (rotated and twisted allowed); a curved element raises.

The default is off: no ``[AddedMass]`` block, or ``model`` unset, leaves the
structural mass unchanged.

Public API
----------
STRIP_HEAVE_FORMULA, STRIP_HEAVE_SOURCE, AddedMass, AddedMassGaussPointData,
AddedMassObject, strip_heave_added_mass_per_length, element_added_mass,
added_mass_from_settings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import numpy.typing as npt

from onedfem.pre_processing.element_library.beam_load_interpolation import beam_field_interpolation
from onedfem.pre_processing.element_library.element_frame import node_dof_rotation

__all__ = [
    "STRIP_HEAVE_FORMULA",
    "STRIP_HEAVE_SOURCE",
    "AddedMass",
    "AddedMassGaussPointData",
    "AddedMassObject",
    "strip_heave_added_mass_per_length",
    "element_added_mass",
    "added_mass_from_settings",
]

STRIP_HEAVE_MODEL = "strip_heave"

STRIP_HEAVE_FORMULA = "rho_fluid * pi * (chord/2)**2"

STRIP_HEAVE_SOURCE = (
    "Newman, J. N. (1977). Marine Hydrodynamics. MIT Press. "
    "Two-dimensional added mass of a flat plate of chord c = 2b for motion "
    "normal to the plate is rho * pi * b**2 per unit span. "
    "Theodorsen, T. (1935). NACA Report 496: apparent mass of a thin airfoil "
    "of semi-chord b is rho * pi * b**2 per unit span."
)


def strip_heave_added_mass_per_length(rho_fluid: float, chord: float) -> float:
    """Heave added mass per unit span of a 2-D flat plate of chord ``chord``.

    ``m_a' = rho_fluid * pi * (chord / 2)**2``. ``rho_fluid`` and ``chord`` are
    supplied by the caller; this function does not choose seawater density or a
    blade chord.
    """
    rho = float(rho_fluid)
    c = float(chord)
    if rho < 0.0:
        raise ValueError(f"rho_fluid must be >= 0, got {rho_fluid!r}")
    if c < 0.0:
        raise ValueError(f"chord must be >= 0, got {chord!r}")
    return rho * np.pi * (c / 2.0) ** 2


@dataclass(frozen=True)
class AddedMass:
    """Labelled strip-heave added mass: fluid density, chord, and local heave axis.

    ``mass_per_length`` is always the Newman/Theodorsen plate formula from
    ``rho_fluid`` and ``chord``. It is not a measured still-water frequency.
    """

    rho_fluid: float
    chord: float
    heave_local_axis: str = "z"
    model: str = STRIP_HEAVE_MODEL

    def __post_init__(self):
        axis = str(self.heave_local_axis).strip().lower()
        if axis not in ("y", "z"):
            raise ValueError(f"heave_local_axis must be 'y' or 'z', got {self.heave_local_axis!r}")
        object.__setattr__(self, "heave_local_axis", axis)
        model = str(self.model).strip().lower()
        if model != STRIP_HEAVE_MODEL:
            raise ValueError(
                f"added-mass model must be {STRIP_HEAVE_MODEL!r} (the labelled strip-heave "
                f"plate formula); got {self.model!r}"
            )
        object.__setattr__(self, "model", model)
        object.__setattr__(self, "rho_fluid", float(self.rho_fluid))
        object.__setattr__(self, "chord", float(self.chord))
        if self.rho_fluid < 0.0:
            raise ValueError(f"rho_fluid must be >= 0, got {self.rho_fluid!r}")
        if self.chord < 0.0:
            raise ValueError(f"chord must be >= 0, got {self.chord!r}")

    @property
    def mass_per_length(self) -> float:
        """``rho_fluid * pi * (chord/2)**2`` (kg/m)."""
        return strip_heave_added_mass_per_length(self.rho_fluid, self.chord)

    @property
    def formula(self) -> str:
        return STRIP_HEAVE_FORMULA

    @property
    def source(self) -> str:
        return STRIP_HEAVE_SOURCE


@dataclass
class AddedMassGaussPointData:
    """One Gauss station of the strip-heave added-mass operator.

    ``Phi_heave`` (n,) maps global element DOFs to the heave translation;
    ``direction`` (3,) is that heave axis in global coordinates.
    """

    xi: float
    weight: float
    jacobian: float
    direction: npt.NDArray[np.float64]
    Phi_heave: npt.NDArray[np.float64]
    mass_per_length: float


@dataclass
class AddedMassObject:
    """Strip-heave added-mass matrix of one element (global DOFs), with Gauss data."""

    element_id: int
    model: AddedMass
    M_a: npt.NDArray[np.float64]
    gauss_data: List[AddedMassGaussPointData] = field(default_factory=list)


def _heave_index(axis: str) -> int:
    return 1 if axis == "y" else 2


def _twist_rotation(R0: np.ndarray, R1: np.ndarray):
    """``R(t)``, ``t`` in [0, 1]: ``R0`` turned about its ``x`` by ``t`` times the twist to ``R1``."""
    beta = float(np.arctan2(R0[:, 2] @ R1[:, 1], R0[:, 1] @ R1[:, 1]))

    def R(t: float) -> np.ndarray:
        c, s = np.cos(beta * t), np.sin(beta * t)
        return R0 @ np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])

    return R


def element_added_mass(
    element,
    p0: npt.ArrayLike,
    p1: npt.ArrayLike,
    model: AddedMass,
    R0: Optional[np.ndarray] = None,
    R1: Optional[np.ndarray] = None,
    n_dof: Optional[int] = None,
    n_gauss: int = 4,
) -> AddedMassObject:
    """Strip-heave added mass of one straight element from ``p0`` to ``p1`` (global).

    Frames default to the element's own (``_frame_R0`` / ``_frame_R1``), else global
    axes. ``n_dof`` defaults to the element's DOF count; DOFs past the first 12
    get no added mass. Curved elements raise.
    """
    if getattr(element, "_frame_arc", None) is not None:
        raise NotImplementedError(
            f"Element {element.element_id}: strip-heave added mass is not available "
            "for curved elements"
        )
    kind = getattr(element, "bending_interpolation", None)
    if kind not in ("hermite", "linear"):
        raise NotImplementedError(
            f"Element {element.element_id}: strip-heave added mass is a beam model; "
            f"{type(element).__name__} has no bending interpolation"
        )
    p0 = np.asarray(p0, dtype=np.float64)
    p1 = np.asarray(p1, dtype=np.float64)
    if R0 is None:
        R0 = getattr(element, "_frame_R0", None)
    R0 = np.eye(3) if R0 is None else np.asarray(R0, dtype=np.float64)
    if R1 is None:
        R1 = getattr(element, "_frame_R1", None)
    R1 = R0 if R1 is None else np.asarray(R1, dtype=np.float64)
    if n_dof is None:
        n_dof = int(np.asarray(element.assemble_global_dof_indices()).size)
    L = float(np.linalg.norm(p1 - p0))
    T = node_dof_rotation(R0, R1, n_dof)
    Rt = _twist_rotation(R0, R1)
    h = _heave_index(model.heave_local_axis)
    m_a = model.mass_per_length

    M_a = np.zeros((n_dof, n_dof))
    cache: List[AddedMassGaussPointData] = []
    xg, wg = np.polynomial.legendre.leggauss(n_gauss)
    detJ = 0.5 * L
    for xi, w in zip(xg, wg):
        t = 0.5 * (1.0 + xi)
        Rs = Rt(t)
        N = beam_field_interpolation(xi, L, kind)
        F = np.zeros((6, n_dof))
        F[:, :12] = N
        F = np.vstack([R0 @ F[:3], R0 @ F[3:]]) @ T
        direction = Rs[:, h]
        Phi_h = direction @ F[:3]
        M_a += w * detJ * m_a * np.outer(Phi_h, Phi_h)
        cache.append(
            AddedMassGaussPointData(
                xi=float(xi),
                weight=float(w),
                jacobian=detJ,
                direction=direction.copy(),
                Phi_heave=Phi_h.copy(),
                mass_per_length=m_a,
            )
        )
    return AddedMassObject(
        element_id=int(element.element_id),
        model=model,
        M_a=0.5 * (M_a + M_a.T),
        gauss_data=cache,
    )


def added_mass_from_settings(simulation_settings: Dict[str, Any]) -> Optional[AddedMass]:
    """The ``[AddedMass]`` settings as an ``AddedMass``, or None when the model is off.

    Off when the block is absent, or when ``model``, ``rho_fluid`` and ``chord``
    are all unset. ``model = strip_heave`` with ``rho_fluid`` and ``chord`` turns
    it on. No seawater density or chord is filled in.
    """
    cfg = (simulation_settings or {}).get("added_mass") or {}
    model = cfg.get("model")
    rho = cfg.get("rho_fluid")
    chord = cfg.get("chord")
    if model is None and rho is None and chord is None:
        return None
    label = STRIP_HEAVE_MODEL if model is None else str(model).strip().lower()
    if label != STRIP_HEAVE_MODEL:
        raise ValueError(
            f"[AddedMass] model must be {STRIP_HEAVE_MODEL!r} (Newman/Theodorsen "
            f"strip heave); got {model!r}"
        )
    if rho is None or chord is None:
        raise ValueError(
            "[AddedMass] strip_heave needs rho_fluid and chord (the labelled "
            f"formula is {STRIP_HEAVE_FORMULA}; no coefficient is invented)"
        )
    return AddedMass(
        rho_fluid=float(rho),
        chord=float(chord),
        heave_local_axis=str(cfg.get("heave_local_axis", "z")),
        model=label,
    )
