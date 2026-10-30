# pre_processing/element_library/rotating_frame.py
"""Beams in a frame rotating at a constant rate: centrifugal load, spin softening and Coriolis terms.

A structure spinning at ``Omega`` (rad/s, vector ``Omega = omega a`` about an axis through
``origin``) is analysed in the rotating frame. A material point at ``x = X + u`` (from the
origin) carries the inertia load ``-rho (Omega x (Omega x x) + 2 Omega x u_dot)``.

Kinematics
    A beam section is rigid in its own plane. Its mass is lumped at four points per Gauss
    station (``section_point_masses``) that reproduce the section's mass, mass centre and
    second moments of mass, so the rigid-section mass matrix
    (``SectionAnalysis.mass_matrix``) is recovered exactly. A point at section offset ``r``
    moves by ``u_p = u + theta x r + 1/2 theta x (theta x r)`` with the element's own field
    interpolation of ``u`` and ``theta`` (``beam_field_interpolation``: Hermite bending for
    Euler-Bernoulli elements, linear for the others), so ``u_p = Phi U_e`` to first order.

Element terms (``W = [Omega x]``, ``f_p = -m_p W^2 X_p`` the centrifugal force on a point)
    - centrifugal load    ``f_c = sum Phi^T f_p``
    - spin softening      ``K_soft = sum m_p Phi^T W^2 Phi``  (negative semi-definite)
    - rotation load stiffness (the second-order ``theta`` term of ``u_p`` under ``f_p``)
      ``K_rot = sum Theta^T [(f_p . r_p) I - 1/2 (f_p r_p^T + r_p f_p^T)] Theta``, which
      gives the propeller moment of a section turning out of the plane of rotation
    - Coriolis (gyroscopic) ``G = sum 2 m_p Phi^T W Phi``  (skew)

    Each sum runs over the element's Gauss points (weight ``w``, Jacobian ``L/2``) and the
    four point masses there. Every term is formed in that Gauss loop, and the Gauss data
    (positions, masses, ``Phi``) are kept in the ``RotatingFrameObject`` so post-processing
    reads the same operators that built the matrices.

Prestress
    The centrifugal tension stiffens the beam through ``K_sigma`` of the static centrifugal
    state, which the analysis runner adds (``linear_geometric_stiffness_matrix`` of each
    element at that state). Free vibration in the rotating frame is then
    ``M q'' + G q' + (K + K_sigma + K_soft + K_rot) q = 0``.

Limitations
    Straight elements (rotated and twisted allowed): a curved element raises. Warping and
    other higher-order DOFs carry no rotating-frame terms. The section points use the
    element's start frame, turned about the element axis by the twist between its two node
    frames.

Public API
----------
Rotation, RotatingFrameGaussPointData, RotatingFrameObject, section_point_masses,
element_rotating_frame, rotation_from_settings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import numpy.typing as npt

from pre_processing.element_library.beam_load_interpolation import beam_field_interpolation
from pre_processing.element_library.element_frame import node_dof_rotation

__all__ = [
    "Rotation",
    "RotatingFrameGaussPointData",
    "RotatingFrameObject",
    "section_point_masses",
    "element_rotating_frame",
    "rotation_from_settings",
]


def _skew(v: npt.ArrayLike) -> npt.NDArray[np.float64]:
    a, b, c = np.asarray(v, dtype=np.float64)
    return np.array([[0.0, -c, b], [c, 0.0, -a], [-b, a, 0.0]])


@dataclass(frozen=True)
class Rotation:
    """Constant rotation ``omega`` (rad/s) about ``axis`` through ``origin`` (global axes)."""

    omega: float
    axis: Tuple[float, float, float] = (0.0, 0.0, 1.0)
    origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self):
        a = np.asarray(self.axis, dtype=np.float64)
        n = float(np.linalg.norm(a))
        if n == 0.0:
            raise ValueError("rotation axis is zero")
        object.__setattr__(self, "axis", tuple(float(v) for v in a / n))
        object.__setattr__(self, "origin", tuple(float(v) for v in np.asarray(self.origin, dtype=np.float64)))

    @property
    def vector(self) -> npt.NDArray[np.float64]:
        return float(self.omega) * np.asarray(self.axis)

    @property
    def W(self) -> npt.NDArray[np.float64]:
        """``[Omega x]``."""
        return _skew(self.vector)

    def with_omega(self, omega: float) -> "Rotation":
        return Rotation(float(omega), self.axis, self.origin)


@dataclass
class RotatingFrameGaussPointData:
    """One Gauss station of an element in the rotating frame.

    ``positions`` (4, 3) are the point masses from the rotation origin and ``offsets`` (4, 3)
    their section offsets ``r``, both global; ``masses`` (4,) per unit length; ``Phi``
    (4, 3, n) map the element DOFs to each point's displacement; ``Theta`` (3, n) to the
    section rotation; ``forces`` (4, 3) are the centrifugal forces per unit length.
    """

    xi: float
    weight: float
    jacobian: float
    positions: npt.NDArray[np.float64]
    offsets: npt.NDArray[np.float64]
    masses: npt.NDArray[np.float64]
    Phi: npt.NDArray[np.float64]
    Theta: npt.NDArray[np.float64]
    forces: npt.NDArray[np.float64]


@dataclass
class RotatingFrameObject:
    """Rotating-frame element terms (global DOFs), with their Gauss data."""

    element_id: int
    rotation: Rotation
    f_c: npt.NDArray[np.float64]
    K_soft: npt.NDArray[np.float64]
    K_rot: npt.NDArray[np.float64]
    G: npt.NDArray[np.float64]
    gauss_data: List[RotatingFrameGaussPointData] = field(default_factory=list)


def section_point_masses(mass: Dict[str, float]) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Four point masses ``(yz (4, 2), m (4,))`` with the section's mass, mass centre and second moments.

    ``mass`` has ``m``, ``y_m``, ``z_m`` and the second moments about the section origin
    ``I_yy = int rho z^2``, ``I_zz = int rho y^2``, ``I_yz = int rho y z`` (as
    ``SectionAnalysis.mass``). The points sit at the mass centre plus ``+-a`` along each
    principal axis of the moments about the mass centre (eigenvalue ``lambda``),
    ``a = sqrt(2 lambda / m)``.
    """
    m = float(mass["m"])
    if m <= 0.0:
        return np.zeros((4, 2)), np.zeros(4)
    c = np.array([float(mass.get("y_m", 0.0)), float(mass.get("z_m", 0.0))])
    # second moments about the mass centre in (y, z): [[int y^2, int yz], [int yz, int z^2]]
    J = np.array([[float(mass.get("I_zz", 0.0)), float(mass.get("I_yz", 0.0))],
                  [float(mass.get("I_yz", 0.0)), float(mass.get("I_yy", 0.0))]]) - m * np.outer(c, c)
    lam, V = np.linalg.eigh(J)
    lam = np.clip(lam, 0.0, None)
    a = np.sqrt(2.0 * lam / m)
    pts = np.array([c + a[0] * V[:, 0], c - a[0] * V[:, 0], c + a[1] * V[:, 1], c - a[1] * V[:, 1]])
    return pts, np.full(4, 0.25 * m)


def _section_mass(element) -> Dict[str, float]:
    """Section mass properties of an element (second moments about the section origin)."""
    sa = getattr(element, "section_analysis", None)
    if getattr(element, "has_thin_walled_section", False) and sa is not None and getattr(sa, "mass", None):
        mp = dict(sa.mass)
        return {k: float(mp.get(k, 0.0)) for k in ("m", "y_m", "z_m", "I_yy", "I_zz", "I_yz")}
    rho = float(element.material_array[3])
    return {"m": rho * float(element.A), "y_m": 0.0, "z_m": 0.0,
            "I_yy": rho * float(element.I_y), "I_zz": rho * float(element.I_z), "I_yz": 0.0}


def _bending_kind(element) -> str:
    return "hermite" if "EulerBernoulli" in type(element).__name__ else "linear"


def _twist_rotation(R0: np.ndarray, R1: np.ndarray):
    """``R(t)``, ``t`` in [0, 1]: ``R0`` turned about its ``x`` by ``t`` times the twist to ``R1``."""
    beta = float(np.arctan2(R0[:, 2] @ R1[:, 1], R0[:, 1] @ R1[:, 1]))

    def R(t: float) -> np.ndarray:
        c, s = np.cos(beta * t), np.sin(beta * t)
        return R0 @ np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])

    return R


def element_rotating_frame(element, p0: npt.ArrayLike, p1: npt.ArrayLike, rotation: Rotation,
                           R0: Optional[np.ndarray] = None, R1: Optional[np.ndarray] = None,
                           n_dof: Optional[int] = None, n_gauss: int = 4) -> RotatingFrameObject:
    """Rotating-frame terms of one straight element from ``p0`` to ``p1`` (global), node frames ``R0``, ``R1``.

    Frames default to the element's own (``_frame_R0`` / ``_frame_R1``), else global axes.
    ``n_dof`` defaults to the element's DOF count; DOFs past the first 12 get no terms.
    """
    if getattr(element, "_frame_arc", None) is not None:
        raise NotImplementedError(
            f"Element {element.element_id}: rotating-frame terms are not available for curved elements")
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
    kind = _bending_kind(element)
    yz, mp = section_point_masses(_section_mass(element))
    W = rotation.W
    W2 = W @ W
    origin = np.asarray(rotation.origin)

    f_c = np.zeros(n_dof)
    K_soft = np.zeros((n_dof, n_dof))
    K_rot = np.zeros((n_dof, n_dof))
    G = np.zeros((n_dof, n_dof))
    cache: List[RotatingFrameGaussPointData] = []
    xg, wg = np.polynomial.legendre.leggauss(n_gauss)
    detJ = 0.5 * L
    for xi, w in zip(xg, wg):
        t = 0.5 * (1.0 + xi)
        Rs = Rt(t)
        # global [u, theta] at the station from global element DOFs: rotate the local fields
        # (element axes R0, as the element kernels) back to global axes
        N = beam_field_interpolation(xi, L, kind)          # (6, 12) local fields of local DOFs
        F = np.zeros((6, n_dof))
        F[:, :12] = N
        F = np.vstack([R0 @ F[:3], R0 @ F[3:]]) @ T        # global fields of global DOFs
        U_ax, Theta = F[:3], F[3:]
        axis_pt = p0 + t * (p1 - p0) - origin
        offs = (Rs[:, 1:3] @ yz.T).T                       # (4, 3) section offsets, global
        pos = axis_pt[None, :] + offs
        Phi = np.stack([U_ax - _skew(r) @ Theta for r in offs])   # theta x r = -[r x] theta
        forces = -(mp[:, None] * (W2 @ pos.T).T)          # centrifugal force per length
        for k in range(4):
            Pk = Phi[k]
            f_c += w * detJ * Pk.T @ forces[k]
            K_soft += w * detJ * mp[k] * Pk.T @ W2 @ Pk
            G += w * detJ * 2.0 * mp[k] * Pk.T @ W @ Pk
            r, f = offs[k], forces[k]
            Kt = (f @ r) * np.eye(3) - 0.5 * (np.outer(f, r) + np.outer(r, f))
            K_rot += w * detJ * Theta.T @ Kt @ Theta
        cache.append(RotatingFrameGaussPointData(
            xi=float(xi), weight=float(w), jacobian=detJ, positions=pos, offsets=offs, masses=mp.copy(),
            Phi=Phi, Theta=Theta, forces=forces))
    return RotatingFrameObject(
        element_id=int(element.element_id), rotation=rotation, f_c=f_c,
        K_soft=0.5 * (K_soft + K_soft.T), K_rot=0.5 * (K_rot + K_rot.T), G=0.5 * (G - G.T),
        gauss_data=cache)


def rotation_from_settings(simulation_settings: Dict[str, Any]) -> Optional[Rotation]:
    """The ``[Rotation]`` settings as a ``Rotation``, or None when absent or ``omega`` is zero."""
    cfg = (simulation_settings or {}).get("rotation") or {}
    omega = cfg.get("omega")
    if omega is None and cfg.get("rpm") is not None:
        omega = float(cfg["rpm"]) * 2.0 * np.pi / 60.0
    if omega is None or float(omega) == 0.0:
        return None
    return Rotation(float(omega), tuple(cfg.get("axis", (0.0, 0.0, 1.0))), tuple(cfg.get("origin", (0.0, 0.0, 0.0))))
