# post_processing/limit_states.py
"""Limit-state results a caller reads from a solved beam.

A static solve already returns nodal displacements. ``deflection_profile`` samples
that field along the beam. ``extreme_fibre_profile`` turns the same motion into
axial strain and homogeneous beam stress on the outer profile. Buckling load
factors and blade vibration modes are the eigenpairs the buckling and eigen
solvers already return: ``buckling_result`` and ``blade_vibration_modes`` only
read them. They do not solve again, and they do not assemble a rotor.

Public API
----------
DeflectionProfile, deflection_profile, ExtremeFibreProfile, extreme_fibre_profile,
BucklingResult, buckling_result, BladeVibrationModes, blade_vibration_modes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Union

import numpy as np

from onedfem.post_processing.graphical_visualisers.beam_field_sampling import interpolation_kind
from onedfem.post_processing.stations import axis_point
from onedfem.pre_processing.element_library.beam_load_interpolation import beam_field_interpolation
from onedfem.pre_processing.element_library.framed_element import local_state
from onedfem.stage_log import stage

__all__ = [
    "DeflectionProfile",
    "deflection_profile",
    "ExtremeFibreProfile",
    "extreme_fibre_profile",
    "BucklingResult",
    "buckling_result",
    "BladeVibrationModes",
    "blade_vibration_modes",
]

_H = 1e-6


@dataclass
class DeflectionProfile:
    """Displacement of a static solution sampled along the beam axis.

    Attributes
    ----------
    s : ndarray (n,)
        Arc length from the first node, metres.
    position : ndarray (n, 3)
        Reference-axis coordinates of each sample.
    displacement : ndarray (n, 6)
        Right-handed ``[u, v, w, theta_x, theta_y, theta_z]`` at each sample.
    """

    s: np.ndarray
    position: np.ndarray
    displacement: np.ndarray


@dataclass
class ExtremeFibreProfile:
    """Axial strain and beam stress on the outer profile, at one station per element.

    Strain is the right-handed fibre extension
    ``eps = du/dx + z d(theta_y)/dx - y d(theta_z)/dx``. Stress is the homogeneous
    Navier stress ``E eps`` from that element's modulus. Rows follow the elements.
    Columns follow ``outer_profile``. ``extreme_strain`` and ``extreme_stress`` are
    ``(min, max)`` over the profile at that station.

    Attributes
    ----------
    s : ndarray (n_stations,)
        Arc length of each station.
    position : ndarray (n_stations, 3)
        Axis coordinates of each station.
    y, z : ndarray (n_points,)
        Outer-profile coordinates in the section frame.
    strain, stress : ndarray (n_stations, n_points)
    extreme_strain, extreme_stress : ndarray (n_stations, 2)
        Minimum and maximum over the outer profile.
    """

    s: np.ndarray
    position: np.ndarray
    y: np.ndarray
    z: np.ndarray
    strain: np.ndarray
    stress: np.ndarray
    extreme_strain: np.ndarray
    extreme_stress: np.ndarray


@dataclass
class BucklingResult:
    """Linear buckling eigenpairs already computed for a reference load.

    Attributes
    ----------
    load_factors : ndarray (n_modes,)
        Multipliers of the reference load. The lowest positive factor is the
        critical load relative to that reference.
    mode_shapes : ndarray
        Buckling modes in the solver's layout, one mode per column when the
        array is two-dimensional.
    """

    load_factors: np.ndarray
    mode_shapes: np.ndarray


@dataclass
class BladeVibrationModes:
    """Undamped vibration modes of the blade beam already computed by an eigen solve.

    These are the modes of the beam that was solved. A rotor is not assembled
    and the modes are not computed again.

    Attributes
    ----------
    frequencies_hz : ndarray (n_modes,)
    mode_shapes : ndarray
        Mode shapes in the solver's layout, one mode per column when the array
        is two-dimensional.
    """

    frequencies_hz: np.ndarray
    mode_shapes: np.ndarray


def _coordinates(coordinates) -> np.ndarray:
    coords = np.asarray(coordinates, dtype=np.float64)
    if coords.ndim == 1:
        coords = np.column_stack([coords, np.zeros(coords.shape[0]), np.zeros(coords.shape[0])])
    if coords.ndim != 2 or coords.shape[1] != 3 or coords.shape[0] < 2:
        raise ValueError("coordinates must be (n_nodes,) or (n_nodes, 3), with at least two nodes")
    return coords


def _nodal_motion(U, n_nodes: int) -> np.ndarray:
    flat = np.asarray(U, dtype=np.float64).ravel()
    if flat.size % n_nodes != 0:
        raise ValueError(f"displacement length {flat.size} is not a multiple of {n_nodes} nodes")
    per_node = flat.size // n_nodes
    if per_node < 6:
        raise ValueError(f"each node needs 6 rigid DOFs, got {per_node}")
    return flat.reshape(n_nodes, per_node)[:, :6]


def deflection_profile(
    coordinates,
    U,
    element_type: str,
    *,
    n_per_element: int = 11,
) -> DeflectionProfile:
    """Sample a static displacement along the beam with the element interpolation.

    ``coordinates`` are the node positions in order. ``U`` is the static solution
    in node-major order (a vector of length ``n_nodes * dof_per_node``, or that
    shape). ``element_type`` selects Hermite or linear bending interpolation.
    Samples include both element ends; the shared node is kept once. Extra DOFs
    (warping, cross-section modes) are not part of the rigid motion.

    Parameters
    ----------
    coordinates : array
        Node coordinates, ``(n_nodes, 3)`` or span positions ``(n_nodes,)``.
    U : array
        Static nodal displacement.
    element_type : str
        Element type name, as ``beam_curve`` takes it.
    n_per_element : int
        Samples on each element, including the two nodes. At least 2.
    """
    if int(n_per_element) < 2:
        raise ValueError("n_per_element must be at least 2")
    coords = _coordinates(coordinates)
    nodal = _nodal_motion(U, coords.shape[0])
    kind = interpolation_kind(element_type)
    xi = np.linspace(-1.0, 1.0, int(n_per_element))
    positions = []
    motion = []
    for e in range(coords.shape[0] - 1):
        p0, p1 = coords[e], coords[e + 1]
        length = float(np.linalg.norm(p1 - p0))
        if length <= 0.0:
            raise ValueError(f"element {e} has zero length")
        Ue = np.concatenate([nodal[e], nodal[e + 1]])
        samples = range(xi.size) if e == 0 else range(1, xi.size)
        for i in samples:
            t = 0.5 * (1.0 + float(xi[i]))
            positions.append((1.0 - t) * p0 + t * p1)
            motion.append(beam_field_interpolation(float(xi[i]), length, kind) @ Ue)
    position = np.vstack(positions)
    displacement = np.vstack(motion)
    step = np.linalg.norm(np.diff(position, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(step)])
    return DeflectionProfile(s=s, position=position, displacement=displacement)


def _rigid_motion(element, U_local: np.ndarray) -> np.ndarray:
    """First six DOFs of each of the two nodes, node-major."""
    U_local = np.asarray(U_local, dtype=np.float64).ravel()
    per_node = int(getattr(element, "dof_per_node", 6))
    if per_node < 6:
        per_node = 6
    if U_local.size >= 2 * per_node:
        return np.concatenate([U_local[:6], U_local[per_node:per_node + 6]])
    if U_local.size >= 12:
        return U_local[:12]
    raise ValueError(f"element displacement has length {U_local.size}, need 12 rigid DOFs")


def _fibre_strain(element, U_local: np.ndarray, xi: float, y: float, z: float) -> float:
    """Axial strain at section point ``(y, z)`` from the rigid-section motion."""
    Ue = _rigid_motion(element, U_local)
    length = float(element.L)
    kind = interpolation_kind(type(element).__name__)
    def field(t: float) -> np.ndarray:
        return beam_field_interpolation(t, length, kind) @ Ue

    rate = (field(xi + _H) - field(xi - _H)) / (2.0 * _H) * (2.0 / length)
    return float(rate[0] + z * rate[4] - y * rate[5])


def _outer_points(outer_profile) -> np.ndarray:
    points = np.asarray(outer_profile, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 1:
        raise ValueError("outer_profile must have shape (n_points, 2) in section axes (y, z)")
    if not np.all(np.isfinite(points)):
        raise ValueError("outer_profile must be finite")
    return points


def extreme_fibre_profile(elements, U_global, outer_profile, *, xi: float = 0.0) -> ExtremeFibreProfile:
    """Extreme-fibre beam strain and stress along the outer profile of a static solution.

    One station per element, at natural coordinate ``xi`` (the element mid-point
    by default). ``outer_profile`` is the section-frame coordinates ``(y, z)`` of
    the outer contour, the same contour on every element. The strain uses the
    element's own displacement interpolation, in its section axes. The stress is
    ``E`` times that strain.

    Parameters
    ----------
    elements : sequence
        Solved beam elements, in mesh order.
    U_global : array
        Static global displacement, the same vector the element DOF map indexes.
    outer_profile : array (n_points, 2)
        Outer-profile points ``(y, z)`` in the section frame.
    xi : float
        Natural coordinate of the station on each element, in ``[-1, 1]``.
    """
    with stage("onedfem", "extreme_fibre"):
        return _extreme_fibre_profile(elements, U_global, outer_profile, xi=xi)


def _extreme_fibre_profile(elements, U_global, outer_profile, *, xi: float = 0.0) -> ExtremeFibreProfile:
    xi = float(xi)
    if not -1.0 <= xi <= 1.0:
        raise ValueError("xi must be in [-1, 1]")
    points = _outer_points(outer_profile)
    U = np.asarray(U_global, dtype=np.float64).ravel()
    strains = []
    stresses = []
    positions = []
    arc = []
    covered = 0.0
    t = 0.5 * (xi + 1.0)
    for element in elements:
        Ue = U[np.asarray(element.assemble_global_dof_indices(), dtype=int)]
        local, U_local, xi_local = local_state(element, Ue, xi)
        row = np.array([_fibre_strain(local, U_local, float(xi_local), float(y), float(z)) for y, z in points])
        modulus = float(local.E)
        strains.append(row)
        stresses.append(modulus * row)
        positions.append(axis_point(element, xi))
        ends = axis_point(element, 1.0) - axis_point(element, -1.0)
        length = float(np.linalg.norm(ends))
        arc.append(covered + t * length)
        covered += length
    strain = np.vstack(strains)
    stress = np.vstack(stresses)
    return ExtremeFibreProfile(
        s=np.asarray(arc, dtype=np.float64),
        position=np.vstack(positions),
        y=points[:, 0].copy(),
        z=points[:, 1].copy(),
        strain=strain,
        stress=stress,
        extreme_strain=np.column_stack([strain.min(axis=1), strain.max(axis=1)]),
        extreme_stress=np.column_stack([stress.min(axis=1), stress.max(axis=1)]),
    )


def _as_pair(solved, first_names: Sequence[str], second_name: str, label: str):
    if isinstance(solved, tuple) and len(solved) == 2:
        return solved[0], solved[1]
    first = None
    for name in first_names:
        first = getattr(solved, name, None)
        if first is not None:
            break
    if first is None:
        secondary = getattr(solved, "secondary_results", None)
        if isinstance(secondary, dict):
            first = secondary.get("global", {}).get(first_names[-1])
    second = getattr(solved, second_name, None)
    if first is None or second is None:
        raise RuntimeError(f"no {label} on this solve; run that analysis first")
    return first, second


def buckling_result(solved: Union[object, tuple]) -> BucklingResult:
    """Buckling load factors and modes a linear buckling solve already returned.

    ``solved`` is a buckling runner after ``run``, or the ``(load_factors,
    mode_shapes)`` pair from ``solve_linear_buckling_eigenpairs`` or
    ``solve_buckling_linearized_problem``. Nothing is solved here.
    """
    if isinstance(solved, BucklingResult):
        return solved
    factors, shapes = _as_pair(
        solved,
        ("buckling_load_factors",),
        "buckling_mode_shapes",
        "buckling result",
    )
    return BucklingResult(
        load_factors=np.asarray(factors, dtype=np.float64).ravel(),
        mode_shapes=np.asarray(shapes, dtype=np.float64),
    )


def blade_vibration_modes(solved: Union[object, tuple]) -> BladeVibrationModes:
    """Blade vibration modes an eigen solve already returned.

    ``solved`` is an eigen runner after ``run``, or the ``(frequencies_hz,
    mode_shapes)`` pair from ``solve_vibration_eigenpairs``. The modes are the
    beam's. A rotor is not assembled, and the eigenproblem is not solved again.
    """
    if isinstance(solved, BladeVibrationModes):
        return solved
    frequencies, shapes = _as_pair(
        solved,
        ("vibration_frequencies_hz", "frequencies_hz"),
        "vibration_mode_shapes",
        "blade vibration modes",
    )
    return BladeVibrationModes(
        frequencies_hz=np.atleast_1d(np.asarray(frequencies, dtype=np.float64).ravel()),
        mode_shapes=np.asarray(shapes, dtype=np.float64),
    )
