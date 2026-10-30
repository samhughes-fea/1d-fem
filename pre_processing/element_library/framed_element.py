# pre_processing/element_library/framed_element.py
"""Element frames for the linear beams: rotation to global axes, and curved, twisted elements.

``Element1DBase`` calls ``setup_frame`` while it builds an element, and wraps the
element kernels (``element_stiffness_matrix``, ``element_force_vector``,
``element_mass_matrix``, ``linear_geometric_stiffness_matrix``) of every subclass with
``frame_kernel``. Elements whose axes are the global axes (along ``+X``, no twist, no
curvature) run their kernels unchanged.

Straight, rotated element
    The kernel runs in the element's own axes and the result is rotated:
    ``K = T^T K_local T``, ``F = T^T F_local``, ``K_sigma(U) = T^T K_sigma_local(T U) T``,
    with ``T`` from ``element_frame.node_dof_rotation``.

Curved, twisted element (``precurvature.txt``)
    The arc is ``m`` straight sub-elements, each in its own frame along the arc
    (``ArcGeometry.sub_frames``), built from the element's own type and section with a
    sub-element length. The chain is assembled in global axes and its inner nodes are
    condensed out, so ``K`` is the stiffness of the arc between its end nodes. The mass
    and ``K_sigma`` use the same condensation (the static shapes of the chain). A twist
    rate ``k_x0`` also adds the extension-twist coupling of the helical fibres to the
    section (``SectionAnalysis.with_pretwist``).

Loads on a mesh with rotated or curved elements
    Point loads must sit on nodes; they are global components and are added once, to the
    lowest-numbered element at the node. Distributed loads are global components per
    length of the beam axis, sampled along global ``X`` (``distributed_load.txt``); each
    element spreads them along its own geometry (straight, or along the arc's
    sub-elements) with its work-conjugate fields. Wall loads are section loads, built in
    the element's (or each sub-element's) own axes. A precurvature makes no equivalent
    load: the curvature is in the geometry.

Public API
----------
setup_frame, frame_kernel, FRAMED_KERNELS, DEFAULT_SUB_ELEMENTS.
"""

from __future__ import annotations

import functools
from typing import Any, Dict

import numpy as np

from pre_processing.element_library.beam_load_interpolation import beam_field_interpolation
from pre_processing.element_library.element_frame import arc_geometry, node_dof_rotation, straight_frame

FRAMED_KERNELS = (
    "element_stiffness_matrix",
    "element_force_vector",
    "element_mass_matrix",
    "linear_geometric_stiffness_matrix",
)
DEFAULT_SUB_ELEMENTS = 8
_TOL = 1e-12


def _supports_frame(type_name: str) -> bool:
    """Whether the element type takes rotated and curved geometry (else ``precurvature.txt``
    keeps its older meaning for it: a reference strain on a straight element)."""
    import importlib

    from pre_processing.element_library.element_factory import ElementFactory

    module = ElementFactory.ELEMENT_CLASS_MAP.get(str(type_name))
    if module is None:
        return False
    cls = getattr(importlib.import_module(module), str(type_name), None)
    return bool(getattr(cls, "supports_element_frame", False))


def _element_frames(element_dictionary: Dict[str, Any], grid_dictionary: Dict[str, Any], idx: int):
    """``(R0, R1, arc or None)`` of element ``idx``."""
    n0, n1 = element_dictionary["connectivity"][idx]
    xyz = np.asarray(grid_dictionary["coordinates"], dtype=np.float64)
    twist = element_dictionary.get("twist")
    beta = 0.0 if twist is None else float(np.asarray(twist, dtype=np.float64)[idx])
    prec = element_dictionary.get("precurvature_per_element")
    k = np.zeros(3) if prec is None else np.asarray(prec, dtype=np.float64)[idx]
    if np.any(np.abs(k) > 0.0) and _supports_frame(np.asarray(element_dictionary["types"])[idx]):
        arc = arc_geometry(xyz[n0], xyz[n1], k, beta)
        return arc.R0, arc.R1, arc
    R = straight_frame(xyz[n0], xyz[n1], beta)
    return R, R, None


def _is_identity(R0: np.ndarray, R1: np.ndarray) -> bool:
    I = np.eye(3)
    return bool(np.allclose(R0, I, atol=_TOL) and np.allclose(R1, I, atol=_TOL))


def _mesh_cache(element_dictionary: Dict[str, Any], grid_dictionary: Dict[str, Any]) -> Dict[str, Any]:
    """Whether any element is rotated or curved, and the element that takes each node's point loads."""
    cache = element_dictionary.get("_frame_cache")
    xyz = np.asarray(grid_dictionary["coordinates"], dtype=np.float64)
    if cache is not None and cache.get("n_nodes") == xyz.shape[0]:
        return cache
    conn = np.asarray(element_dictionary["connectivity"])
    framed = False
    for i in range(conn.shape[0]):
        R0, R1, arc = _element_frames(element_dictionary, grid_dictionary, i)
        if arc is not None or not _is_identity(R0, R1):
            framed = True
            break
    owner = np.full(xyz.shape[0], -1, dtype=np.int64)
    for i in range(conn.shape[0]):            # lowest-numbered element at each node
        for n in conn[i]:
            if owner[n] < 0:
                owner[n] = i
    cache = {"framed": framed, "owner": owner, "n_nodes": xyz.shape[0], "node_frames": None}
    try:
        element_dictionary["_frame_cache"] = cache
    except TypeError:
        pass
    return cache


def _node_frames(element_dictionary: Dict[str, Any], grid_dictionary: Dict[str, Any]) -> np.ndarray:
    """(n_nodes, 3, 3): each node's own frame, the frame there of the lowest-numbered element
    at the node. Node vectors held in section components (the slopes of the third-order
    beams) are written in it."""
    cache = _mesh_cache(element_dictionary, grid_dictionary)
    if cache.get("node_frames") is None:
        conn = np.asarray(element_dictionary["connectivity"])
        frames = np.tile(np.eye(3), (cache["n_nodes"], 1, 1))
        done = set()
        for i in range(conn.shape[0]):
            R0, R1, _ = _element_frames(element_dictionary, grid_dictionary, i)
            for n, R in zip(conn[i], (R0, R1)):
                if int(n) not in done:
                    frames[n] = R
                    done.add(int(n))
        cache["node_frames"] = frames
    return cache["node_frames"]


def _slope_T(element, R0, R1, N0, N1, n_dof: int) -> np.ndarray:
    """``node_dof_rotation`` with the slopes of a third-order beam held as shear strains.

    The element's slope DOFs ``[v', w']`` (``frame_slope_dofs``) are the slopes along its
    own axis. Across a kink or twist the slopes of two elements differ, so the node holds
    the section shear strains ``[gamma_xy, gamma_xz]`` in its own frame ``N`` instead and
    each element recovers its slopes as ``v' = gamma_xy + theta_z`` and
    ``w' = gamma_xz - theta_y`` (its own components of the node rotation). On an element
    whose frame is the node's, ``gamma`` is the shear strain of the element there.
    """
    slopes = getattr(type(element), "frame_slope_dofs", ())
    T = node_dof_rotation(R0, R1, n_dof, pairs=[(st, nd, N0 if nd == 0 else N1) for st, nd in slopes])
    for st, nd in slopes:
        R = R0 if nd == 0 else R1
        rot = 3 if nd == 0 else 9
        T[st, rot:rot + 3] = R[:, 2]        # + theta_z (element axes)
        T[st + 1, rot:rot + 3] = -R[:, 1]   # - theta_y
    return T


def setup_frame(element, element_dictionary, grid_dictionary, idx: int, ctor_kwargs: Dict[str, Any]) -> None:
    """Set the element's frame attributes and, on a framed mesh, its nodal point loads."""
    R0, R1, arc = _element_frames(element_dictionary, grid_dictionary, idx)
    element._frame_R0, element._frame_R1, element._frame_arc = R0, R1, arc
    element._frame_N0 = element._frame_N1 = None
    element._frame_on = arc is not None or not _is_identity(R0, R1)
    element._frame_mesh = _mesh_cache(element_dictionary, grid_dictionary)["framed"]
    element._frame_nodal_load = None
    element._frame_ctor = None
    element._frame_depth = 0
    if not element._frame_mesh:
        return
    if not getattr(type(element), "supports_element_frame", False):
        raise NotImplementedError(
            f"Element {element.element_id} ({type(element).__name__}): rotated or curved elements "
            "(a beam off the global x axis, a [twist] column or precurvature.txt) are supported by the "
            "linear beams only: Euler-Bernoulli, Timoshenko and Benscoter, GBT, Reddy and Levinson; and by the "
            "corotational Euler-Bernoulli, Timoshenko and Benscoter beams.")
    # distributed loads: global components per length of the beam axis, sampled along global X;
    # spread here along the element's own geometry (the kernels would read them in local x)
    dl = np.asarray(element.distributed_load_array)
    element._frame_distributed = None
    if dl.size and np.any(np.asarray(dl, dtype=np.float64)[:, 3:9] != 0.0):
        element._frame_distributed = _load_table(np.asarray(dl, dtype=np.float64))
    if dl.size:
        element.distributed_load_array = np.zeros((0, dl.shape[1]))
    element._frame_p0 = np.asarray(grid_dictionary["coordinates"], dtype=np.float64)[
        element_dictionary["connectivity"][idx][0]].copy()
    if getattr(type(element), "frame_slope_dofs", ()):
        nf = _node_frames(element_dictionary, grid_dictionary)
        n0, n1 = element_dictionary["connectivity"][idx]
        element._frame_N0, element._frame_N1 = nf[n0], nf[n1]
    # wall loads are section (local) loads: the kernel builds them in the element's axes
    # point loads at this element's nodes, if it owns them (global components)
    cache = _mesh_cache(element_dictionary, grid_dictionary)
    xyz = np.asarray(grid_dictionary["coordinates"], dtype=np.float64)
    nodes = list(element_dictionary["connectivity"][idx])
    load = np.zeros(12, dtype=np.complex128)
    pl = np.asarray(element.point_load_array, dtype=np.float64)
    scale = max(float(np.ptp(xyz, axis=0).max()), 1.0)
    for row in (pl if pl.size else np.zeros((0, 9))):
        d = np.linalg.norm(xyz - row[:3], axis=1)
        n = int(np.argmin(d))
        if d[n] > 1e-6 * scale:
            raise ValueError(
                f"point load at {row[:3].tolist()} is not at a node; on a mesh with rotated or curved "
                "elements point loads must sit on nodes")
        if cache["owner"][n] != idx:
            continue
        phase = float(row[9]) if row.size > 9 else 0.0
        a = nodes.index(n)
        load[6 * a:6 * a + 6] += row[3:9] * np.exp(1j * phase)
    element._frame_nodal_load = load.real if not np.any(load.imag) else load
    element.point_load_array = np.zeros((0, pl.shape[1] if pl.ndim == 2 and pl.shape[1] else 9))
    element._E_0_voigt = np.zeros_like(element._E_0_voigt)   # the curvature is in the geometry now
    if arc is not None:
        element._frame_ctor = ctor_kwargs


# --------------------------------------------------------------- line loads on framed elements
def _bending_kind(element) -> str:
    """Field interpolation the element's loads are spread with (``beam_load_interpolation``)."""
    return "hermite" if "EulerBernoulli" in type(element).__name__ else "linear"


def _load_table(dl: np.ndarray):
    """``q(X)``: the six global load components, linear in the global ``X`` of the samples."""
    x = np.asarray(dl[:, 0], dtype=np.float64)
    order = np.argsort(x, kind="stable")
    xs, q = x[order], np.asarray(dl[order, 3:9])
    span = max(float(xs[-1] - xs[0]), 1.0)

    def table(X):
        X = np.atleast_1d(np.asarray(X, dtype=np.float64))
        if np.any(X < xs[0] - 1e-9 * span) or np.any(X > xs[-1] + 1e-9 * span):
            raise ValueError(
                f"distributed load samples cover X = {xs[0]:g} to {xs[-1]:g}; the element reaches "
                f"X = {X.min():g} to {X.max():g}")
        if xs.size == 1:
            return np.repeat(q[:1], X.size, axis=0)
        return np.column_stack([np.interp(X, xs, q[:, j]) for j in range(6)])

    return table


def _field(element, L: float):
    """``N(xi)`` (6, n): the element's work-conjugate right-handed fields over a length ``L``."""
    own = getattr(type(element), "frame_load_fields", None)
    if own is not None:
        return lambda xi: own(element, xi)
    kind = _bending_kind(element)
    return lambda xi: beam_field_interpolation(xi, L, kind)


def _segment_load(table, p0: np.ndarray, R: np.ndarray, L: float, N_of, T: np.ndarray, n_gauss: int = 6) -> np.ndarray:
    """Global nodal loads of the line load ``table`` on the straight segment from ``p0`` along ``R[:, 0]``.

    ``N_of(xi)`` (6, n) are the segment's fields in its own axes and ``T`` (n, n) its DOF
    rotation; the result has the ``n`` DOFs of ``N_of``.
    """
    xg, wg = np.polynomial.legendre.leggauss(n_gauss)
    pts = p0[None, :] + (0.5 * (xg + 1.0) * L)[:, None] * R[:, 0][None, :]
    q = table(pts[:, 0])
    q_local = np.hstack([q[:, :3] @ R, q[:, 3:] @ R])        # R^T q per point
    F_local = sum(N_of(x).T @ qg * (w * L / 2) for x, w, qg in zip(xg, wg, q_local))
    return T[:F_local.size, :F_local.size].T @ F_local


# ---------------------------------------------------------------------- kernels
def _T(element, n_dof: int, R0=None, R1=None) -> np.ndarray:
    R0 = element._frame_R0 if R0 is None else R0
    R1 = element._frame_R1 if R1 is None else R1
    if element._frame_N0 is not None:
        return _slope_T(element, R0, R1, element._frame_N0, element._frame_N1, n_dof)
    return node_dof_rotation(R0, R1, n_dof)


def _with_K(obj, K):
    if hasattr(obj, "K_e"):
        obj.K_e = K
        return obj
    return K


def frame_kernel(name: str, fn):
    """Wrap a kernel so that it returns global-axis results for a rotated or curved element."""

    @functools.wraps(fn)
    def wrapped(self, *args, **kwargs):
        if not getattr(self, "_frame_mesh", False) or self.__dict__.get("_frame_depth", 0):
            return fn(self, *args, **kwargs)
        self._frame_depth = 1
        try:
            if self._frame_arc is not None:
                return _curved(self, name, fn, *args, **kwargs)
            return _straight(self, name, fn, *args, **kwargs)
        finally:
            self._frame_depth = 0

    wrapped._frame_wrapped = True
    return wrapped


def _straight(self, name, fn, *args, **kwargs):
    if name == "linear_geometric_stiffness_matrix":
        U = np.asarray(args[0], dtype=np.float64)
        T = _T(self, U.size)
        return T.T @ fn(self, T @ U, *args[1:], **kwargs) @ T
    out = fn(self, *args, **kwargs)
    if name == "element_stiffness_matrix":
        T = _T(self, out.K_e.shape[0])
        out.K_e = T.T @ out.K_e @ T
        _rotate_gauss_data(out, T)
    elif name == "element_mass_matrix":
        T = _T(self, out.M_e.shape[0])
        out.M_e = T.T @ out.M_e @ T
    elif name == "element_force_vector":
        T = _T(self, out.F_e.shape[0])
        F = T.T @ out.F_e
        if self._frame_distributed is not None:
            L = float(self.grid_array[1, 0])
            Fd = _segment_load(self._frame_distributed, self._frame_p0, self._frame_R0, L, _field(self, L), T)
            F[:Fd.size] += Fd
        if self._frame_nodal_load is not None:
            F = F.astype(np.result_type(F, self._frame_nodal_load))
            F[:12] += self._frame_nodal_load
        out.F_e = F
        _rotate_gauss_data(out, T)
    return out


def _rotate_gauss_data(obj, T: np.ndarray) -> None:
    """Make cached Gauss-point operators act on global DOFs: ``B T`` and ``T^T N``."""
    n = T.shape[0]
    for gp in getattr(obj, "gauss_data", None) or []:
        for attr in ("B_matrix", "shape_functions", "shape_derivatives"):
            X = getattr(gp, attr, None)
            if X is None:
                continue
            X = np.asarray(X)
            if X.ndim >= 2 and X.shape[-1] == n:
                setattr(gp, attr, X @ T)
            elif X.ndim >= 2 and X.shape[-2] == n:
                setattr(gp, attr, np.swapaxes(np.swapaxes(X, -1, -2) @ T, -1, -2))


# --------------------------------------------------------------- curved elements
class _Chain:
    """The straight sub-elements of a curved element and the condensation of their inner nodes."""

    def __init__(self, element):
        from pre_processing.element_library.element_1D_base import Element1DBase

        arc = element._frame_arc
        m = int(getattr(element, "frame_sub_elements", DEFAULT_SUB_ELEMENTS))
        self.m = m
        kw = dict(element._frame_ctor)
        ed = dict(kw["element_dictionary"])
        idx = int(np.where(np.asarray(ed["ids"]) == element.element_id)[0][0])
        prec = np.array(ed["precurvature_per_element"], dtype=np.float64, copy=True)
        k_twist = float(prec[idx, 0])
        prec[idx] = 0.0
        ed["precurvature_per_element"] = prec
        if "twist" in ed:
            tw = np.array(ed["twist"], dtype=np.float64, copy=True)
            tw[idx] = 0.0
            ed["twist"] = tw
        pt = np.zeros(len(ed["ids"]))
        pt[idx] = k_twist
        ed["pretwist_rate"] = pt
        ed.pop("_frame_cache", None)
        gd = dict(kw["grid_dictionary"])
        xyz = np.array(gd["coordinates"], dtype=np.float64, copy=True)
        n0, n1 = ed["connectivity"][idx]
        h = arc.sub_chord(m)
        xyz[n0] = 0.0                         # only the element's own nodes are read
        xyz[n1] = [h, 0.0, 0.0]
        gd["coordinates"] = xyz
        ed["_frame_cache"] = {"framed": False, "owner": np.zeros(xyz.shape[0], dtype=np.int64),
                              "n_nodes": xyz.shape[0]}
        kw.update(element_dictionary=ed, grid_dictionary=gd,
                  point_load_array=np.zeros((0, 9)), distributed_load_array=np.zeros((0, 9)))
        self.sub: Element1DBase = type(element)(**kw)
        self.frames = arc.sub_frames(m)
        n_dof = int(element._n_dof if hasattr(element, "_n_dof") else 2 * element.dof_per_node)
        self.n_dof = n_dof
        self.dpn = n_dof // 2
        self.n_ext = self.dpn - 6
        self.K_sub = np.asarray(self.sub.element_stiffness_matrix().K_e, dtype=np.float64)
        slopes = getattr(type(element), "frame_slope_dofs", ())
        if slopes:
            # chain node frames: the element's node frames at its ends, the sub-element frames inside
            nodes = [element._frame_N0] + list(self.frames[1:]) + [element._frame_N1]
            self.T = [_slope_T(element, R, R, nodes[i], nodes[i + 1], n_dof) for i, R in enumerate(self.frames)]
        else:
            self.T = [node_dof_rotation(R, R, n_dof) for R in self.frames]
        self.maps = [self._map(i) for i in range(m)]
        n = (m + 1) * self.dpn
        K = np.zeros((n, n))
        for T, g in zip(self.T, self.maps):
            K[np.ix_(g, g)] += T.T @ self.K_sub @ T
        b = self._boundary()
        i = np.setdiff1d(np.arange(n), b)
        self.b, self.i = b, i
        Phi = np.zeros((n, b.size))
        Phi[b, np.arange(b.size)] = 1.0
        if i.size:
            Phi[i] = -np.linalg.solve(K[np.ix_(i, i)], K[np.ix_(i, b)])
        self.Phi = Phi
        self.K = Phi.T @ K @ Phi
        self.K = 0.5 * (self.K + self.K.T)

    def _node_dofs(self, node: int):
        base = node * self.dpn
        return list(range(base, base + 6)), list(range(base + 6, base + self.dpn))

    def _map(self, i: int) -> np.ndarray:
        """Chain DOFs of sub-element ``i`` in element DOF order."""
        s0, e0 = self._node_dofs(i)
        s1, e1 = self._node_dofs(i + 1)
        return np.array(s0 + s1 + e0 + e1)

    def _boundary(self) -> np.ndarray:
        """Chain DOFs of the curved element's own DOFs, in its DOF order."""
        s0, e0 = self._node_dofs(0)
        s1, e1 = self._node_dofs(self.m)
        return np.array(s0 + s1 + e0 + e1)

    def assemble(self, sub_matrices) -> np.ndarray:
        n = (self.m + 1) * self.dpn
        K = np.zeros((n, n))
        for Ks, T, g in zip(sub_matrices, self.T, self.maps):
            K[np.ix_(g, g)] += T.T @ Ks @ T
        return self.Phi.T @ K @ self.Phi

    def line_load(self, element) -> np.ndarray:
        """Nodal loads of the element's distributed and wall loads, condensed from the sub-elements."""
        F = np.zeros((self.m + 1) * self.dpn)
        any_load = False
        F_wall = np.asarray(self.sub.element_force_vector().F_e, dtype=np.float64)   # wall loads, local
        if np.any(F_wall != 0.0):
            any_load = True
            for T, g in zip(self.T, self.maps):
                F[g] += T.T @ F_wall
        table = element._frame_distributed
        if table is not None:
            any_load = True
            arc = element._frame_arc
            h = arc.length / self.m
            p_mid = element._frame_p0 - (arc.position(-0.5 * arc.length, np.zeros(3)))
            chord = arc.sub_chord(self.m)
            N_of = _field(self.sub, chord)
            for i, (R, T, g) in enumerate(zip(self.frames, self.T, self.maps)):
                p0 = arc.position(-0.5 * arc.length + i * h, p_mid)
                # per length of the arc: each chord carries its arc's share
                Fd = _segment_load(table, p0, R, chord, N_of, T) * (h / chord)
                F[g[:Fd.size]] += Fd
        return self.Phi.T @ F if any_load else np.zeros(self.n_dof)

    def sub_displacements(self, U_e: np.ndarray):
        """Each sub-element's displacements in its own axes."""
        U = self.Phi @ np.asarray(U_e, dtype=np.float64)
        return [T @ U[g] for T, g in zip(self.T, self.maps)]


def chain(element) -> _Chain:
    c = element.__dict__.get("_frame_chain")
    if c is None:
        c = _Chain(element)
        element._frame_chain = c
    return c


def _curved(self, name, fn, *args, **kwargs):
    from pre_processing.element_library.gauss_point_data import ElementObject, ForceObject, MassObject

    c = chain(self)
    if name == "element_stiffness_matrix":
        return ElementObject(element_id=self.element_id, element_type=type(self).__name__, K_e=c.K.copy(),
                             gauss_data=[], integration_scheme=f"curved: {c.m} condensed sub-elements")
    if name == "element_force_vector":
        F = np.zeros(c.n_dof, dtype=np.result_type(float, self._frame_nodal_load))
        F[:12] += self._frame_nodal_load
        F += c.line_load(self)
        return ForceObject(element_id=self.element_id, element_type=type(self).__name__, F_e=F, gauss_data=[])
    if name == "element_mass_matrix":
        M_sub = np.asarray(c.sub.element_mass_matrix().M_e, dtype=np.float64)
        return MassObject(element_id=self.element_id, element_type=type(self).__name__,
                          M_e=c.assemble([M_sub] * c.m))
    if name == "linear_geometric_stiffness_matrix":
        subs = c.sub_displacements(args[0])
        return c.assemble([c.sub.linear_geometric_stiffness_matrix(u) for u in subs])
    raise AssertionError(name)


def local_state(element, U_e, xi: float):
    """``(element, U_local, xi_local)`` for recovering the section state at ``xi`` of ``element``.

    For a curved element, the straight sub-element at ``xi`` with its own displacements.
    """
    U_e = np.asarray(U_e, dtype=np.float64)
    if not getattr(element, "_frame_mesh", False):
        return element, U_e, xi
    if element._frame_arc is None:
        return element, _T(element, U_e.size) @ U_e, xi
    c = chain(element)
    t = 0.5 * (float(xi) + 1.0) * c.m
    i = min(int(np.floor(t)), c.m - 1)
    return c.sub, c.sub_displacements(U_e)[i], 2.0 * (t - i) - 1.0


__all__ = ["setup_frame", "frame_kernel", "FRAMED_KERNELS", "DEFAULT_SUB_ELEMENTS", "local_state", "chain"]
