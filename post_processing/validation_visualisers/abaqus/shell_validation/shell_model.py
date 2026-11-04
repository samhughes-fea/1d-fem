# post_processing/validation_visualisers/abaqus/shell_validation/shell_model.py
"""Shell models of prismatic thin-walled beams, written as Abaqus input files.

The section is the same description that 1d-fem's thin-walled section analysis uses:
points ``{id: (y, z)}``, walls ``(wall_id, start, end, laminate_id)`` and laminates as
plies ``(material, angle_deg, thickness)`` listed from the bottom of the wall (``n = -h/2``)
to the top. The beam runs along ``x`` and every wall mid-line is swept along it into a
strip of shell elements. Walls that share a section point share its nodes, so the shells
meet on their mid-surfaces, as the section analysis assumes.

Element nodes are ordered ``(x_i, s_a), (x_i+1, s_a), (x_i+1, s_b), (x_i, s_b)`` with ``s``
running from the wall's start point to its end point. The shell normal is then
``x cross s``, which is the wall normal ``n`` of the section analysis, and each ply's
orientation puts its fibre at ``angle`` from ``x`` toward ``s``. So a laminate means the
same thing in both models.

Supports and loads are given per end section:

- ``clamp``: every node of the section is fixed (warping restrained).
- ``fork``: a reference node is tied to the section by a kinematic coupling of the
  in-plane DOFs 2 and 3 only, so the section stays rigid in its plane while it warps
  freely. The reference node's U2, U3 and UR1 are fixed.
- ``diaphragm``: the same coupling, with the reference node free, to apply shear forces
  and a torque without restraining the warping.
- ``rigid``: a kinematic coupling of all six DOFs to the reference node.
- ``free``: nothing.

Axial forces and bending moments at an end are applied as nodal forces with the linear
stress distribution of a plane section, weighted by each node's share of the membrane
axial stiffness ``t ds / a11`` (``a11`` from the inverse membrane stiffness of its walls),
so they do not restrain the warping either.

The file has no Abaqus CAE dependency: run it with ``abaqus job=<name> input=<name>.inp``.
With ``dialect="calculix"`` the same model is written for CalculiX (output requests
differ), which is how the models are smoke-tested without Abaqus.

Public API
----------
Ply, Laminate, ShellSection, ShellBeamModel, EndCondition, NodalLoad, Step, write_inp.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class Material:
    """Orthotropic ply material in its own axes (1 = fibre)."""

    name: str
    E1: float
    E2: float
    G12: float
    nu12: float
    G13: float
    G23: float
    rho: float = 0.0

    def engineering_constants(self) -> Tuple[float, ...]:
        """``E1, E2, E3, nu12, nu13, nu23, G12, G13, G23`` with transverse isotropy in 2-3."""
        nu23 = float(np.clip(self.E2 / (2.0 * self.G23) - 1.0, 0.0, 0.49))
        return (self.E1, self.E2, self.E2, self.nu12, self.nu12, nu23, self.G12, self.G13, self.G23)

    def reduced_stiffness(self) -> np.ndarray:
        nu21 = self.nu12 * self.E2 / self.E1
        d = 1.0 - self.nu12 * nu21
        return np.array([[self.E1 / d, self.nu12 * self.E2 / d, 0.0],
                         [self.nu12 * self.E2 / d, self.E2 / d, 0.0],
                         [0.0, 0.0, self.G12]])


@dataclass(frozen=True)
class Ply:
    material: Material
    angle_deg: float
    thickness: float


@dataclass(frozen=True)
class Laminate:
    """Plies from the bottom (``n = -h/2``) to the top of the wall."""

    plies: Tuple[Ply, ...]

    @property
    def thickness(self) -> float:
        return float(sum(p.thickness for p in self.plies))

    def membrane_stiffness(self) -> np.ndarray:
        """``A`` (3, 3) in wall axes ``[eps_xx, eps_ss, gamma_xs]``."""
        A = np.zeros((3, 3))
        for p in self.plies:
            c, s = np.cos(np.deg2rad(p.angle_deg)), np.sin(np.deg2rad(p.angle_deg))
            T = np.array([[c * c, s * s, c * s], [s * s, c * c, -c * s], [-2 * c * s, 2 * c * s, c * c - s * s]])
            A += T.T @ p.material.reduced_stiffness() @ T * p.thickness
        return A

    def axial_stiffness(self) -> float:
        """Membrane axial stiffness per unit width with ``N_ss = N_xs = 0``: ``1 / inv(A)[0, 0]``."""
        return float(1.0 / np.linalg.inv(self.membrane_stiffness())[0, 0])


@dataclass
class ShellSection:
    """A thin-walled section: points, walls and laminates, with elements per wall."""

    points: Dict[int, Tuple[float, float]]
    walls: List[Tuple[int, int, int, int]]          # (wall_id, start point, end point, laminate_id)
    laminates: Dict[int, Laminate]
    element_size: float                             # target element size along the walls, m
    min_elements_per_wall: int = 2

    def wall_divisions(self) -> Dict[int, int]:
        out = {}
        for wid, a, b, _ in self.walls:
            length = float(np.hypot(*(np.subtract(self.points[b], self.points[a]))))
            out[wid] = max(self.min_elements_per_wall, int(np.ceil(length / self.element_size - 1e-9)))
        return out


@dataclass(frozen=True)
class EndCondition:
    """Support or load introduction at ``x = 0`` (``"root"``) or ``x = L`` (``"tip"``)."""

    kind: str = "free"                               # clamp | fork | diaphragm | rigid | free
    reference: Tuple[float, float] = (0.0, 0.0)      # (y, z) of the reference node
    fix_axial: bool = False                          # fork: also fix U1 of the section node nearest the reference


@dataclass(frozen=True)
class NodalLoad:
    """A load per unit load factor: an end resultant or a force on a set of section nodes.

    ``where`` is ``"root"`` or ``"tip"`` for an end load, or an ``x`` coordinate. ``resultant``
    is ``[F_x, F_y, F_z, M_x, M_y, M_z]`` about the section origin. At an end, ``F_x``,
    ``M_y`` and ``M_z`` go on the section nodes as a plane-section stress distribution and
    ``F_y``, ``F_z``, ``M_x`` on the end's reference node (which needs a ``fork``,
    ``diaphragm`` or ``rigid`` end for them). At an interior ``x`` only ``F_y`` and ``F_z``
    are allowed, spread evenly over the section nodes with ``z`` (for ``F_y``) or ``y`` (for
    ``F_z``) equal to ``at`` within the mesh tolerance.
    """

    where: object
    resultant: Tuple[float, float, float, float, float, float]
    at: Optional[float] = None


@dataclass(frozen=True)
class Step:
    """``kind`` is ``"buckle"`` (``n_modes``) or ``"static"`` (``nlgeom``, ``load_points``)."""

    kind: str
    n_modes: int = 10
    nlgeom: bool = True
    load_points: Tuple[float, ...] = (1.0,)          # load factors with output
    max_increment: float = 0.05
    min_increment: float = 1e-6


@dataclass
class ShellBeamModel:
    name: str
    section: ShellSection
    length: float
    n_x: int                                          # elements along the span
    root: EndCondition
    tip: EndCondition
    loads: List[NodalLoad]
    step: Step
    stations: Tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)   # x / L with displacement output
    element_type: str = "S4R"
    description: str = ""

    # -------------------------------------------------------------- mesh
    def mesh(self) -> "ShellMesh":
        return ShellMesh.build(self)


@dataclass
class ShellMesh:
    """Nodes on a grid of section nodes ``j`` and stations ``k``: label ``k * n_sec + j + 1``.

    Quadratic elements (``S8R``) add a section node in the middle of every segment and a
    station between every pair, and skip the grid nodes in the middle of each element.
    """

    section_yz: np.ndarray                            # (n_sec, 2) section node coordinates
    section_weight: np.ndarray                        # (n_sec,) axial stiffness share, N
    midside: np.ndarray                               # (n_sec,) True for a segment's middle node
    segments: List[Tuple[int, int, int, int]]         # (node a, middle node or -1, node b, wall_id)
    x: np.ndarray                                     # station coordinates
    order: int = 1
    elements: Dict[int, List[Tuple[int, ...]]] = field(default_factory=dict)   # wall_id -> node tuples

    @property
    def n_sec(self) -> int:
        return self.section_yz.shape[0]

    @property
    def last(self) -> int:
        """Station index of the tip."""
        return self.x.size - 1

    def label(self, k: int, j: int) -> int:
        """Node label of section node ``j`` at station ``k`` (1-based)."""
        return k * self.n_sec + j + 1

    def exists(self, k: int, j: int) -> bool:
        return not (self.order == 2 and k % 2 == 1 and self.midside[j])

    def section_nodes(self, k: int) -> List[int]:
        """Section-node indices ``j`` with a node at station ``k``."""
        return [j for j in range(self.n_sec) if self.exists(k, j)]

    @classmethod
    def build(cls, model: ShellBeamModel) -> "ShellMesh":
        order = {"S4R": 1, "S4": 1, "S8R": 2}.get(model.element_type.upper())
        if order is None:
            raise ValueError(f"element type must be S4R, S4 or S8R, got {model.element_type!r}")
        sec = model.section
        divisions = sec.wall_divisions()
        yz: List[Tuple[float, float]] = []
        mid: List[bool] = []
        point_node: Dict[int, int] = {}

        def new(p, is_mid=False):
            yz.append(tuple(map(float, p)))
            mid.append(is_mid)
            return len(yz) - 1

        def point(pid):
            if pid not in point_node:
                point_node[pid] = new(sec.points[pid])
            return point_node[pid]

        weight: Dict[int, float] = {}

        def add(j, v):
            weight[j] = weight.get(j, 0.0) + v

        segments = []
        for wid, a, b, lid in sec.walls:
            n = divisions[wid]
            pa, pb = np.asarray(sec.points[a], float), np.asarray(sec.points[b], float)
            chain = [point(a)] + [new(pa + (pb - pa) * i / n) for i in range(1, n)] + [point(b)]
            ds = float(np.linalg.norm(pb - pa)) / n
            k_ax = sec.laminates[lid].axial_stiffness() * ds
            for i, (na, nb) in enumerate(zip(chain[:-1], chain[1:])):
                if order == 2:   # Simpson weights: the consistent loads of a quadratic edge
                    nm = new(pa + (pb - pa) * (i + 0.5) / n, True)
                    add(na, k_ax / 6), add(nm, 4 * k_ax / 6), add(nb, k_ax / 6)
                else:
                    nm = -1
                    add(na, k_ax / 2), add(nb, k_ax / 2)
                segments.append((na, nm, nb, wid))
        x = np.linspace(0.0, model.length, order * model.n_x + 1)
        mesh = cls(np.array(yz), np.array([weight[j] for j in range(len(yz))]), np.array(mid), segments, x, order)
        lab = mesh.label
        for na, nm, nb, wid in segments:
            elems = mesh.elements.setdefault(wid, [])
            for e in range(model.n_x):
                k = order * e
                if order == 1:
                    elems.append((lab(k, na), lab(k + 1, na), lab(k + 1, nb), lab(k, nb)))
                else:
                    elems.append((lab(k, na), lab(k + 2, na), lab(k + 2, nb), lab(k, nb),
                                  lab(k + 1, na), lab(k + 2, nm), lab(k + 1, nb), lab(k, nm)))
        return mesh

    def station_index(self, x_over_L: float) -> int:
        """The element-corner station nearest ``x / L``."""
        corners = np.arange(0, self.x.size, self.order)
        return int(corners[np.argmin(np.abs(self.x[corners] / self.x[-1] - x_over_L))])

    def plane_section_forces(self, Fx: float, My: float, Mz: float) -> np.ndarray:
        """Nodal ``f_i = w_i (a + b y_i + c z_i)`` with ``sum f = Fx``, ``sum f z = My``, ``-sum f y = Mz``."""
        y, z = self.section_yz.T
        w = self.section_weight
        G = np.array([[w.sum(), w @ y, w @ z], [w @ y, w @ y**2, w @ (y * z)], [w @ z, w @ (y * z), w @ z**2]])
        coef = np.linalg.solve(G, np.array([Fx, -Mz, My]))
        return w * (coef[0] + coef[1] * y + coef[2] * z)


# ------------------------------------------------------------------ writer
def _fmt(v: float) -> str:
    return f"{v:.10g}"


def _lines(values: Sequence, per_line: int = 8) -> str:
    vals = [str(v) for v in values]
    return "".join(", ".join(vals[i:i + per_line]) + "\n" for i in range(0, len(vals), per_line))


def write_inp(model: ShellBeamModel, path: str, dialect: str = "abaqus") -> ShellMesh:
    """Write ``model`` to ``path``; returns the mesh (for post-processing)."""
    if dialect not in ("abaqus", "calculix"):
        raise ValueError(f"dialect must be 'abaqus' or 'calculix', got {dialect!r}")
    mesh = model.mesh()
    sec = model.section
    out: List[str] = []
    w = out.append
    w("*Heading\n")
    w(f"** {model.name}: {model.description or 'thin-walled beam shell model'}\n")
    w("** Generated by 1d-fem post_processing/validation_visualisers/abaqus/shell_validation.\n")
    w("** Units: N, m. Beam axis x; section coordinates (y, z) as in 1d-fem.\n")

    # nodes
    w("*Node\n")
    for k, xk in enumerate(mesh.x):
        for j in mesh.section_nodes(k):
            y, z = mesh.section_yz[j]
            w(f"{mesh.label(k, j)}, {_fmt(xk)}, {_fmt(y)}, {_fmt(z)}\n")
    n_nodes = mesh.x.size * mesh.n_sec
    ref_labels: Dict[str, int] = {}
    for i, (end, cond) in enumerate((("root", model.root), ("tip", model.tip))):
        if cond.kind in ("fork", "diaphragm", "rigid"):
            ref_labels[end] = n_nodes + 1 + i
            xk = 0.0 if end == "root" else model.length
            w(f"{ref_labels[end]}, {_fmt(xk)}, {_fmt(cond.reference[0])}, {_fmt(cond.reference[1])}\n")
    # CalculiX: a rigid end is a *RIGID BODY, whose rotations live on a second node's translations
    rot_labels: Dict[str, int] = {}
    if dialect == "calculix":
        for i, (end, cond) in enumerate((("root", model.root), ("tip", model.tip))):
            if cond.kind == "rigid":
                rot_labels[end] = n_nodes + 3 + i
                xk = 0.0 if end == "root" else model.length
                w(f"{rot_labels[end]}, {_fmt(xk)}, {_fmt(cond.reference[0])}, {_fmt(cond.reference[1])}\n")

    # elements, one element set per wall
    eid = 0
    for wid, *_ in sec.walls:
        w(f"*Element, type={model.element_type}, elset=WALL{wid}\n")
        for conn in mesh.elements[wid]:
            eid += 1
            w(f"{eid}, " + ", ".join(map(str, conn)) + "\n")

    # materials
    mats = {p.material.name: p.material for lam in sec.laminates.values() for p in lam.plies}
    for name, m in mats.items():
        w(f"*Material, name={name}\n*Elastic, type=ENGINEERING CONSTANTS\n")
        c = m.engineering_constants()
        w(", ".join(_fmt(v) for v in c[:8]) + "\n" + _fmt(c[8]) + "\n")
        if m.rho > 0:
            w(f"*Density\n{_fmt(m.rho)}\n")

    # one orientation per (wall, ply angle): fibre at angle from x toward s
    for wid, a, b, lid in sec.walls:
        es = np.subtract(sec.points[b], sec.points[a]).astype(float)
        es /= np.linalg.norm(es)
        s3 = np.array([0.0, es[0], es[1]])
        x3 = np.array([1.0, 0.0, 0.0])
        for ang in sorted({p.angle_deg for p in sec.laminates[lid].plies}):
            c, s = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
            e1, e2 = c * x3 + s * s3, -s * x3 + c * s3
            w(f"*Orientation, name={_orientation_name(wid, ang)}, system=RECTANGULAR\n")
            w(", ".join(_fmt(v) for v in (*e1, *e2)) + "\n")
            if dialect == "abaqus":
                w("3, 0.\n")

    # shell sections (CalculiX takes a composite section on quadratic elements only)
    for wid, _, _, lid in sec.walls:
        plies = sec.laminates[lid].plies
        if dialect == "calculix" and len(plies) == 1:
            p = plies[0]
            w(f"*Shell Section, elset=WALL{wid}, material={p.material.name}, "
              f"orientation={_orientation_name(wid, p.angle_deg)}\n{_fmt(p.thickness)}\n")
            continue
        w(f"*Shell Section, elset=WALL{wid}, composite\n")
        for p in plies:
            w(f"{_fmt(p.thickness)}, 3, {p.material.name}, {_orientation_name(wid, p.angle_deg)}\n")

    # node sets: end sections and output stations
    for end, k in (("ROOT", 0), ("TIP", mesh.last)):
        w(f"*Nset, nset={end}\n" + _lines([mesh.label(k, j) for j in mesh.section_nodes(k)]))
    station_ids = sorted({mesh.station_index(s) for s in model.stations})
    for k in station_ids:
        w(f"*Nset, nset=ST{k}\n" + _lines([mesh.label(k, j) for j in mesh.section_nodes(k)]))
    w("*Nset, nset=STATIONS\n" + _lines([f"ST{k}" for k in station_ids]))
    if ref_labels:
        w("*Nset, nset=REFS\n" + _lines(sorted(ref_labels.values())))
        for end, lab in ref_labels.items():
            w(f"*Nset, nset=RP_{end.upper()}\n{lab}\n")

    # couplings
    for end, cond in (("root", model.root), ("tip", model.tip)):
        if end not in ref_labels:
            continue
        if end in rot_labels:
            w(f"*Rigid Body, nset={end.upper()}, ref node={ref_labels[end]}, rot node={rot_labels[end]}\n")
            continue
        w(f"*Surface, type=NODE, name=S_{end.upper()}\n{end.upper()}\n")
        w(f"*Coupling, constraint name=C_{end.upper()}, ref node={ref_labels[end]}, surface=S_{end.upper()}\n")
        w("*Kinematic\n" + ("1, 6\n" if cond.kind == "rigid" else "2, 3\n"))

    # boundary conditions
    w("*Boundary\n")
    for end, cond, k in (("root", model.root, 0), ("tip", model.tip, mesh.last)):
        if cond.kind == "clamp":
            w(f"{end.upper()}, 1, 6\n")
        if cond.kind in ("fork", "diaphragm"):
            # only the in-plane DOFs are coupled: the reference node's out-of-plane rotations carry nothing
            w(f"RP_{end.upper()}, 5, 6\n")
        if cond.kind == "fork":
            w(f"RP_{end.upper()}, 2, 4\n")
            if cond.fix_axial:
                j = int(np.argmin(np.hypot(*(mesh.section_yz - np.asarray(cond.reference)).T)))
                w(f"{mesh.label(k, j)}, 1, 1\n")

    # step
    st = model.step
    if st.kind == "static":
        w("*Time Points, name=LEVELS\n" + _lines([_fmt(t) for t in st.load_points]))
    abq = dialect == "abaqus"
    if st.kind == "buckle":
        w("*Step, name=BUCKLE, perturbation\n*Buckle\n" if abq else "*Step\n*Buckle\n")
        w(f"{st.n_modes}, , {max(2 * st.n_modes, st.n_modes + 8)}, 500\n")
    elif st.kind == "static":
        if abq:
            w(f"*Step, name=LOAD, nlgeom={'YES' if st.nlgeom else 'NO'}, inc=10000\n*Static\n")
        else:
            w(f"*Step{', nlgeom' if st.nlgeom else ''}, inc=10000\n*Static\n")
        first = min(st.max_increment, min(st.load_points))
        w(f"{_fmt(first)}, 1., {_fmt(st.min_increment)}, {_fmt(st.max_increment)}\n")
    else:
        raise ValueError(f"unknown step kind {st.kind!r}")

    cl: Dict[Tuple[int, int], float] = {}

    def add(label, dof, value):
        rot = {ref_labels[e]: rot_labels[e] for e in rot_labels}
        if dof > 3 and label in rot:
            label, dof = rot[label], dof - 3
        if value != 0.0:
            cl[(label, dof)] = cl.get((label, dof), 0.0) + float(value)

    for load in model.loads:
        Fx, Fy, Fz, Mx, My, Mz = map(float, load.resultant)
        if load.where in ("root", "tip"):
            k = 0 if load.where == "root" else mesh.last
            for j, f in enumerate(mesh.plane_section_forces(Fx, My, Mz)):
                add(mesh.label(k, j), 1, f)
            if Fy or Fz or Mx:
                if load.where not in ref_labels:
                    raise ValueError(f"{model.name}: shear or torque at the {load.where} needs a coupled end")
                lab = ref_labels[load.where]
                ry, rz = (model.root if load.where == "root" else model.tip).reference
                add(lab, 2, Fy)
                add(lab, 3, Fz)
                add(lab, 4, Mx - (ry * Fz - rz * Fy))   # move the torque from the origin to the reference
        else:
            if Fx or Mx or My or Mz:
                raise ValueError("interior loads take F_y and F_z only")
            k = mesh.station_index(float(load.where) / model.length)
            y, z = mesh.section_yz.T
            for dof, F, coord in ((2, Fy, z), (3, Fz, y)):
                if F:
                    sel = np.flatnonzero(np.abs(coord - load.at) < 1e-9 * max(1.0, abs(load.at)))
                    if sel.size == 0:
                        raise ValueError(f"no section nodes at {load.at} for an interior load")
                    for j in sel:
                        add(mesh.label(k, j), dof, F / sel.size)
    w("*Cload\n")
    for (lab, dof), v in sorted(cl.items()):
        w(f"{lab}, {dof}, {_fmt(v)}\n")

    # output
    if dialect == "abaqus":
        if st.kind == "static":
            w("*Output, field, time points=LEVELS\n")
        else:
            w("*Output, field\n")
        w("*Node Output, nset=STATIONS\nU, UR\n")
        if ref_labels:
            w("*Node Output, nset=REFS\nU, UR, RF, RM\n")
        w("*Output, history, frequency=0\n")
        if st.kind == "buckle":
            w("*Node File, nset=STATIONS\nU\n")   # for *Imperfection in a follow-on job
    else:
        if st.kind == "static":
            w("*Node Print, nset=STATIONS, time points=LEVELS\nU\n")
            if ref_labels:
                w("*Node Print, nset=REFS, time points=LEVELS\nU\n")
        else:
            w("*Node Print, nset=STATIONS\nU\n")
    w("*End Step\n")
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(out))
    return mesh


def _orientation_name(wall_id: int, angle: float) -> str:
    tag = f"{angle:+g}".replace("+", "P").replace("-", "M").replace(".", "D")
    return f"ORI_W{wall_id}_{tag}"
