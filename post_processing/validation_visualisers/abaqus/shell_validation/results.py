# post_processing/validation_visualisers/abaqus/shell_validation/results.py
"""Read shell results and reduce each output station to beam displacements.

Both solvers' results become one table: frames (buckling modes or load levels), each
with a value (the buckling factor or the load factor) and the displacements of the
station nodes. Abaqus results come from ``extract_odb.py`` (a CSV), CalculiX results from
the ``.dat`` file of ``*NODE PRINT``.

A station's beam displacements are the rigid in-plane motion that fits its nodes best
(a weighted Procrustes fit): the twist ``phi`` and the motion ``(v, w)`` of the section
origin. The axial displacement ``u`` of the origin is the least-squares plane through the
nodes' ``u_x``. ``residual`` is what the rigid motion leaves, relative to the largest
in-plane displacement, and tells a global mode from a local (plate) one.

Public API
----------
ShellFrames, read_abaqus_csv, read_calculix_dat, section_motion, station_motions.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

from post_processing.validation_visualisers.abaqus.shell_validation.shell_model import ShellMesh


@dataclass
class ShellFrames:
    kind: str                                         # "buckle" or "static"
    values: List[float] = field(default_factory=list)          # eigenvalue or load factor per frame
    displacements: List[Dict[int, np.ndarray]] = field(default_factory=list)   # label -> (3,)


def read_abaqus_csv(path: str) -> ShellFrames:
    """Rows ``kind, frame, value, label, u1, u2, u3`` written by ``extract_odb.py``."""
    frames: Dict[int, Dict[int, np.ndarray]] = {}
    values: Dict[int, float] = {}
    kind = "static"
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            kind = row["kind"]
            f = int(row["frame"])
            values[f] = float(row["value"])
            frames.setdefault(f, {})[int(row["label"])] = np.array([float(row[c]) for c in ("u1", "u2", "u3")])
    order = sorted(frames)
    return ShellFrames(kind, [values[f] for f in order], [frames[f] for f in order])


_HEADER = re.compile(r"displacements \(vx,vy,vz\) for set (\S+) and time\s+(\S+)")


def read_calculix_dat(path: str, kind: str, node_set: str = "STATIONS") -> ShellFrames:
    """``*NODE PRINT`` blocks of ``node_set`` in a CalculiX ``.dat`` file.

    For a buckling step the factors come from the ``B U C K L I N G   F A C T O R`` table
    and the blocks are the modes in order; for a static step each block is one time point.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    out = ShellFrames(kind)
    factors: List[float] = []
    in_mode = kind != "buckle"     # a buckling step also prints its base state, before the modes
    i = 0
    while i < len(lines):
        line = lines[i]
        if "E I G E N V A L U E" in line:
            in_mode = True
        if "B U C K L I N G   F A C T O R" in line:
            i += 1
            while i < len(lines) and not re.match(r"^\s*\d+\s+\S+E[+-]\d+\s*$", lines[i]):
                i += 1
            while i < len(lines) and re.match(r"^\s*\d+\s+\S+E[+-]\d+\s*$", lines[i]):
                factors.append(float(lines[i].split()[1]))
                i += 1
            continue
        m = _HEADER.search(line)
        if m and m.group(1).upper() == node_set.upper() and not in_mode:
            i += 1
            continue
        if m and m.group(1).upper() == node_set.upper():
            t = float(m.group(2))
            disp: Dict[int, np.ndarray] = {}
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            while i < len(lines) and lines[i].strip():
                parts = lines[i].split()
                disp[int(parts[0])] = np.array([float(v) for v in parts[1:4]])
                i += 1
            out.displacements.append(disp)
            out.values.append(t)
            if kind == "buckle":
                in_mode = False
            continue
        i += 1
    if kind == "buckle":
        if len(factors) < len(out.displacements):
            raise ValueError(f"{path}: {len(out.displacements)} mode blocks but {len(factors)} buckling factors")
        out.values = factors[:len(out.displacements)]
    return out


@dataclass
class SectionMotion:
    u: float
    v: float
    w: float
    phi: float
    residual: float

    def as_list(self) -> List[float]:
        return [self.u, self.v, self.w, self.phi]


def section_motion(mesh: ShellMesh, k: int, disp: Dict[int, np.ndarray], linear: bool = False) -> SectionMotion:
    """Beam displacements of station ``k`` from its node displacements (see the module notes).

    ``linear`` fits the infinitesimal motion ``u_y = v - phi z``, ``u_z = w + phi y``
    instead of a finite rotation: use it for buckling modes and linear solutions, whose
    displacements are not rigid motions at finite size.
    """
    js = mesh.section_nodes(k)
    P = mesh.section_yz[js]
    U = np.array([disp[mesh.label(k, j)] for j in js])
    wts = mesh.section_weight[js] / mesh.section_weight[js].sum()
    scale = max(float(np.abs(U[:, 1:3]).max()), 1e-300)
    if linear:
        n = len(js)
        A = np.zeros((2 * n, 3))
        A[0::2, 0] = 1.0
        A[0::2, 2] = -P[:, 1]
        A[1::2, 1] = 1.0
        A[1::2, 2] = P[:, 0]
        sw = np.repeat(np.sqrt(wts), 2)
        v, w, phi = np.linalg.lstsq(A * sw[:, None], U[:, 1:3].ravel() * sw, rcond=None)[0]
        residual = float(np.abs(A @ np.array([v, w, phi]) - U[:, 1:3].ravel()).max() / scale)
        origin = np.array([v, w])
        phi = float(phi)
    else:
        Q = P + U[:, 1:3]
        cP, cQ = wts @ P, wts @ Q
        H = (P - cP).T @ (wts[:, None] * (Q - cQ))
        phi = float(np.arctan2(H[0, 1] - H[1, 0], H[0, 0] + H[1, 1]))
        c, s = np.cos(phi), np.sin(phi)
        R = np.array([[c, -s], [s, c]])
        fit = (P - cP) @ R.T + cQ
        origin = R @ (-cP) + cQ
        residual = float(np.abs(Q - fit).max() / scale)
    A = np.column_stack([np.ones(len(js)), P])
    coef = np.linalg.lstsq(A * np.sqrt(wts)[:, None], U[:, 0] * np.sqrt(wts), rcond=None)[0]
    return SectionMotion(float(coef[0]), float(origin[0]), float(origin[1]), phi, residual)


def station_motions(mesh: ShellMesh, frames: ShellFrames, stations, linear: bool = False
                    ) -> Dict[str, List[SectionMotion]]:
    """``{x/L: [motion per frame]}`` for the model's output stations."""
    out: Dict[str, List[SectionMotion]] = {}
    for s in stations:
        k = mesh.station_index(s)
        out[f"{mesh.x[k] / mesh.x[-1]:.6g}"] = [section_motion(mesh, k, d, linear) for d in frames.displacements]
    return out
