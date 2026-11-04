"""Interlaminar stresses of the 2D section against the CalculiX prism of the same mesh.

The prism (``post_processing/validation_visualisers/iea15mw_sections/make_interlaminar_reference.py``)
gives the stresses at the ``3 x 3`` in-plane Gauss points of every brick in a slice at
mid-length, which carries the example's section forces. Each interface value is taken from
those points as ``SolidSectionAnalysis.interface_stress`` does: a quadratic in the distance
from the interface through the Gauss lines of the ply's two rows nearest it, taken to the
interface and averaged over the two plies, at each of the three Gauss stations ``xi`` along
it. The 2D section is evaluated at the same points, so both sides are fitted from the same
kind of data.
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional

import numpy as np

from pre_processing.section_integration.solid.material import bond_matrix, material_axes

REFERENCE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "post_processing",
                         "validation_visualisers", "iea15mw_sections", "reference", "prism_20_interlaminar.json")


def load_reference(path: str = REFERENCE) -> Optional[dict]:
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


class _PrismField:
    """Interlaminar stresses of the prism at the in-plane Gauss points of each element."""

    def __init__(self, mesh, ref: dict):
        if ref["n_elements_2d"] != len(mesh.elements):
            raise ValueError("the prism reference was made from another mesh")
        self.mesh = mesh
        self.s = np.asarray(ref["gauss_stress"])                   # (n_el, n_ip, 6), section axes
        self.a = float(ref["ip_scale"])
        self.ip = {(int(i), int(j)): k for k, (i, j) in enumerate(ref["ip_xi_eta"])}
        self.levels = sorted({int(j) for _, j in ref["ip_xi_eta"]})
        self._T = {}

    def interlaminar(self, e: int, xi: float, eta: float) -> np.ndarray:
        al = float(self.mesh.alpha[e])
        T = self._T.get(al)
        if T is None:
            T = self._T[al] = np.linalg.inv(bond_matrix(material_axes(0.0, al)))
        sig = self.s[e, self.ip[(int(round(xi / self.a)), int(round(eta / self.a)))]]
        return (T @ sig)[2:5]                                       # sigma_nn, tau_sn, tau_xn


def _side(sa, value, e: int, eta_edge: float, xi: float, forces, etas):
    rows = [e]
    f = sa._row_neighbours.get((e, -int(eta_edge)))
    if f is not None and sa._same_ply(e, f):
        rows.append(f)
    p0 = sa.stress_at(e, xi, eta_edge, forces=forces)
    d, v = [], []
    for r in rows:
        for eta in etas:
            q = sa.stress_at(r, xi, eta, forces=forces)
            d.append(float(np.hypot(q.y - p0.y, q.z - p0.z)))
            v.append(value(r, xi, eta, q))
    coef = np.polyfit(np.asarray(d), np.asarray(v), min(len(d) - 1, 2))
    return coef[-1], (p0.y, p0.z)


def compare(sa, mesh, forces, ref: dict, te: Optional[tuple] = None, te_radius: float = 0.02) -> Dict:
    """Interface stresses ``[sigma_nn, tau_sn, tau_xn]`` of both models on every ply interface."""
    prism = _PrismField(mesh, ref)

    def two_d(r, xi, eta, q):
        return np.asarray(q.interlaminar)

    def pr(r, xi, eta, q):
        return prism.interlaminar(r, xi, eta)

    stations = [prism.a * j for j in prism.levels]
    rows: List[dict] = []
    for e, f in sa.ply_interfaces():
        same = int(mesh.elements[e][3]) == int(mesh.elements[f][0])
        for xi in stations:
            vals = {}
            for name, fn in (("2d", two_d), ("prism", pr)):
                a, x0 = _side(sa, fn, e, 1.0, xi, forces, stations)
                b, _ = _side(sa, fn, f, -1.0, xi if same else -xi, forces, stations)
                vals[name] = 0.5 * (a + b)
            at_te = bool(te is not None and np.hypot(x0[0] - te[0][0], x0[1] - te[0][1]) < te_radius * te[1])
            rows.append({"below": int(e), "above": int(f), "y": float(x0[0]), "z": float(x0[1]),
                         "2d": vals["2d"].tolist(), "prism": vals["prism"].tolist(), "at_trailing_edge": at_te})
    return summarise(rows)


def summarise(rows: List[dict]) -> Dict:
    away = [r for r in rows if not r["at_trailing_edge"]]
    shear = lambda v: float(np.hypot(v[1], v[2]))                  # noqa: E731
    w2 = max(away, key=lambda r: shear(r["2d"]))
    wp = max(away, key=lambda r: shear(r["prism"]))
    P = np.array([r["prism"] for r in away])
    D = np.array([r["2d"] for r in away])
    scale = np.abs(P).max(axis=0)
    return {
        "rows": rows,
        "worst_2d": {"y": w2["y"], "z": w2["z"], "shear": shear(w2["2d"]), "prism_there": shear(w2["prism"]),
                     "below": w2["below"], "above": w2["above"]},
        "worst_prism": {"y": wp["y"], "z": wp["z"], "shear": shear(wp["prism"]), "2d_there": shear(wp["2d"]),
                        "below": wp["below"], "above": wp["above"]},
        "same_location": bool(w2["below"] == wp["below"] and w2["above"] == wp["above"]),
        # differences relative to the largest prism value of each component
        "rms_difference_over_max": (np.sqrt(np.mean((D - P) ** 2, axis=0)) / scale).tolist(),
        "max_difference_over_max": (np.abs(D - P).max(axis=0) / scale).tolist(),
        "n_points": len(away),
    }


def figure(result: Dict, path: str) -> None:
    """Parity of the three interface stresses, one marker per interface point (no fitting)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    away = [r for r in result["rows"] if not r["at_trailing_edge"]]
    P = np.array([r["prism"] for r in away]) / 1e6
    D = np.array([r["2d"] for r in away]) / 1e6
    names = [r"$\sigma_{nn}$", r"$\tau_{sn}$", r"$\tau_{xn}$"]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))
    for k, ax in enumerate(axes):
        lim = 1.05 * max(np.abs(P[:, k]).max(), np.abs(D[:, k]).max())
        ax.plot([-lim, lim], [-lim, lim], color="#b5b3ad", lw=0.8, zorder=0)
        ax.scatter(P[:, k], D[:, k], s=6, color="#2a78d6", alpha=0.7, linewidths=0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{names[k]}, CalculiX prism (MPa)", fontsize=8)
        ax.set_ylabel(f"{names[k]}, 2D section (MPa)", fontsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=7, colors="#52514e")
    fig.suptitle(f"Interlaminar stresses on {len(away)} interface points, trailing edge excluded", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
