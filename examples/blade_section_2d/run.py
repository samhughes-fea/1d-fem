"""2D section stresses of an IEA 15 MW blade section: a cap-web junction and a sandwich panel.

Run from the repository root (needs implicit-geometry-kernel, ``1d-fem[geometry]``)::

    python examples/blade_section_2d/run.py

The 70 % span station of the IEA Wind 15 MW blade (``examples/iea15mw``) is cut by
implicit-geometry-kernel and meshed with layered ``Q8`` elements, two rows per ply. It
carries the section forces of the outboard 20 m of the blade under a flapwise tip load
and a torque. The example shows:

1. **Axial stress** over the section, against the thin-walled model's ply stresses at the
   middle of the suction-side spar cap.
2. **Interlaminar stresses on the ply interfaces** (``interface_stresses``), largest where
   the shear web meets the spar cap, where a thin-walled model has none to give.
3. **Through-thickness profiles in a sandwich panel** (gelcoat, glass face, foam core,
   glass face): axial stress and in-plane shear, against the thin-walled model's ply
   stresses at the same point.

Writes figures and ``summary.json`` next to this file.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (REPO, os.path.join(REPO, "examples", "iea15mw")):
    if p not in sys.path:
        sys.path.insert(0, p)
import stations as S  # noqa: E402

INDEX = 20                     # 70 % span
ARM = 20.0                     # m from the tip load to the section
P_FLAP, TORQUE = 200e3, 100e3  # N and N m
PANEL_X = 0.65                 # chord fraction of the suction-side panel probe
COL = {"face": "#2a78d6", "core": "#eda100", "tw": "#0b0b0b"}


def section_forces() -> np.ndarray:
    """``FORCES`` order (N, M_y, M_z, V_y, V_z, T): flapwise shear V_z with its moment
    M_y = -ARM V_z, and the torque."""
    return np.array([0.0, -ARM * P_FLAP, 0.0, 0.0, P_FLAP, TORQUE])


def build(rows_per_ply: int = 2, element_size=None):
    from pre_processing.section_integration.solid.from_geometry import section_mesh_from_cut
    from pre_processing.section_integration.solid.section_analysis import SolidSectionAnalysis
    from pre_processing.section_integration.thin_walled.from_geometry import thin_walled_from_cut

    cut = S.geometry(INDEX).section_at(0.5)
    lams = S.laminates(INDEX)
    mesh = section_mesh_from_cut(cut, lams, rows_per_ply=rows_per_ply, element_size=element_size)
    gs = thin_walled_from_cut(cut, lams)
    return cut, mesh, SolidSectionAnalysis(mesh), gs, gs.analysis()


def element_centres(mesh) -> np.ndarray:
    return np.array([mesh.nodes[c[:4]].mean(axis=0) for c in mesh.elements])


def panel_column(mesh, sa, chord: float, pitch: float):
    """The stack of elements through the suction-side skin at ``PANEL_X`` of the chord, inside first."""
    y0 = (PANEL_X - pitch) * chord
    start = None
    for e, lab in enumerate(mesh.labels):
        X = mesh.nodes[mesh.elements[e][:4]]
        if lab.startswith("skin_suction") and X[:, 0].min() <= y0 < X[:, 0].max():
            start = e
            break
    nb = sa._row_neighbours
    stack = [start]
    for side in (+1, -1):
        e = start
        while (e, side) in nb:
            e = nb[(e, side)]
            stack.append(e)
    cen = element_centres(mesh)
    return sorted(stack, key=lambda e: cen[e, 1])


def panel_profile(sa, mesh, column, forces):
    """Wall-axis stresses through the panel: rows (n, sigma_xx, tau_xs, tau_xn, material id),
    ``n`` the height from the inner face."""
    g = np.linspace(-1, 1, 5)
    rows, base = [], 0.0
    for e in column:
        c = mesh.elements[e]
        lo, hi = mesh.nodes[c[:2]].mean(axis=0), mesh.nodes[c[2:4]].mean(axis=0)
        h = float(np.linalg.norm(hi - lo))
        flip = hi[1] < lo[1]                       # the element's eta runs outside-in
        for eta in g:
            p = sa.stress_at(e, 0.0, eta, forces=forces)
            f = 0.5 * (1.0 + (-eta if flip else eta))
            rows.append((base + f * h, p.stress_wall[0], p.stress_wall[5], p.stress_wall[4], int(mesh.material_ids[e])))
        base += h
    return rows


def thin_walled_panel(tw, gs, y_target: float, forces):
    """Thin-walled ply stresses at the suction-side skin point nearest ``y_target``."""
    from pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

    sec = tw.section
    best = None
    for w in sec.walls:
        if not gs.wall_names[w.wall_id].startswith("skin_suction"):
            continue
        a, b = (np.asarray(p, float) for p in sec.endpoints(w))
        if b[1] < 0 and a[1] < 0:
            continue
        if min(a[0], b[0]) <= y_target <= max(a[0], b[0]):
            s = abs(y_target - a[0]) / abs(b[0] - a[0]) * np.linalg.norm(b - a)
            best = (w.wall_id, s)
            break
    rec = SaintVenantRecovery(tw).at(best[0], best[1], forces, points_per_ply=5)
    lam = tw.laminates[sec.wall(best[0]).laminate_id]
    h = sum(p.thickness for p in lam.plies)
    return rec, [(pp.n + 0.5 * h, pp.stress_wall[0], pp.stress_wall[2], pp.tau_xn) for pp in rec.plies]


def cap_check(sa, mesh, tw, gs, forces):
    """Mean sigma_xx over the suction-side cap carbon: 2D (area-weighted over the cap
    elements) against thin-walled (the carbon plies' stresses, length-weighted over the cap walls)."""
    from pre_processing.section_integration.solid.elements import shape
    from pre_processing.section_integration.thin_walled.stress_recovery import SaintVenantRecovery

    g, wg = np.polynomial.legendre.leggauss(2)
    num = den = 0.0
    for e, lab in enumerate(mesh.labels):
        if not lab.startswith("spar_cap_suction"):
            continue
        X = mesh.nodes[mesh.elements[e]]
        for a, wa in zip(g, wg):
            for b, wb in zip(g, wg):
                da = wa * wb * np.linalg.det(shape(mesh.kinds[e], a, b)[1] @ X)
                num += da * sa.stress_at(e, a, b, forces=forces).stress[0]
                den += da
    s2d = num / den
    sec = tw.section
    rec = SaintVenantRecovery(tw)
    num = den = 0.0
    for w in sec.walls:
        if not gs.wall_names[w.wall_id].startswith("spar_cap_suction"):
            continue
        a, b = (np.asarray(p, float) for p in sec.endpoints(w))
        L = float(np.linalg.norm(b - a))
        r = rec.at(w.wall_id, 0.5 * L, forces, points_per_ply=3)
        plies = tw.laminates[w.laminate_id].plies
        c = [pp.stress_wall[0] for pp in r.plies if plies[pp.ply].material.E1 > 1e11]
        num += L * float(np.mean(c))
        den += L
    stw = num / den
    return {"solid_2d": float(s2d), "thin_walled": float(stw), "difference": float(s2d / stw - 1.0)}


def junction_interfaces(sa, mesh, forces):
    """Interface stresses on all ply interfaces, and the largest interlaminar shear among them."""
    res = sa.interface_stresses(forces=forces, xi=(0.0,))
    lab = mesh.labels
    out = []
    for e, f, y, z, s in res:
        out.append({"below": int(e), "above": int(f), "y": y, "z": z, "sigma_nn": float(s[0]),
                    "tau_sn": float(s[1]), "tau_xn": float(s[2]), "labels": [lab[e], lab[f]]})
    worst = max(out, key=lambda r: np.hypot(r["tau_xn"], r["tau_sn"])) if out else None
    return out, worst


# ---------------------------------------------------------------------------
def figures(ctx, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    mesh, sa, forces = ctx["mesh"], ctx["sa"], ctx["forces"]
    polys = [mesh.nodes[c[:4]] for c in mesh.elements]
    sxx = np.array([sa.stress_at(e, 0.0, 0.0, forces=forces).stress[0] for e in range(len(mesh.elements))])

    # 1. axial stress
    fig, ax = plt.subplots(figsize=(10, 3.6))
    lim = np.abs(sxx).max() / 1e6
    pc = PolyCollection(polys, array=sxx / 1e6, cmap="RdBu_r", edgecolors="none")
    pc.set_clim(-lim, lim)
    ax.add_collection(pc)
    ax.autoscale_view()
    ax.set_aspect("equal")
    cb = fig.colorbar(pc, ax=ax, shrink=0.8)
    cb.set_label("axial stress σxx, MPa", fontsize=8)
    ax.set_title(f"IEA 15 MW, r = {S.station(INDEX)['r']:.1f} m: σxx under a {P_FLAP / 1e3:.0f} kN flap load "
                 f"{ARM:.0f} m outboard and a {TORQUE / 1e3:.0f} kN m torque", fontsize=9, loc="left")
    _style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "axial_stress.png"), dpi=150)
    plt.close(fig)

    # 2. interface shear, whole section and the zoomed junction
    rows, worst = ctx["interfaces"], ctx["worst"]
    Y = np.array([r["y"] for r in rows])
    Z = np.array([r["z"] for r in rows])
    T = np.array([np.hypot(r["tau_xn"], r["tau_sn"]) for r in rows]) / 1e6
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={"width_ratios": [2.2, 1]})
    for ax, zoom in zip(axes, (False, True)):
        ax.add_collection(PolyCollection(polys, facecolors="#e9e8e3", edgecolors="#fcfcfb", linewidths=0.2))
        sc = ax.scatter(Y, Z, c=T, s=6 if not zoom else 22, cmap="viridis", vmin=0, vmax=T.max(), zorder=3)
        ax.set_aspect("equal")
        if zoom and worst is not None:
            h = 0.12
            ax.set_xlim(worst["y"] - h, worst["y"] + h)
            ax.set_ylim(worst["z"] - h, worst["z"] + h)
            ax.set_title("cap-web junction", fontsize=9, loc="left")
        else:
            ax.autoscale_view()
            ax.set_title("interlaminar shear |τ| on every ply interface", fontsize=9, loc="left")
        _style(ax)
    cb = fig.colorbar(sc, ax=axes, shrink=0.85)
    cb.set_label("√(τxn² + τsn²), MPa", fontsize=8)
    fig.savefig(os.path.join(out_dir, "interface_shear.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 3. sandwich panel through-thickness
    prof, twp = ctx["panel"], ctx["panel_tw"]
    names = {mid: n for mid, n in _material_names(mesh).items()}
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), sharey=True)
    for ax, k, lab in zip(axes, (1, 2), ("axial stress σxx, MPa", "in-plane shear τxs, MPa")):
        n = np.array([r[0] for r in prof]) * 1e3
        v = np.array([r[k] for r in prof]) / 1e6
        core = np.array(["foam" in names.get(r[4], "") for r in prof])
        ax.plot(v[~core], n[~core], "o", ms=3, color=COL["face"], label="2D, faces")
        ax.plot(v[core], n[core], "o", ms=3, color=COL["core"], label="2D, core")
        ax.plot([r[k] / 1e6 for r in twp], [r[0] * 1e3 for r in twp], "-", color=COL["tw"], lw=1.2,
                label="thin-walled")
        ax.set_xlabel(lab, fontsize=8)
        ax.axvline(0, color="#c3c2b7", lw=0.6)
        _style(ax)
    axes[0].set_ylabel("height through the panel from the inner face, mm", fontsize=8)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"Suction-side sandwich panel at {100 * PANEL_X:.0f} % chord", fontsize=9, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "sandwich_panel.png"), dpi=150)
    plt.close(fig)


def _material_names(mesh):
    out = {}
    for mid, mat in mesh.materials.items():
        for name, p in S.DATA["materials"].items():
            if abs(mat.E1 - p["E1"]) < 1e-6 * p["E1"] and abs(mat.E2 - p["E2"]) < 1e-6 * p["E2"]:
                out[mid] = name
    return out


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=7, colors="#52514e")


def main(plot: bool = True, rows_per_ply: int = 2, element_size=None, out_dir: str = HERE) -> dict:
    cut, mesh, sa, gs, tw = build(rows_per_ply, element_size)
    m = S.station_model(INDEX)
    forces = section_forces()
    interfaces, worst = junction_interfaces(sa, mesh, forces)
    column = panel_column(mesh, sa, m.chord, m.pitch_axis)
    panel = panel_profile(sa, mesh, column, forces)
    y_panel = float(np.mean([r for r in element_centres(mesh)[column, 0]]))
    _, panel_tw = thin_walled_panel(tw, gs, y_panel, forces)
    cap = cap_check(sa, mesh, tw, gs, forces)
    names = _material_names(mesh)
    core_txn = [abs(r[3]) for r in panel if "foam" in names.get(r[4], "")]
    tw_core_txn = [abs(r[3]) for r in panel_tw]
    summary = {
        "station": INDEX, "r": S.station(INDEX)["r"], "section_forces": forces.tolist(),
        "n_elements": len(mesh.elements), "rows_per_ply": rows_per_ply,
        "cap_sigma_xx": cap,
        "max_interface_shear": None if worst is None else worst,
        "n_interface_points": len(interfaces),
        "panel": {"y": y_panel, "n_elements": len(column),
                  "core_tau_xn_max_2d": float(max(core_txn)) if core_txn else None,
                  "tau_xn_max_thin_walled": float(max(tw_core_txn)) if tw_core_txn else None},
    }
    if plot:
        figures({"mesh": mesh, "sa": sa, "forces": forces, "interfaces": interfaces, "worst": worst,
                 "panel": panel, "panel_tw": panel_tw}, out_dir)
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    return summary


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
