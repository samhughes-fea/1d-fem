# post_processing/sensitivity/layup_sensitivity.py
"""Sensitivities of beam responses to ply angles and thicknesses (``layup_sensitivity.txt``).

The section derivatives are exact (``thin_walled.layup_sensitivity``). The element
derivatives ``dK_e`` and ``dM_e`` are the element's own matrices formed with the derivative
section, which is exact because the Timoshenko and Benscoter elements are linear in the
section's stiffness (6 x 6, or the Benscoter 8 x 8), its mass and its warping inertia. The global responses then follow from the solved job:

- **static:** ``K U = F`` with ``F`` independent of the layup, so for ``J = g^T U``
  (a displacement) ``dJ/dp = -lambda^T dK U`` with ``K lambda = g`` on the free DOFs (the
  adjoint), and the compliance ``F^T U`` has ``dC/dp = -U^T dK U``;
- **eigen:** ``d(omega^2) = phi^T (dK - omega^2 dM) phi / phi^T M phi`` for each mode (simple
  eigenvalues), and ``df = d(omega^2) / (8 pi^2 f)``.
- **linear buckling:** ``(K + lambda K_g(U_0)) phi = 0`` with the prestress ``K U_0 = F``, so
  ``d lambda = -phi^T (dK + lambda dK_g) phi / phi^T K_g phi``. ``K_g`` changes with the
  section at fixed ``U_0`` (its forces ``D e`` and second-order terms ``H``: the element's
  ``K_g`` formed with ``dD`` and with ``dH``) and through ``U_0``: ``K_g`` is linear in
  ``U_0``, so ``phi^T K_g(dU_0) phi = g^T dU_0`` with ``g_i = phi^T K_g(e_i) phi``, and one
  adjoint solve ``K mu = g`` gives ``g^T dU_0 = -mu^T dK U_0``. Timoshenko elements only.

``write_layup_sensitivity`` writes ``secondary_results/sensitivity/layup_sensitivity.csv``
with one row per parameter and response: the parameter value, the response, its derivative
(per degree or per metre) and, for thicknesses, the logarithmic derivative ``(t / J) dJ/dt``.
The products are formed element by element (``lambda_e^T dK_e U_e``), so no global
derivative matrix is assembled.

Public API
----------
LAYUP_SENSITIVITY_COLUMNS, element_derivatives, layup_sensitivities, write_layup_sensitivity.
"""

from __future__ import annotations

import csv
import dataclasses
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from pre_processing.parsing.layup_sensitivity_parser import (
    LAYUP_SENSITIVITY_FILE,
    Response,
    parse_layup_sensitivity_file,
)
from pre_processing.section_integration.thin_walled.layup_sensitivity import (
    LayupParameter,
    element_section_derivative,
    section_derivative,
)

__all__ = ["LAYUP_SENSITIVITY_COLUMNS", "element_derivatives", "layup_sensitivities", "write_layup_sensitivity"]

LAYUP_SENSITIVITY_COLUMNS = ("parameter", "laminate_id", "ply", "kind", "value", "response", "response_value",
                             "derivative", "relative_derivative")

#: elements whose stiffness and mass are linear in the section's 6 x 6 stiffness and mass
SUPPORTED = ("LinearTimoshenkoBeamElement3D", "LinearBenscoterBeamElement3D")
BENSCOTER = ("LinearBenscoterBeamElement3D",)
#: elements whose geometric stiffness is linear in the section stiffness and second-order terms
BUCKLING = ("LinearTimoshenkoBeamElement3D",)


def _check_element(el) -> None:
    name = type(el).__name__
    if name not in SUPPORTED:
        raise NotImplementedError(
            f"layup sensitivities need {', '.join(SUPPORTED)} elements; element {el.element_id} is {name}")
    sa = getattr(el, "_assigned_section_analysis", None)
    if sa is not None and getattr(sa, "tapered", False) and getattr(sa, "taper_parts", None) is None:
        raise NotImplementedError(f"element {el.element_id}: layup sensitivities of a tapered element need its "
                                  "taper split by laminate ([Section_Taper_Laminate] in section_taper.txt)")


def _curved(el) -> bool:
    return bool(getattr(el, "_frame_mesh", False)) and getattr(el, "_frame_arc", None) is not None


class _Swapped:
    """``el`` with ``d_section`` as its assigned section (and no pretwisted copy of the old one)."""

    def __init__(self, el, d_section):
        self.el, self.d = el, d_section

    def __enter__(self):
        el = self.el
        self.saved = (el._assigned_section_analysis, el.logger_operator,
                      el.__dict__.pop("_pretwisted_section_analysis", None))
        el._assigned_section_analysis = self.d
        el.logger_operator = None
        return el

    def __exit__(self, *exc):
        el = self.el
        el.__dict__.pop("_pretwisted_section_analysis", None)
        el._assigned_section_analysis, el.logger_operator, pre = self.saved
        if pre is not None:
            el._pretwisted_section_analysis = pre
        return False


def element_derivatives(el, d_section, mass: bool = False) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """``dK_e`` (and ``dM_e``) of element ``el`` for the derivative section ``d_section``.

    The element forms its own matrices with the derivative section in place of its own; a
    pretwisted element's helical coupling (``with_pretwist``) is linear in the section, so it
    is formed the same way. A curved element is a chain of straight sub-elements condensed by
    their static shapes ``Phi``: ``K = Phi^T K_c Phi`` is stationary in ``Phi``, so
    ``dK = Phi^T dK_c Phi``, while its mass takes the change of the shapes too,
    ``dPhi_i = -K_ii^-1 (dK_c Phi)_i``.
    """
    _check_element(el)
    if _curved(el):
        return _curved_derivatives(el, d_section, mass)
    with _Swapped(el, d_section) as e:
        dK = e.element_stiffness_matrix().K_e
        dM = e.element_mass_matrix().M_e if mass else None
    return dK, dM


def _curved_derivatives(el, d_section, mass: bool):
    from pre_processing.element_library.framed_element import chain

    c = chain(el)
    with _Swapped(c.sub, d_section) as sub:
        dK_sub = np.asarray(sub.element_stiffness_matrix().K_e, dtype=np.float64)
        dM_sub = np.asarray(sub.element_mass_matrix().M_e, dtype=np.float64) if mass else None
    n = (c.m + 1) * c.dpn

    def chain_matrix(A):
        out = np.zeros((n, n))
        for T, g in zip(c.T, c.maps):
            out[np.ix_(g, g)] += T.T @ A @ T
        return out

    dKc = chain_matrix(dK_sub)
    dK = c.Phi.T @ dKc @ c.Phi
    if not mass:
        return 0.5 * (dK + dK.T), None
    Kc = chain_matrix(c.K_sub)
    Mc = chain_matrix(np.asarray(c.sub.element_mass_matrix().M_e, dtype=np.float64))
    dPhi = np.zeros_like(c.Phi)
    if c.i.size:
        dPhi[c.i] = -np.linalg.solve(Kc[np.ix_(c.i, c.i)], (dKc @ c.Phi)[c.i])
    dM = c.Phi.T @ chain_matrix(dM_sub) @ c.Phi + dPhi.T @ Mc @ c.Phi + c.Phi.T @ Mc @ dPhi
    return 0.5 * (dK + dK.T), 0.5 * (dM + dM.T)


def element_geometric_derivative(el, d_section, U_e) -> np.ndarray:
    """The element's geometric stiffness about ``U_e`` formed with the section ``d_section``."""
    with _Swapped(el, d_section) as e:
        return e.linear_geometric_stiffness_matrix(U_e)


def _geometric_gradient(el, phi_e) -> np.ndarray:
    """``g_k = phi_e^T K_g(e_k) phi_e`` over the element DOFs (``K_g`` is linear in ``U_e``)."""
    n = el._n_dof
    g = np.zeros(n)
    for k in range(n):
        e = np.zeros(n)
        e[k] = 1.0
        g[k] = phi_e @ el.linear_geometric_stiffness_matrix(e) @ phi_e
    return g


def layup_sensitivities(elements, thin_walled, params: Sequence[LayupParameter], *, K, U=None, F=None,
                        fixed_dofs=(), M=None, frequencies_hz=None, mode_shapes=None,
                        responses: Sequence[Response] = (), dof_index=None, load_factors=None,
                        prestress=None) -> List[dict]:
    """Rows of ``LAYUP_SENSITIVITY_COLUMNS`` for a solved static (``U``, ``F``), eigen
    (``M``, ``frequencies_hz``, ``mode_shapes``) or linear buckling (``load_factors``,
    ``mode_shapes``, ``prestress`` ``U_0``, ``fixed_dofs``) job.

    ``thin_walled`` is the per-element ``ThinWalledSectionAnalysis`` (``None`` where the
    element has no thin-walled section); ``dof_index(node_id, dof)`` maps a response DOF to
    its global index.
    """
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    K = sp.csr_matrix(K)
    n = K.shape[0]
    buckling = load_factors is not None
    eigen = mode_shapes is not None and not buckling
    for el in elements:
        _check_element(el)
        if buckling and (type(el).__name__ not in BUCKLING or _curved(el)):
            raise NotImplementedError(f"buckling-load sensitivities need straight {', '.join(BUCKLING)} elements; "
                                      f"element {el.element_id} is {type(el).__name__}"
                                      + (", curved or pretwisted" if _curved(el) else ""))
    cache: Dict[Tuple[int, LayupParameter], object] = {}
    rows: List[dict] = []

    dofs = [el.assemble_global_dof_indices() for el in elements]
    if buckling:
        U0 = np.asarray(prestress, dtype=float).ravel()
        Phi = np.asarray(mode_shapes, dtype=float)
        if Phi.ndim == 1:
            Phi = Phi[:, None]
        lams = np.atleast_1d(np.asarray(load_factors, dtype=float))
        fixed = np.unique(np.asarray(fixed_dofs, dtype=int))
        free = np.setdiff1d(np.arange(n), fixed)
        lu = spla.splu(sp.csc_matrix(K[free][:, free]))
        modes = {}
        for r in responses:
            if r.type != "buckling":
                continue
            k = r.mode - 1
            if k >= Phi.shape[1]:
                raise ValueError(f"asked for buckling mode {r.mode}; the job has {Phi.shape[1]}")
            phi = Phi[:, k].copy()
            phi[fixed] = 0.0
            g = np.zeros(n)
            kg = 0.0
            for el, d in zip(elements, dofs):
                g[d] += _geometric_gradient(el, phi[d])
                kg += float(phi[d] @ el.linear_geometric_stiffness_matrix(U0[d]) @ phi[d])
            mu = np.zeros(n)
            mu[free] = lu.solve(g[free])
            modes[r] = (float(lams[k]), phi, mu, kg)
    elif not eigen:
        U = np.asarray(U, dtype=float).ravel()
        F = np.asarray(F, dtype=float).ravel()
        fixed = np.unique(np.asarray(fixed_dofs, dtype=int))
        free = np.setdiff1d(np.arange(n), fixed)
        lu = spla.splu(sp.csc_matrix(K[free][:, free]))
        adj = {}
        for r in responses:
            if r.type == "displacement":
                g = np.zeros(n)
                g[dof_index(r.node_id, r.dof)] = 1.0
                lam = np.zeros(n)
                lam[free] = lu.solve(g[free])
                adj[r] = (float(U[dof_index(r.node_id, r.dof)]), lam)
            elif r.type == "compliance":
                adj[r] = (float(F @ U), U)                # self-adjoint: dC = -U^T dK U
    else:
        Phi = np.asarray(mode_shapes, dtype=float)
        if Phi.ndim == 1:
            Phi = Phi[:, None]
        f = np.atleast_1d(np.asarray(frequencies_hz, dtype=float))

    def form(A_e, x, y):
        return sum(float(x[g] @ A @ y[g]) for g, A in zip(dofs, A_e))

    for p in params:
        dKe, dMe, dGe = [], [], []
        for el, tw in zip(elements, thin_walled):
            if tw is None:
                dKe.append(np.zeros((el._n_dof, el._n_dof)))
                dMe.append(np.zeros((el._n_dof, el._n_dof)))
                dGe.append((None, None, None))
                continue
            bens = type(el).__name__ in BENSCOTER
            key = (id(tw), p, bens)
            if key not in cache:
                cache[key] = section_derivative(tw, p, benscoter=bens, geometric=buckling)
            dsa = element_section_derivative(el.section_analysis, cache[key])
            a, b = element_derivatives(el, dsa, mass=eigen)
            dKe.append(a)
            dMe.append(b)
            if buckling:
                from pre_processing.section_integration.section_analysis import ELEMENT_SIGN

                dH = dataclasses.replace(el.section_analysis, geometric_terms=cache[key].geometric_terms @ ELEMENT_SIGN)
                dGe.append((el, dsa, dH))
        value = _parameter_value(thin_walled, p)
        if buckling:
            for r, (lam, phi, mu, kg) in modes.items():
                num = form(dKe, phi, phi) - lam * form(dKe, mu, U0)
                for (el, dsa, dH), d in zip(dGe, dofs):
                    if el is None:
                        continue
                    u0 = U0[d]
                    dKg = element_geometric_derivative(el, dsa, u0) + element_geometric_derivative(el, dH, u0)
                    num += lam * float(phi[d] @ dKg @ phi[d])
                rows.append(_row(p, value, r, lam, -num / kg))
        elif not eigen:
            for r, (J, lam) in adj.items():
                rows.append(_row(p, value, r, J, -form(dKe, lam, U)))
        else:
            for r in responses:
                if r.type != "frequency":
                    continue
                k = r.mode - 1
                if k >= Phi.shape[1]:
                    raise ValueError(f"asked for mode {r.mode}; the job has {Phi.shape[1]}")
                phi = Phi[:, k]
                w2 = (2 * np.pi * f[k]) ** 2
                dw2 = (form(dKe, phi, phi) - w2 * form(dMe, phi, phi)) / float(phi @ np.asarray(M @ phi).ravel())
                rows.append(_row(p, value, r, float(f[k]), dw2 / (8 * np.pi**2 * f[k])))
    return rows


def _parameter_value(thin_walled, p: LayupParameter) -> float:
    for tw in thin_walled:
        if tw is not None and p.laminate_id in tw.laminates:
            ply = tw.laminates[p.laminate_id].plies[p.ply]
            return float(ply.angle_deg if p.kind == "angle" else ply.thickness)
    raise ValueError(f"no section uses laminate {p.laminate_id}")


def _row(p: LayupParameter, value: float, r: Response, J: float, dJ: float) -> dict:
    rel = value / J * dJ if (p.kind == "thickness" and J != 0.0) else ""
    return {"parameter": p.label, "laminate_id": p.laminate_id, "ply": p.ply + 1, "kind": p.kind, "value": value,
            "response": r.label, "response_value": J, "derivative": dJ, "relative_derivative": rel}


def write_layup_sensitivity(job_dir: str, job_results_dir: str, element_dictionary, elements, runner,
                            solver_type: str) -> Optional[str]:
    """Write ``layup_sensitivity.csv`` for a job with ``layup_sensitivity.txt``; the path, or ``None``."""
    params, responses = parse_layup_sensitivity_file(os.path.join(job_dir, LAYUP_SENSITIVITY_FILE))
    if not params:
        return None
    from pre_processing.element_library.dof_layout import dof_names, mesh_dof_per_node

    tw = element_dictionary.get("thin_walled_analysis_per_element")
    if tw is None:
        raise ValueError("layup sensitivities need thin-walled sections (laminate.txt, section_geometry.txt)")
    dpn = mesh_dof_per_node(element_dictionary)
    names = dof_names(dpn)

    def dof_index(node_id, dof):
        name = "W" if (dof == "Q1" and "W" in names) else dof
        if name not in names:
            raise ValueError(f"unknown DOF {dof!r}; the mesh has {names}")
        return int(node_id) * dpn + names.index(name)

    if solver_type == "static":
        rows = layup_sensitivities(elements, tw, params, K=runner.K_global, U=runner.U_global, F=runner.F_global,
                                   fixed_dofs=runner.fixed_dofs, responses=responses, dof_index=dof_index)
    elif solver_type == "eigen":
        rows = layup_sensitivities(elements, tw, params, K=runner.K_global, M=runner.M_global,
                                   frequencies_hz=runner.vibration_frequencies_hz,
                                   mode_shapes=runner.vibration_mode_shapes, responses=responses)
    elif solver_type == "buckling":
        rows = layup_sensitivities(elements, tw, params, K=runner.buckling_stiffness,
                                   fixed_dofs=runner.buckling_constrained_dofs,
                                   load_factors=runner.buckling_load_factors,
                                   mode_shapes=runner.buckling_mode_shapes,
                                   prestress=runner.buckling_prestress_u, responses=responses)
    else:
        raise ValueError(f"layup sensitivities need a static, eigen or linear buckling job, not {solver_type!r}")
    out_dir = os.path.join(job_results_dir, "secondary_results", "sensitivity")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "layup_sensitivity.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=LAYUP_SENSITIVITY_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    return path
