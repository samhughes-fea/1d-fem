# post_processing/sensitivity/layup_sensitivity.py
"""Sensitivities of beam responses to ply angles and thicknesses (``layup_sensitivity.txt``).

The section derivatives are exact (``thin_walled.layup_sensitivity``). The element
derivatives ``dK_e`` and ``dM_e`` are the element's own matrices formed with the derivative
section, which is exact because the Timoshenko family is linear in the section's 6 x 6
stiffness and mass. The global responses then follow from the solved job:

- **static:** ``K U = F`` with ``F`` independent of the layup, so for ``J = g^T U``
  (a displacement) ``dJ/dp = -lambda^T dK U`` with ``K lambda = g`` on the free DOFs (the
  adjoint), and the compliance ``F^T U`` has ``dC/dp = -U^T dK U``;
- **eigen:** ``d(omega^2) = phi^T (dK - omega^2 dM) phi / phi^T M phi`` for each mode (simple
  eigenvalues), and ``df = d(omega^2) / (8 pi^2 f)``.

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
SUPPORTED = ("LinearTimoshenkoBeamElement3D",)


def _check_element(el) -> None:
    name = type(el).__name__
    if name not in SUPPORTED:
        raise NotImplementedError(
            f"layup sensitivities need {', '.join(SUPPORTED)} elements; element {el.element_id} is {name}")
    if getattr(el, "_pretwist_rate", 0.0):
        raise NotImplementedError(f"element {el.element_id}: layup sensitivities of pretwisted elements")
    if getattr(el, "_section_taper", None) is not None:
        raise NotImplementedError(f"element {el.element_id}: layup sensitivities of tapered elements")


def element_derivatives(el, d_section, mass: bool = False) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """``dK_e`` (and ``dM_e``) of element ``el`` for the derivative section ``d_section``.

    The element forms its own matrices with the derivative section in place of its own.
    """
    _check_element(el)
    saved_sa, saved_log = el._assigned_section_analysis, el.logger_operator
    el._assigned_section_analysis = d_section
    el.logger_operator = None
    try:
        dK = el.element_stiffness_matrix().K_e
        dM = el.element_mass_matrix().M_e if mass else None
    finally:
        el._assigned_section_analysis, el.logger_operator = saved_sa, saved_log
    return dK, dM


def layup_sensitivities(elements, thin_walled, params: Sequence[LayupParameter], *, K, U=None, F=None,
                        fixed_dofs=(), M=None, frequencies_hz=None, mode_shapes=None,
                        responses: Sequence[Response] = (), dof_index=None) -> List[dict]:
    """Rows of ``LAYUP_SENSITIVITY_COLUMNS`` for a solved static (``U``, ``F``) or eigen
    (``M``, ``frequencies_hz``, ``mode_shapes``) job.

    ``thin_walled`` is the per-element ``ThinWalledSectionAnalysis`` (``None`` where the
    element has no thin-walled section); ``dof_index(node_id, dof)`` maps a response DOF to
    its global index.
    """
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    K = sp.csr_matrix(K)
    n = K.shape[0]
    eigen = mode_shapes is not None
    for el in elements:
        _check_element(el)
    cache: Dict[Tuple[int, LayupParameter], object] = {}
    rows: List[dict] = []

    if not eigen:
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

    dofs = [el.assemble_global_dof_indices() for el in elements]

    def form(A_e, x, y):
        return sum(float(x[g] @ A @ y[g]) for g, A in zip(dofs, A_e))

    for p in params:
        dKe, dMe = [], []
        for el, tw in zip(elements, thin_walled):
            if tw is None:
                dKe.append(np.zeros((el._n_dof, el._n_dof)))
                dMe.append(np.zeros((el._n_dof, el._n_dof)))
                continue
            key = (id(tw), p)
            if key not in cache:
                cache[key] = section_derivative(tw, p)
            dsa = element_section_derivative(el.section_analysis, cache[key])
            a, b = element_derivatives(el, dsa, mass=eigen)
            dKe.append(a)
            dMe.append(b)
        value = _parameter_value(thin_walled, p)
        if not eigen:
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
    else:
        raise ValueError(f"layup sensitivities need a static or eigen job, not {solver_type!r}")
    out_dir = os.path.join(job_results_dir, "secondary_results", "sensitivity")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "layup_sensitivity.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=LAYUP_SENSITIVITY_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    return path
