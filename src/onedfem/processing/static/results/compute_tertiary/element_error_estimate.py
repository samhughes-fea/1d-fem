# processing/static/results/compute_tertiary/element_error_estimate.py
"""Element error indicators from Gauss-level resultants. No adaptivity loop.

The estimator is the Zienkiewicz–Zhu energy-norm indicator on the beam
stress-resultant field already stored at Gauss points.

Formula
    Recovered resultants ``S*`` at a Gauss point are the linear interpolant of
    the nodal-projected secondary stresses at the element's two nodes
    ``S*(xi) = N1(xi) S_node0 + N2(xi) S_node1``, ``N1 = (1 - xi)/2``,
    ``N2 = (1 + xi)/2``. The finite-element field ``S_h`` is the Gauss stress
    from ``S_h = D @ B @ U_e``. The element indicator is

        eta_e**2 = sum_g w_g * |J| * (S* - S_h).T @ pinv(D) @ (S* - S_h)

    with ``D`` and ``|J|`` from the formulation cache at that Gauss point.

Source
    Zienkiewicz, O. C. and Zhu, J. Z. (1987). "A simple error estimator and
    adaptive procedure for practical engineering analysis." *International
    Journal for Numerical Methods in Engineering* 24:337–357. Superconvergent
    patch recovery as later stated by the same authors is not run here: the
    recovered field is the suite's existing nodal projection of Gauss
    resultants. This module does not refine the mesh.

Public API
----------
ElementErrorEstimate, element_error_indicators.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
import numpy.typing as npt

__all__ = ["ElementErrorEstimate", "element_error_indicators"]

ZZ_SOURCE = (
    "Zienkiewicz, O. C. and Zhu, J. Z. (1987). International Journal for "
    "Numerical Methods in Engineering 24:337-357. Energy-norm indicator from "
    "Gauss resultants versus the nodal-projected recovered field. No adaptivity."
)


@dataclass
class ElementErrorEstimate:
    """Per-element ZZ energy-norm indicators. ``eta_e`` has one entry per element."""

    eta_e: npt.NDArray[np.float64]
    eta_rel_e: npt.NDArray[np.float64]
    strain_energy_e: npt.NDArray[np.float64]
    eta_global: float
    source: str = ZZ_SOURCE


def _connectivity_row(element_dictionary, element_id: int) -> Sequence[int]:
    conn = np.asarray(element_dictionary["connectivity"])
    ids = element_dictionary.get("ids")
    if ids is None:
        return conn[int(element_id)]
    ids = np.asarray(ids)
    matches = np.flatnonzero(ids.astype(int) == int(element_id))
    if matches.size == 0:
        return conn[int(element_id)]
    return conn[int(matches[0])]


def element_error_indicators(
    gaussian_results,
    formulation_cache,
    nodal_results,
    element_dictionary,
) -> Optional[ElementErrorEstimate]:
    """ZZ energy-norm indicators from existing Gauss and nodal-projected resultants.

    Returns None when Gauss stress, nodal stress, or the formulation cache is
    missing. Does not change the mesh.
    """
    if gaussian_results is None or formulation_cache is None or nodal_results is None:
        return None
    if gaussian_results.stress is None or nodal_results.stress is None:
        return None
    objects = list(formulation_cache.element_objects)
    if not objects:
        return None
    nodal = np.asarray(nodal_results.stress, dtype=np.float64)
    n_e = len(objects)
    eta2 = np.zeros(n_e, dtype=np.float64)
    energy = np.zeros(n_e, dtype=np.float64)
    gauss_stress = gaussian_results.stress
    gauss_energy = gaussian_results.internal_energy_density

    for i, elem_obj in enumerate(objects):
        n0, n1 = (int(v) for v in _connectivity_row(element_dictionary, elem_obj.element_id))
        s0 = np.asarray(nodal[n0], dtype=np.float64).ravel()
        s1 = np.asarray(nodal[n1], dtype=np.float64).ravel()
        stresses = gauss_stress[i]
        densities = None if gauss_energy is None else gauss_energy[i]
        for k, gp in enumerate(elem_obj.gauss_data):
            xi = float(gp.xi)
            n1s, n2s = 0.5 * (1.0 - xi), 0.5 * (1.0 + xi)
            s_h = np.asarray(stresses[k], dtype=np.float64).ravel()
            n_comp = min(s_h.size, s0.size, s1.size)
            s_star = n1s * s0[:n_comp] + n2s * s1[:n_comp]
            delta = s_star - s_h[:n_comp]
            D = np.asarray(gp.D_matrix, dtype=np.float64)
            n_d = min(D.shape[0], D.shape[1], n_comp)
            if n_d == 0:
                continue
            D_use = D[:n_d, :n_d]
            residual = delta[:n_d]
            density = float(residual @ np.linalg.pinv(D_use) @ residual)
            wJ = float(gp.weight) * float(gp.jacobian)
            eta2[i] += wJ * density
            if densities is not None:
                energy[i] += wJ * float(np.asarray(densities[k]).reshape(-1)[0])
        eta2[i] = max(eta2[i], 0.0)

    eta = np.sqrt(eta2)
    eta_rel = np.zeros_like(eta)
    for i, (e, u) in enumerate(zip(eta, energy)):
        if u > 0.0:
            eta_rel[i] = e / np.sqrt(2.0 * u)
        elif e == 0.0:
            eta_rel[i] = 0.0
        else:
            eta_rel[i] = np.nan
    return ElementErrorEstimate(
        eta_e=eta,
        eta_rel_e=eta_rel,
        strain_energy_e=energy,
        eta_global=float(np.sqrt(np.sum(eta2))),
    )
