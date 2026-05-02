"""
Shared helpers for optional Vlasov warping on beam elements.

**Contract (single source of truth)**
--------------------------------------

1. **Mesh DOFs:** ``mesh_uses_warping_dof(element_dictionary)`` is true if any row turns on
   warping: explicit ``[warping]`` column (preferred; ``ElementParser`` always emits a
   ``warping`` array) **or** legacy type names containing ``"Warping"`` when ``warping`` is absent.
   When true, the static runner allocates **7 DOF per node** (χ warping at local index 6).

2. **Warping stiffness:** ``element_warping_stiffness_on(element_dictionary, idx, type_str)`` is
   true iff that **element row** assembles ``E·Γ`` on ``D[6,6]`` (or the EB 7×7 equivalent).
   Use ``[warping]=1`` on the row or legacy ``*Warping*`` type when no column exists.

3. **Section Γ:** From ``section.txt`` tier **11** only (key ``"Gamma"`` in ``section_dictionary``).
   If the tier has no Γ column, ``section_array`` has no index 9 — treat Γ as **0** at parse time.
   ``effective_warping_gamma(Γ_section, stiffness_on)`` zeroes Γ for ``D`` when stiffness is off.

Use :func:`beam_warping_policy` to bundle mesh/stiffness/Γ for an element row after
``section_array`` is available.

Strict Γ validation (optional): set ``simulation_settings["warping"]["strict_gamma"]`` or environment
variable ``FEM_WARPING_STRICT_GAMMA=1``. When enabled, ``[warping]`` on with ``Γ≤0`` raises
:class:`ValueError` at element construction.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np

logger = logging.getLogger(__name__)

STRICT_GAMMA_ENV_VAR = "FEM_WARPING_STRICT_GAMMA"
BEAM_SECTION_HINTS_ENV_VAR = "FEM_BEAM_SECTION_HINTS"


def parse_warping_cell(token: str) -> bool:
    """Parse a warping cell from ``element.txt`` (0/1, true/false, case-insensitive)."""
    t = str(token).strip().lower()
    if t in ("1", "true", "yes", "on"):
        return True
    if t in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"Invalid [warping] value {token!r}; use 0/1 or true/false.")


def mesh_uses_warping_dof(element_dictionary: Dict[str, Any]) -> bool:
    """
    True if the global mesh allocates a seventh DOF (χ) per node for beam elements.

    When ``element_dictionary`` contains ``"warping"``, this is true if any row is on.
    Otherwise, legacy behaviour: any element *type* string contains ``"Warping"``.
    """
    wf = element_dictionary.get("warping")
    if wf is not None:
        arr = np.asarray(wf).ravel()
        if arr.size == 0:
            return False
        if arr.dtype == object or arr.dtype.kind in ("U", "S", "O"):
            return any(parse_warping_cell(str(x)) for x in arr)
        return bool(np.any(arr))
    types = element_dictionary.get("types")
    if types is None:
        return False
    return any("Warping" in str(t) for t in np.asarray(types).ravel())


def element_warping_stiffness_on(
    element_dictionary: Dict[str, Any],
    element_idx: int,
    element_type_str: str,
) -> bool:
    """
    True if this element row assembles Vlasov warping stiffness (``E·Γ`` on ``D[6,6]``).

    When ``element_dictionary`` contains ``"warping"``, use that row; otherwise use whether
    the type name contains ``"Warping"`` (legacy).
    """
    wf = element_dictionary.get("warping")
    if wf is not None:
        raw = np.asarray(element_dictionary["warping"])[element_idx]
        if raw.dtype == object or raw.dtype.kind in ("U", "S", "O"):
            return parse_warping_cell(str(raw))
        return bool(raw)
    return "Warping" in str(element_type_str)


def effective_warping_gamma(section_gamma: float, stiffness_on: bool) -> float:
    """Γ used in ``D[6,6]``: zero if warping assembly is off for this element."""
    return float(section_gamma) if stiffness_on else 0.0


def warping_strict_gamma_enabled(element_dictionary: Optional[Dict[str, Any]] = None) -> bool:
    """
    True if strict validation should raise when warping stiffness is requested but Γ≤0.

    Order of precedence: ``element_dictionary['_warping_config']['strict_gamma']`` (set by the job
    runner from ``simulation_settings['warping']['strict_gamma']``), then environment variable
    ``FEM_WARPING_STRICT_GAMMA``.
    """
    if element_dictionary:
        cfg = element_dictionary.get("_warping_config")
        if isinstance(cfg, dict) and cfg.get("strict_gamma"):
            return True
    env = os.environ.get(STRICT_GAMMA_ENV_VAR, "").strip().lower()
    return env in ("1", "true", "yes", "on")


def enforce_strict_section_gamma(
    *,
    element_dictionary: Optional[Dict[str, Any]],
    element_id: int,
    stiffness_on: bool,
    section_gamma: float,
) -> None:
    """Raise ``ValueError`` if strict mode is on and warping stiffness is on but Γ≤0."""
    if not warping_strict_gamma_enabled(element_dictionary):
        return
    if stiffness_on and section_gamma <= 0.0:
        raise ValueError(
            f"element {element_id}: warping stiffness requested ([warping]=1 or legacy Warping type) "
            f"but section warping constant Gamma is missing or <= 0 (strict mode; "
            f"set {STRICT_GAMMA_ENV_VAR} or simulation_settings warping.strict_gamma)"
        )


def maybe_warn_timoshenko_default_kappa(section_array: np.ndarray, element_id: int) -> None:
    """
    If ``FEM_BEAM_SECTION_HINTS`` is enabled and ``section_array`` has no explicit κ tier (length 5),
    log that Timoshenko shear stiffness uses the default rectangular κ = 5/6.

    Thin-walled open sections should set tier 8+ with an explicit ``[kappa]`` in ``section.txt``.
    """
    env = os.environ.get(BEAM_SECTION_HINTS_ENV_VAR, "").strip().lower()
    if env not in ("1", "true", "yes", "on"):
        return
    if section_array.size > 5:
        return
    logger.warning(
        "element %s: section tier has no explicit shear correction κ — Timoshenko uses default 5/6 for rectangular "
        "sections; for thin-walled open sections set tier 8+ with [kappa] (see docs/conventions/BEAM_SHEAR_CORRECTION_AND_THINWALL.md)",
        element_id,
    )


def warn_if_degenerate_warping_stiffness(
    *,
    stiffness_on: bool,
    section_gamma: float,
    element_id: int,
) -> None:
    """Emit one warning per element when warping is on but section ``Gamma<=0`` (Option A)."""
    if stiffness_on and section_gamma <= 0.0:
        logger.warning(
            "element %s: warping stiffness requested but Gamma<=0 in section; D[6,6] has no restoring "
            "stiffness (degenerate warping channel)",
            element_id,
        )


def section_gamma_from_section_array(section_array: np.ndarray) -> float:
    """
    Γ (m⁶) from the element's ``section_array`` slice, or 0 if the tier has no ``Gamma`` column.

    Matches ``Element1DBase`` ordering: index **9** only when 11-column ``section.txt`` was parsed
    (optional shear centre columns occupy 7–8).
    """
    if section_array.size >= 10:
        return float(section_array[9])
    return 0.0


@dataclass(frozen=True)
class BeamWarpingPolicy:
    """Resolved warping flags and Γ for one element row (linear and nonlinear beams)."""

    mesh_allocates_chi_dof: bool
    warping_stiffness_on: bool
    gamma_section: float
    gamma_effective: float


def beam_warping_policy(
    element_dictionary: Dict[str, Any],
    element_idx: int,
    element_type_str: str,
    gamma_section: float,
) -> BeamWarpingPolicy:
    """
    Single bundle for ``mesh_uses_warping_dof``, ``element_warping_stiffness_on``, and
    ``effective_warping_gamma`` — use after ``section_array`` / Γ is known.

    Parameters
    ----------
    gamma_section
        Raw Γ from section row (0 if ``section.txt`` has no ``[Gamma]`` column).
    """
    mesh = mesh_uses_warping_dof(element_dictionary)
    stiff = element_warping_stiffness_on(element_dictionary, element_idx, element_type_str)
    geff = effective_warping_gamma(gamma_section, stiff)
    return BeamWarpingPolicy(
        mesh_allocates_chi_dof=mesh,
        warping_stiffness_on=stiff,
        gamma_section=float(gamma_section),
        gamma_effective=float(geff),
    )
