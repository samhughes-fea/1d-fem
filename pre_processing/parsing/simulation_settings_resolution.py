# pre_processing/parsing/simulation_settings_resolution.py
"""
Canonical taxonomy §1–§5 resolution for simulation_settings after parsing simulation_settings.txt.

Merges legacy [Modal]/[Dynamic] into eigen/buckling/transient dicts, validates ``enabled`` flags,
and normalizes ``type`` to one of: static, eigen, transient, harmonic, buckling.
"""

from __future__ import annotations

import logging
import os
from typing import Any, Dict

logger = logging.getLogger(__name__)

_SILENCE_LEGACY_WARN_ENV = "FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS"
_LEGACY_MODAL_ERROR_ENV = "FEM_LEGACY_MODAL_ERROR"


def _legacy_warnings_enabled() -> bool:
    return os.environ.get(_SILENCE_LEGACY_WARN_ENV, "").strip().lower() not in (
        "1",
        "true",
        "yes",
        "on",
    )


CANONICAL_SIMULATION_TYPES = frozenset({"static", "eigen", "transient", "harmonic", "buckling"})

# Accepted on the legacy [Simulation] → [Type] line before normalization
LEGACY_TYPE_LINE_ALIASES = frozenset({"modal", "dynamic", "static_nonlinear"})

VALID_TYPE_LINE_INPUTS = CANONICAL_SIMULATION_TYPES | LEGACY_TYPE_LINE_ALIASES


def _legacy_modal_error_strict() -> bool:
    """When true, legacy ``[Modal]`` / ``[Type] modal`` inputs raise instead of warning."""
    return os.environ.get(_LEGACY_MODAL_ERROR_ENV, "").strip().lower() in (
        "1",
        "true",
        "yes",
        "on",
    )


def _enforce_legacy_modal_strict_errors(
    *,
    modal_header_seen: bool,
    type_line_explicit: bool,
    parsed_type_before_finalize: str,
) -> None:
    if not _legacy_modal_error_strict():
        return
    if modal_header_seen:
        raise ValueError(
            f"Legacy [Modal] section is not allowed when {_LEGACY_MODAL_ERROR_ENV}=1; "
            "use [Eigen] for vibration or [Buckling] for linear buckling "
            "(see docs/conventions/SIMULATION_SETTINGS_TAXONOMY.md)."
        )
    pt = str(parsed_type_before_finalize).strip().lower()
    if type_line_explicit and pt == "modal":
        raise ValueError(
            f"Legacy [Type] modal is not allowed when {_LEGACY_MODAL_ERROR_ENV}=1; "
            "use [Type] eigen or buckling with matching [Eigen]/[Buckling] sections."
        )


def _emit_legacy_input_warnings(
    *,
    type_line_explicit: bool,
    parsed_type_before_finalize: str,
    modal_header_seen: bool,
) -> None:
    """Log deprecation notices for [Modal] and legacy [Type] modal/dynamic (opt-out via env)."""
    if not _legacy_warnings_enabled():
        return
    if modal_header_seen:
        logger.warning(
            "[Modal] in simulation_settings.txt is deprecated: use [Eigen] for vibration "
            "(with [Type] eigen) or [Buckling] for linear buckling (with [Type] buckling). "
            "Legacy keys remain merged on read. Silence warnings with %s=1.",
            _SILENCE_LEGACY_WARN_ENV,
        )
    if type_line_explicit and parsed_type_before_finalize in LEGACY_TYPE_LINE_ALIASES:
        if parsed_type_before_finalize == "modal":
            logger.warning(
                "[Type] modal is deprecated; prefer [Type] eigen or buckling (and matching "
                "[Eigen]/[Buckling] sections). Silence warnings with %s=1.",
                _SILENCE_LEGACY_WARN_ENV,
            )
        elif parsed_type_before_finalize == "dynamic":
            logger.warning(
                "[Type] dynamic is deprecated; prefer [Type] transient (and optional [Transient] "
                "section). Silence warnings with %s=1.",
                _SILENCE_LEGACY_WARN_ENV,
            )
        elif parsed_type_before_finalize == "static_nonlinear":
            logger.warning(
                "[Type] static_nonlinear is deprecated; use [Type] static with nonlinear "
                "settings ([Newton], [Nonlinear]). Silence warnings with %s=1.",
                _SILENCE_LEGACY_WARN_ENV,
            )

def effective_transient_config(settings: Dict[str, Any]) -> Dict[str, Any]:
    """Merge [Transient] over [Dynamic]; transient keys win."""
    d = dict(settings.get("dynamic") or {})
    t = dict(settings.get("transient") or {})
    out = {**d, **t}
    return out


def effective_eigen_config(settings: Dict[str, Any]) -> Dict[str, Any]:
    """Vibration eigenproblem settings (legacy ``modal`` + ``eigen`` section)."""
    m = settings.get("modal") or {}
    e = settings.get("eigen") or {}
    return {"num_modes": int(e.get("num_modes", m.get("num_modes", 10)))}


def effective_buckling_config(settings: Dict[str, Any]) -> Dict[str, Any]:
    """Linear buckling settings (legacy ``modal`` + ``buckling`` section)."""
    m = settings.get("modal") or {}
    b = settings.get("buckling") or {}
    return {
        "num_modes": int(b.get("num_modes", m.get("num_modes", 10))),
        "buckling_prestress": str(
            b.get("buckling_prestress", m.get("buckling_prestress", "linear_static"))
        ).lower(),
        "buckling_load_factor": float(b.get("buckling_load_factor", m.get("buckling_load_factor", 1.0))),
        "buckling_nonlinear_prestress_twins": bool(
            b.get(
                "buckling_nonlinear_prestress_twins",
                m.get("buckling_nonlinear_prestress_twins", False),
            )
        ),
    }


def _hydrate_taxonomy_from_legacy_sections(settings: Dict[str, Any]) -> None:
    """Fill eigen/buckling/transient defaults from legacy modal/dynamic dicts."""
    m = settings.get("modal") or {}
    modal_explicit = bool(settings.get("_modal_section_in_input"))
    if not modal_explicit:
        m = {}

    if modal_explicit and str(m.get("analysis", "vibration")).lower() == "buckling":
        bk = settings.setdefault("buckling", {})
        # Parser defaults may already populate taxonomy dicts; modal keys override when set.
        nm = m.get("num_modes")
        if nm is not None:
            bk["num_modes"] = int(nm)
        bp = m.get("buckling_prestress")
        if bp is not None:
            bk["buckling_prestress"] = str(bp).lower()
        blf = m.get("buckling_load_factor")
        if blf is not None:
            bk["buckling_load_factor"] = float(blf)
        if "buckling_nonlinear_prestress_twins" in m:
            bk["buckling_nonlinear_prestress_twins"] = bool(
                m["buckling_nonlinear_prestress_twins"]
            )
    elif modal_explicit:
        eg = settings.setdefault("eigen", {})
        nm = m.get("num_modes")
        if nm is not None:
            eg["num_modes"] = int(nm)

    d = settings.get("dynamic") or {}
    if not settings.get("_dynamic_section_in_input"):
        d = {}
    t = settings.setdefault("transient", {})
    for key in ("time_step", "end_time", "scheme"):
        if key in d and d[key] is not None:
            t[key] = d[key]


def _taxonomy_enabled_sections(settings: Dict[str, Any]) -> list[str]:
    out = []
    for name in ("static", "eigen", "transient", "harmonic", "buckling"):
        sec = settings.get(name)
        if isinstance(sec, dict) and sec.get("enabled", False):
            out.append(name)
    return out


def _normalize_type_token(raw: str, settings: Dict[str, Any]) -> str:
    t = raw.strip().lower()
    if t == "dynamic":
        return "transient"
    if t == "modal":
        if str(settings.get("modal", {}).get("analysis", "vibration")).lower() == "buckling":
            return "buckling"
        return "eigen"
    if t == "static_nonlinear":
        settings["_resolved_static_kind"] = "nonlinear"
        return "static"
    if t == "static":
        settings.setdefault("_resolved_static_kind", "linear")
        return "static"
    if t in CANONICAL_SIMULATION_TYPES:
        if t == "static":
            settings.setdefault("_resolved_static_kind", "linear")
        return t
    raise ValueError(f"Invalid simulation type after normalization: {raw!r}")


def finalize_simulation_settings(settings: Dict[str, Any], *, type_line_explicit: bool) -> None:
    """
    Mutates ``settings`` in place: hydrate taxonomy dicts, validate ``enabled`` sections,
    and set canonical ``settings['type']``.
    """
    parsed_type_before_finalize = str(settings.get("type", "static")).lower()
    modal_header_seen = bool(settings.get("_modal_section_in_input"))

    _hydrate_taxonomy_from_legacy_sections(settings)

    _enforce_legacy_modal_strict_errors(
        modal_header_seen=modal_header_seen,
        type_line_explicit=type_line_explicit,
        parsed_type_before_finalize=parsed_type_before_finalize,
    )

    try:
        enabled_secs = _taxonomy_enabled_sections(settings)
        raw_type = str(settings.get("type", "static")).lower()

        if len(enabled_secs) > 1:
            raise ValueError(
                "Only one taxonomy section may have enabled=true among "
                "[Static], [Eigen], [Transient], [Harmonic], [Buckling]; "
                f"found enabled on: {enabled_secs}"
            )

        if len(enabled_secs) == 1:
            sec = enabled_secs[0]
            derived = {
                "static": "static",
                "eigen": "eigen",
                "transient": "transient",
                "harmonic": "harmonic",
                "buckling": "buckling",
            }[sec]
            if type_line_explicit:
                try:
                    canonical_from_line = _normalize_type_token(raw_type, settings)
                except ValueError:
                    canonical_from_line = raw_type
                if canonical_from_line != derived:
                    raise ValueError(
                        f"[Simulation]/[Type] resolves to {canonical_from_line!r} but taxonomy section "
                        f"[{sec.title()}] has enabled=true (expects type {derived!r})."
                    )
            settings["type"] = derived
            logger.debug(
                "Resolved simulation type from enabled taxonomy section [%s] -> %s", sec, derived
            )
            return

        # Legacy: no taxonomy enabled flags — use [Type] line (or default)
        settings["type"] = _normalize_type_token(raw_type, settings)
        if not type_line_explicit and settings["type"] == "static":
            logger.warning("Simulation type not specified, using default: 'static'")

        _emit_legacy_input_warnings(
            type_line_explicit=type_line_explicit,
            parsed_type_before_finalize=parsed_type_before_finalize,
            modal_header_seen=modal_header_seen,
        )
    finally:
        settings.pop("_modal_section_in_input", None)
        settings.pop("_dynamic_section_in_input", None)
