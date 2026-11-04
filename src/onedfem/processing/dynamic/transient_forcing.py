# processing/dynamic/transient_forcing.py
"""Time-varying nodal loads for §3 transient integration (file or analytic envelope)."""

from __future__ import annotations

import logging
import os
from typing import Callable, Optional

import numpy as np

_logger = logging.getLogger(__name__)


def _resolve_path(job_dir: Optional[str], rel_or_abs: str) -> str:
    p = str(rel_or_abs).strip()
    if os.path.isabs(p):
        return p
    if job_dir:
        return os.path.normpath(os.path.join(job_dir, p))
    return os.path.normpath(p)


def _interp_timeseries(t: float, times: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Linear interpolation of vector ``values`` over ``times`` (broadcast-safe)."""
    t = float(t)
    if times.size == 0:
        return values[0].copy() if values.ndim > 1 else np.array([values[0]], dtype=np.float64)
    if t <= float(times[0]):
        return np.asarray(values[0], dtype=np.float64).ravel().copy()
    if t >= float(times[-1]):
        return np.asarray(values[-1], dtype=np.float64).ravel().copy()
    idx = int(np.searchsorted(times, t, side="right") - 1)
    idx = max(0, min(idx, times.size - 2))
    t0, t1 = float(times[idx]), float(times[idx + 1])
    w = 0.0 if t1 <= t0 else (t - t0) / (t1 - t0)
    v0 = np.asarray(values[idx], dtype=np.float64).ravel()
    v1 = np.asarray(values[idx + 1], dtype=np.float64).ravel()
    return (1.0 - w) * v0 + w * v1


def _load_force_time_series(
    path: str, total_dof: int
) -> tuple[str, np.ndarray, np.ndarray]:
    """
    Load ``(times, data)`` from whitespace-delimited numeric file.

    * Two columns ``time scale``: scalar multiplier per row (applied to ``F_ref``).
    * ``1 + total_dof`` columns: ``time`` plus full global force vector at each time.
    """
    raw = np.loadtxt(path, dtype=np.float64)
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    if raw.shape[1] < 2:
        raise ValueError(f"force_time_series_file must have >= 2 columns, got shape {raw.shape}")
    times = raw[:, 0]
    if raw.shape[1] == 2:
        return "scale", times, raw[:, 1:2]
    if raw.shape[1] == total_dof + 1:
        return "full", times, raw[:, 1:]
    raise ValueError(
        f"force_time_series_file: expected 2 columns (t, scale) or {total_dof + 1} columns "
        f"(t + full F), got {raw.shape[1]}"
    )


def _analytic_envelope(
    t: float,
    *,
    kind: str,
    amp: float,
    freq_hz: float,
    phase_rad: float,
    t_start: Optional[float],
    t_end: Optional[float],
) -> float:
    k = str(kind or "none").strip().lower()
    if k in ("", "none", "off"):
        return 1.0
    if k == "sin":
        return 1.0 + float(amp) * np.sin(2.0 * np.pi * float(freq_hz) * t + float(phase_rad))
    if k == "sin_burst":
        t0 = float(t_start if t_start is not None else 0.0)
        t1 = float(t_end if t_end is not None else t0 + 1.0)
        if t < t0 or t > t1:
            return 1.0
        return 1.0 + float(amp) * np.sin(2.0 * np.pi * float(freq_hz) * t + float(phase_rad))
    raise ValueError(f"Unknown force_analytic kind {kind!r} (use none, sin, sin_burst)")


def build_transient_force_func(
    F_ref: np.ndarray,
    dyn_config: dict,
    *,
    total_dof: int,
    job_dir: Optional[str],
    end_time: float,
) -> Callable[[float], np.ndarray]:
    """
    Build ``F(t)`` from assembled reference load ``F_ref`` and merged transient settings.

    Combines ``load_scale``, optional ``load_ramp``, optional ``force_time_series_file``,
    and optional ``force_analytic`` envelope (multiplicative on the scaled reference pattern).
    """
    F_ref = np.asarray(F_ref, dtype=np.float64).ravel()
    scale = float(dyn_config.get("load_scale", 1.0))
    ramp = bool(dyn_config.get("load_ramp", False))
    denom = max(float(end_time), 1e-15)

    ts_path = dyn_config.get("force_time_series_file")
    mode = "scale"
    times = np.array([0.0, denom], dtype=np.float64)
    series = np.ones((2, 1), dtype=np.float64)
    full_F: Optional[np.ndarray] = None

    if ts_path:
        abs_p = _resolve_path(job_dir, str(ts_path))
        if not os.path.isfile(abs_p):
            raise FileNotFoundError(f"force_time_series_file not found: {abs_p}")
        mode, times, series = _load_force_time_series(abs_p, total_dof)
        _logger.info("Loaded transient force time series from %s mode=%s", abs_p, mode)
        if mode == "full":
            full_F = series
            series = np.ones((times.size, 1), dtype=np.float64)

    fa = dyn_config.get("force_analytic")
    fa_amp = float(dyn_config.get("force_analytic_amplitude", 0.0) or 0.0)
    fa_f = float(dyn_config.get("force_analytic_frequency_hz", 0.0) or 0.0)
    fa_ph = float(dyn_config.get("force_analytic_phase_rad", 0.0) or 0.0)
    fa_t0 = dyn_config.get("force_analytic_t_start")
    fa_t1 = dyn_config.get("force_analytic_t_end")
    fa_t0f = float(fa_t0) if fa_t0 is not None else None
    fa_t1f = float(fa_t1) if fa_t1 is not None else None

    def F_func(t: float) -> np.ndarray:
        tr = float(t)
        base = F_ref * scale
        if ramp:
            base = base * min(1.0, tr / denom)
        if mode == "full" and full_F is not None:
            fvec = _interp_timeseries(tr, times, full_F)
            if ramp:
                fvec = fvec * min(1.0, tr / denom)
            out = fvec * scale
        else:
            s = float(_interp_timeseries(tr, times, series)[0])
            out = base * s
        env = _analytic_envelope(
            tr,
            kind=str(fa or "none"),
            amp=fa_amp,
            freq_hz=fa_f,
            phase_rad=fa_ph,
            t_start=fa_t0f,
            t_end=fa_t1f,
        )
        return np.asarray(out * env, dtype=np.float64).ravel()

    return F_func
