"""
Harmonic FRF mesh convergence: peak amplitude vs n_elements.

Reads primary_results/harmonic/frequency_response.csv from result directories under
results/ that correspond to job_benchmark_harmonic_sdof_n* jobs.

The primary quantity of interest (QoI) is the peak amplitude (max |H(ω)|) over the
frequency grid. Convergence is measured relative to the Richardson-extrapolated value
from the three finest available meshes (n=32, 64, 128 when present), or, if fewer
than three meshes are available, the finest mesh is treated as the reference.

Outputs:
  - harmonic_plots/harmonic_convergence.png  (peak amplitude vs n, and error vs n)
  - harmonic_plots/harmonic_convergence.csv
"""
from __future__ import annotations

import csv
import glob
import math
import re
import sys
from pathlib import Path
from typing import Final

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"
_N_RE: Final[re.Pattern] = re.compile(r"benchmark_harmonic_sdof_n(\d+)")

# Three-grid GCI settings (refinement ratio R = 2 for n = 32, 64, 128)
_N_FINE, _N_MED, _N_COARSE = 128, 64, 32
_R = 2.0
_F_S = 1.25


def _read_peak_amplitude(csv_path: Path) -> float | None:
    try:
        max_amp = 0.0
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                real = float(row["response_real"])
                imag = float(row["response_imag"])
                amp = math.sqrt(real * real + imag * imag)
                if amp > max_amp:
                    max_amp = amp
        return max_amp if max_amp > 0.0 else None
    except Exception:
        return None


def _richardson_extrap(phi_fine: float, phi_med: float, phi_coarse: float) -> float:
    """Richardson extrapolate from three meshes (R=2). Returns extrapolated value."""
    e21 = phi_med - phi_fine
    e32 = phi_coarse - phi_med
    if abs(e21) < 1e-30 or abs(e32) < 1e-30:
        return phi_fine
    ratio = e32 / e21
    if ratio <= 0 or not math.isfinite(ratio):
        return phi_fine
    p = math.log(abs(ratio)) / math.log(_R)
    denom = (_R ** p) - 1.0
    if denom <= 0 or not math.isfinite(denom):
        return phi_fine
    return phi_fine + (phi_fine - phi_med) / denom


def run_harmonic_convergence() -> None:
    out_dir = SCRIPT_DIR / "harmonic_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(
        _RESULTS_DIR
        / "*benchmark_harmonic_sdof_n*"
        / "primary_results"
        / "harmonic"
        / "frequency_response.csv"
    )
    found = sorted(glob.glob(pattern))
    if not found:
        print("No FEM harmonic convergence results found — run benchmark n-variant jobs first.")
        return

    n_to_amp: dict[int, float] = {}
    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        m = _N_RE.search(csv_path.parents[2].name)
        if not m:
            continue
        n = int(m.group(1))
        amp = _read_peak_amplitude(csv_path)
        if amp is not None:
            n_to_amp[n] = amp

    if not n_to_amp:
        print("No valid FRF amplitude values found in harmonic convergence results.")
        return

    ns = sorted(n_to_amp.keys())
    amps = [n_to_amp[n] for n in ns]

    # Compute reference for relative error:
    # Prefer Richardson extrapolate from n=32, 64, 128 if all three are available.
    # Otherwise use the finest available mesh amplitude.
    amp_ref: float | None = None
    ref_label = ""
    if _N_FINE in n_to_amp and _N_MED in n_to_amp and _N_COARSE in n_to_amp:
        amp_ref = _richardson_extrap(n_to_amp[_N_FINE], n_to_amp[_N_MED], n_to_amp[_N_COARSE])
        ref_label = f"Richardson extrapolate (n={_N_COARSE}/{_N_MED}/{_N_FINE})"
    else:
        n_finest = max(n_to_amp.keys())
        amp_ref = n_to_amp[n_finest]
        ref_label = f"finest available mesh (n={n_finest})"

    errs = [100.0 * abs(a - amp_ref) / amp_ref for a in amps] if amp_ref else [0.0] * len(amps)

    fig, (ax_amp, ax_err) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Harmonic FRF peak-amplitude mesh convergence", fontsize=12, fontweight="bold")

    ax_amp.semilogx(ns, [a * 1e6 for a in amps], "^-", color="#8064A2", linewidth=1.5, markersize=6)
    ax_amp.set_xlabel("Number of elements $n$")
    ax_amp.set_ylabel("Peak amplitude [μm/N]")
    ax_amp.set_title("Peak |H(ω)| vs n")
    ax_amp.grid(True, which="both", ls="--", alpha=0.5)

    # Only plot error when there are ≥2 points with meaningful relative error
    meaningful = [(n, e) for n, e in zip(ns, errs) if e > 1e-12]
    if len(meaningful) >= 2:
        ns_m, errs_m = zip(*meaningful)
        ax_err.loglog(ns_m, errs_m, "^-", color="#8064A2", linewidth=1.5, markersize=6,
                      label=f"vs {ref_label}")
        ax_err.set_xlabel("Number of elements $n$")
        ax_err.set_ylabel(r"$|A^{\rm FEM} - A^{\rm ref}| / A^{\rm ref}$ [%]")
        ax_err.set_title("Relative error vs reference")
        ax_err.grid(True, which="both", ls="--", alpha=0.5)
        ax_err.legend(fontsize="small")
    else:
        ax_err.text(0.5, 0.5, "Insufficient data\nfor error plot", ha="center", va="center",
                    transform=ax_err.transAxes, fontsize=11)
        ax_err.set_title("Relative error vs reference")

    fig.tight_layout()
    out_png = out_dir / "harmonic_convergence.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_png}")

    csv_out = out_dir / "harmonic_convergence.csv"
    with open(csv_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["n_elements", "peak_amplitude_m_per_n", "reference_m_per_n", "abs_rel_error_pct"])
        for n, amp in zip(ns, amps):
            err = 100.0 * abs(amp - amp_ref) / amp_ref if amp_ref else float("nan")
            writer.writerow([n, f"{amp:.10e}", f"{amp_ref:.10e}" if amp_ref else "nan", f"{err:.6f}"])
    print(f"Saved: {csv_out}")

    print(f"\n--- Harmonic peak-amplitude convergence  (ref: {ref_label}) ---")
    for n, amp in zip(ns, amps):
        err = 100.0 * abs(amp - amp_ref) / amp_ref if amp_ref else float("nan")
        print(f"  n={n:4d}:  peak_amp = {amp:.4e} m/N   |error| = {err:.4f}%")


if __name__ == "__main__":
    run_harmonic_convergence()
