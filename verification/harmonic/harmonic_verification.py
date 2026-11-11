"""
Harmonic frequency-response verification: FEM vs EB SDOF analytical vs Timoshenko SDOF.

Reads primary_results/harmonic/frequency_response.csv from result directories under
results/. Loads properties from the pinned benchmark job at
jobs/job_benchmark_harmonic_sdof/.

Outputs:
  - harmonic_plots/harmonic_frf_verification.png  (FRF magnitude + phase overlay)
  - harmonic_plots/harmonic_frf_verification.csv
"""
from __future__ import annotations

import csv
import glob
import math
import sys
from pathlib import Path
from typing import Final

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
PROJECT_ROOT: Final[Path] = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

_BENCHMARK_JOB_DIR: Final[Path] = PROJECT_ROOT / "jobs" / "job_benchmark_harmonic_sdof"
_RESULTS_DIR: Final[Path] = PROJECT_ROOT / "results"


def _read_frequency_response(csv_path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            rows.append((
                float(row["frequency_hz"]),
                float(row["response_real"]),
                float(row["response_imag"]),
            ))
    return rows


def run_harmonic_verification() -> None:
    from verification.harmonic.analytical_harmonic import (
        harmonic_sdof_frequency_response,
        load_harmonic_sdof_benchmark,
    )
    from verification.harmonic.timoshenko_analytical_harmonic import (
        load_timoshenko_harmonic_benchmark,
        timoshenko_sdof_frequency_response,
        timoshenko_sdof_natural_frequency_hz,
    )

    out_dir = SCRIPT_DIR / "harmonic_plots"
    out_dir.mkdir(exist_ok=True)

    pattern = str(_RESULTS_DIR / "*" / "primary_results" / "harmonic" / "frequency_response.csv")
    found = sorted(glob.glob(pattern))
    if not found:
        print("No FEM harmonic results found. Run the benchmark job first.")
        return

    eb_benchmark = load_harmonic_sdof_benchmark(_BENCHMARK_JOB_DIR)
    t_benchmark = load_timoshenko_harmonic_benchmark(_BENCHMARK_JOB_DIR)

    eb_rows = harmonic_sdof_frequency_response(eb_benchmark)
    t_rows = timoshenko_sdof_frequency_response(t_benchmark)
    f_n_eb = eb_benchmark.natural_frequency_hz
    f_n_t = timoshenko_sdof_natural_frequency_hz(t_benchmark)

    eb_freqs = np.array([r[0] for r in eb_rows])
    eb_amp = np.array([abs(complex(r[1], r[2])) for r in eb_rows])
    eb_phase = np.degrees([math.atan2(-r[2], r[1]) for r in eb_rows])

    t_freqs = np.array([r[0] for r in t_rows])
    t_amp = np.array([abs(complex(r[1], r[2])) for r in t_rows])
    t_phase = np.degrees([math.atan2(-r[2], r[1]) for r in t_rows])

    csv_out_rows: list[list] = []

    for csv_path_str in found:
        csv_path = Path(csv_path_str)
        result_name = csv_path.parents[2].name
        fem_rows = _read_frequency_response(csv_path)
        if not fem_rows:
            print(f"  {result_name}: empty frequency_response.csv, skipping.")
            continue

        fem_freqs = np.array([r[0] for r in fem_rows])
        fem_amp = np.array([abs(complex(r[1], r[2])) for r in fem_rows])
        fem_phase = np.degrees([math.atan2(-r[2], r[1]) for r in fem_rows])

        # Check if frequency grids are compatible for comparison
        grids_match = (
            len(fem_rows) == len(eb_rows)
            and np.allclose(fem_freqs, eb_freqs, atol=1e-6)
        )

        fig, (ax_amp, ax_phase) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        fig.suptitle(
            f"Harmonic FRF verification — {result_name}\n"
            f"EB SDOF: f_n={f_n_eb:.2f} Hz,  "
            f"Timoshenko SDOF: f_n={f_n_t:.2f} Hz,  "
            f"ζ={eb_benchmark.damping_ratio:.3f}",
            fontsize=11,
        )

        if grids_match:
            ax_amp.semilogy(fem_freqs, fem_amp, "b--", linewidth=1.5, label="FEM", zorder=3)
        ax_amp.semilogy(eb_freqs, eb_amp, "r-", linewidth=1.2, label=f"EB SDOF (f_n={f_n_eb:.2f} Hz)")
        ax_amp.semilogy(t_freqs, t_amp, "g-.", linewidth=1.2, label=f"Timoshenko SDOF (f_n={f_n_t:.2f} Hz)")
        ax_amp.set_ylabel("|H(f)| [m]")
        ax_amp.set_title("FRF magnitude")
        ax_amp.legend()
        ax_amp.grid(True, which="both", ls="--", alpha=0.5)

        if grids_match:
            ax_phase.plot(fem_freqs, fem_phase, "b--", linewidth=1.5, label="FEM", zorder=3)
        ax_phase.plot(eb_freqs, eb_phase, "r-", linewidth=1.2, label=f"EB SDOF (f_n={f_n_eb:.2f} Hz)")
        ax_phase.plot(t_freqs, t_phase, "g-.", linewidth=1.2, label=f"Timoshenko SDOF (f_n={f_n_t:.2f} Hz)")
        ax_phase.set_xlabel("Frequency [Hz]")
        ax_phase.set_ylabel("Phase [deg]")
        ax_phase.set_title("FRF phase")
        ax_phase.legend()
        ax_phase.grid(ls="--", alpha=0.5)

        if not grids_match:
            msg = (
                f"FEM grid ({len(fem_rows)} pts, "
                f"{fem_freqs[0]:.1f}–{fem_freqs[-1]:.1f} Hz) differs from analytical grid "
                f"({len(eb_rows)} pts, {eb_freqs[0]:.1f}–{eb_freqs[-1]:.1f} Hz) — "
                "FEM curve omitted from overlay."
            )
            fig.text(0.5, 0.02, msg, ha="center", fontsize=8, color="orange")

        fig.tight_layout(rect=(0, 0.04, 1, 1))
        out_png = out_dir / f"harmonic_frf_verification_{result_name}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"  Saved: {out_png}")

        # Save analytical rows to CSV (FEM only when grid matches)
        for i in range(len(eb_rows)):
            row = [
                result_name,
                eb_freqs[i],
                fem_amp[i] if grids_match else float("nan"),
                eb_amp[i],
                t_amp[i],
                fem_phase[i] if grids_match else float("nan"),
                eb_phase[i],
                t_phase[i],
            ]
            csv_out_rows.append(row)

        note = "" if grids_match else " [grid mismatch — analytical only]"
        print(
            f"  {result_name}: FEM pts={len(fem_rows)}, analytical pts={len(eb_rows)}"
            f"{note}"
        )

    if csv_out_rows:
        csv_out = out_dir / "harmonic_frf_verification.csv"
        header = [
            "result", "frequency_hz",
            "amplitude_fem_m", "amplitude_eb_m", "amplitude_timoshenko_m",
            "phase_fem_deg", "phase_eb_deg", "phase_timoshenko_deg",
        ]
        with open(csv_out, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(csv_out_rows)
        print(f"  Saved: {csv_out}")

    print("Harmonic verification done.")


if __name__ == "__main__":
    run_harmonic_verification()
