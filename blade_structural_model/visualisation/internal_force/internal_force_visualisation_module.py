# blade_structural_model\internal_force\internal_force_visualisation_module.py

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from labellines import labelLines

class InternalForceVisualisationModule:
    def __init__(self, tsr_names: list, rR_locii: list[float], output_dir: Path):
        self.fig, self.axs = plt.subplots(3, 3, figsize=(13, 10), sharex="col")
        self.tsr_names = tsr_names
        self.rR_locii = rR_locii
        self.output_dir = output_dir

        self.COLORS = {
            "load":    "#4F81BD",
            "shear":   "#9BBB59",
            "bending": "#C0504D",
            "torsion": "#8064A2",
            "Mz": "#C0504D",
            "My": "#C0504D",
            "T":  "#8064A2",
        }

        self.worst_curves = {
            "Mz": {"tsr": None, "rR": None, "values": None, "max_val": -np.inf},
            "My": {"tsr": None, "rR": None, "values": None, "max_val": -np.inf},
            "T":  {"tsr": None, "rR": None, "values": None, "max_val": -np.inf},
        }

        self.plot_single()

    def plot_single(self, tsr: str, rR, loads: dict, resultants: dict):
        ax = self.axs
        ax[0, 0].plot(rR, loads["F_y"], label=tsr, color=self.COLORS["load"])
        ax[1, 0].plot(rR, resultants["V_y"], label=tsr, color=self.COLORS["shear"])
        ax[2, 0].plot(rR, resultants["M_z"], label=tsr, color=self.COLORS["bending"])

        ax[0, 1].plot(rR, loads["F_z"], label=tsr, color=self.COLORS["load"])
        ax[1, 1].plot(rR, resultants["V_z"], label=tsr, color=self.COLORS["shear"])
        ax[2, 1].plot(rR, resultants["M_y"], label=tsr, color=self.COLORS["bending"])

        ax[0, 2].plot(rR, loads["M_x"], label=tsr, color=self.COLORS["load"])
        ax[1, 2].plot(rR, resultants["T"], label=tsr, color=self.COLORS["torsion"])

        # Update worst-case tracking
        self._update_worst_case(tsr, rR, resultants)

    def _update_worst_case(self, tsr, rR, resultants):
        for key_csv, key_plot in zip(["M_z", "M_y", "T"], ["Mz", "My", "T"]):
            series = resultants[key_csv]
            max_val = np.max(np.abs(series))
            if max_val > self.worst_curves[key_plot]["max_val"]:
                self.worst_curves[key_plot] = {
                    "tsr": tsr,
                    "rR": rR,
                    "values": series,
                    "max_val": max_val
                }

    def _fill_nested_segment(self, ax, x, y_stack, curve_order, color_map, alpha=0.25):
        for i in range(len(y_stack) - 1):
            lower = np.array(y_stack[i])
            upper = np.array(y_stack[i + 1])
            key = curve_order[i]
            mask_up = upper >= lower
            mask_dn = upper < lower
            ax.fill_between(x, lower, upper, where=mask_up, color=color_map[key], alpha=alpha, linewidth=0, zorder=2)
            ax.fill_between(x, lower, upper, where=mask_dn, color=color_map[key], alpha=alpha, linewidth=0, zorder=2)

    def _plot_worst_case_envelope(self):
        ax = self.axs[2, 2]
        curves = {key: self.worst_curves[key]["values"] for key in ["Mz", "My", "T"]}
        rR = self.worst_curves["Mz"]["rR"]

        # Nested segment fill
        order = ["T", "My", "Mz"]
        for i in range(len(rR) - 1):
            xseg = [rR[i], rR[i + 1]]
            y_stack = [
                [0, 0],
                [curves["T"][i], curves["T"][i + 1]],
                [curves["My"][i], curves["My"][i + 1]],
                [curves["Mz"][i], curves["Mz"][i + 1]],
            ]
            self._fill_nested_segment(ax, xseg, y_stack, order, self.COLORS, alpha=0.25)

        # Overlay curves and markers
        for key in order:
            data = self.worst_curves[key]
            ax.plot(data["rR"], data["values"], label=f"${key}(x)$ – {data['tsr']}", color=self.COLORS[key], zorder=5)

        # Mark extrema
        ax.plot(rR[np.argmax(curves["My"])], np.max(curves["My"]), marker='o', color=self.COLORS["My"], label="_nolegend_", zorder=6)
        ax.plot(rR[np.argmax(curves["T"])], np.max(curves["T"]), marker='o', color=self.COLORS["T"], label="_nolegend_", zorder=6)
        ax.plot(rR[np.argmin(curves["Mz"])], np.min(curves["Mz"]), marker='o', color=self.COLORS["Mz"], label="_nolegend_", zorder=6)

        # Label and style
        ax.set_title("Worst-Case Internal Action Envelope")
        ax.set_xlim(-0.02, 1.02)
        ax.grid(True)
        labelLines(ax.get_lines(), zorder=6)
        for loc in self.rR_locii:
            ax.axvline(loc, linestyle="--", color="black", linewidth=1)

    def finalize(self, filename="internal_force_equilibrium.png"):
        ylabels = [
            ("$f_y$ [N/m]", "$f_z$ [N/m]", "$m_x$ [Nm/m]"),
            ("$V_y(x)$ [N]", "$V_z(x)$ [N]", "$T(x)$ [Nm]"),
            ("$M_z(x)$ [Nm]", "$M_y(x)$ [Nm]", ""),
        ]

        # Label and style axes
        for r in range(3):
            for c in range(3):
                ax = self.axs[r][c]
                label = ylabels[r][c]
                if label == "":
                    ax.axis("off")
                    continue
                ax.set_ylabel(label)
                ax.grid(True)
                ax.set_xlim(-0.02, 1.02)
                if (r < 2) and not (r == 1 and c == 2):
                    ax.tick_params(labelbottom=False)

        self.axs[2, 0].set_xlabel("r / R")
        self.axs[2, 1].set_xlabel("r / R")
        self.axs[1, 2].set_xlabel("r / R")
        self.axs[1, 2].tick_params(labelbottom=True)

        # Add vertical reference lines
        for ax in self.axs.flat:
            if ax.get_visible() and ax != self.axs[2, 2]:
                for loc in self.rR_locii:
                    ax.axvline(loc, linestyle="--", color="black", linewidth=1)
                if ax.lines:
                    labelLines(ax.get_lines(), zorder=3)

        # Plot worst-case in final subplot
        self._plot_worst_case_envelope()

        # ───── Uniform Y-Limits Across Rows ─────
        for row_idx in range(3):
            row_axes = [self.axs[row_idx, col] for col in range(3) if self.axs[row_idx, col].get_visible()]
            y_mins, y_maxs = [], []
            for ax in row_axes:
                for line in ax.get_lines():
                    y = line.get_ydata()
                    y_mins.append(np.nanmin(y))
                    y_maxs.append(np.nanmax(y))
            if y_mins and y_maxs:
                y_min, y_max = min(y_mins), max(y_maxs)
                padding = 0.05 * (y_max - y_min) if (y_max - y_min) > 0 else 1.0
                for ax in row_axes:
                    ax.set_ylim(y_min - padding, y_max + padding)

        # Final formatting and save
        self.fig.suptitle("Internal Force Equilibrium Diagram", fontsize=15)
        self.fig.tight_layout(rect=[0, 0, 1, 0.95])

        out_path = self.output_dir / filename
        self.fig.savefig(out_path, dpi=300)
        print(f"[SAVED] Plot saved to: {out_path}")