import matplotlib.pyplot as plt
from labellines import labelLines

class VisualisationModule:
    def __init__(self, tsr_names: list, colors: dict):
        self.fig, self.axs = plt.subplots(3, 3, figsize=(13, 10), sharex="col")
        self.tsr_names = tsr_names
        self.colors = colors

    def plot_single(self, tsr: str, rR, loads: dict, resultants: dict):
        ax = self.axs
        ax[0, 0].plot(rR, loads["F_y"], label=tsr, color=self.colors["load"])
        ax[1, 0].plot(rR, resultants["V_y"], label=tsr, color=self.colors["shear"])
        ax[2, 0].plot(rR, resultants["M_z"], label=tsr, color=self.colors["bending"])

        ax[0, 1].plot(rR, loads["F_z"], label=tsr, color=self.colors["load"])
        ax[1, 1].plot(rR, resultants["V_z"], label=tsr, color=self.colors["shear"])
        ax[2, 1].plot(rR, resultants["M_y"], label=tsr, color=self.colors["bending"])

        ax[0, 2].plot(rR, loads["M_x"], label=tsr, color=self.colors["load"])
        ax[1, 2].plot(rR, resultants["T"], label=tsr, color=self.colors["torsion"])

    def finalize(self, out_path):
        ylabels = [
            ("$f_y$ [N/m]", "$f_z$ [N/m]", "$m_x$ [Nm/m]"),
            ("$V_y(x)$ [N]", "$V_z(x)$ [N]", "$T(x)$ [Nm]"),
            ("$M_z(x)$ [Nm]", "$M_y(x)$ [Nm]", ""),
        ]

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

        for ax in self.axs.flat:
            if ax.get_visible() and ax != self.axs[2, 2]:
                ax.axvline(0.0, linestyle="--", color="black", linewidth=1)
                ax.axvline(0.125, linestyle="--", color="black", linewidth=1)
                ax.axvline(1.0, linestyle="--", color="black", linewidth=1)
                if ax.lines:
                    labelLines(ax.get_lines(), zorder=3)

        self.fig.suptitle("Internal Force Equilibrium Diagram – TSR4 – TSR8", fontsize=15)
        self.fig.tight_layout(rect=[0, 0, 1, 0.95])
        self.fig.savefig(out_path, dpi=300)
        plt.show()
