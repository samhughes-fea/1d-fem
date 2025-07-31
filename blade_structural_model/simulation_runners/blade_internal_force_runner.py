from pathlib import Path
import pandas as pd

from blade_structural_model.containers.internal_force.blade_internal_force_row import BladeInternalForceRow
from blade_structural_model.containers.internal_force.blade_internal_force_bin import BladeInternalForceBin
from internal_force.blade_load_profile_reader import LoadProfileReader
from internal_force.blade_internal_force_processor import InternalForceProcessor
from internal_force.blade_internal_force_visualisation_module import InternalForceVisualisationModule

class BladeInternalForceRunner:
    def __init__(self, load_dir: Path, out_dir: Path, tsr_names: list, R: float, L: float, colors: dict):
        self.load_dir = load_dir
        self.out_dir = out_dir
        self.tsr_names = tsr_names
        self.R = R
        self.L = L
        self.colors = colors
        self.result_bin = BladeInternalForceBin()

    def run(self):
        plotter = InternalForceVisualisationModule(self.tsr_names, self.colors)
        self.out_dir.mkdir(parents=True, exist_ok=True)

        for tsr in self.tsr_names:
            path = self.load_dir / f"{tsr}.csv"
            if not path.is_file():
                print(f"[WARNING] Skipping {tsr}: file not found → {path}")
                continue

            profile = LoadProfileReader(path, self.R, self.L)
            profile.read()

            calc = InternalForceProcessor(profile.data, profile.dx)
            calc.compute()

            row = BladeInternalForceRow(
                tsr=tsr,
                r_over_R=profile.rR,
                f_y=profile.data["F_y"],
                f_z=profile.data["F_z"],
                m_x=profile.data["M_x"],
                V_y=calc.results["V_y"],
                M_z=calc.results["M_z"],
                V_z=calc.results["V_z"],
                M_y=calc.results["M_y"],
                T=calc.results["T"],
            )

            self.result_bin.rows.append(row)
            plotter.plot_single(tsr, row.r_over_R, profile.data, calc.results)
            row.to_dataframe().to_csv(self.out_dir / f"{tsr}_processed.csv", index=False)

        plotter.finalize(self.out_dir / "blade_internal_actions.png")
        print(f"[DONE] Figures and CSVs written to: {self.out_dir}")