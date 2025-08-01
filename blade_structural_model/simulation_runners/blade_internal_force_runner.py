# blade_structural_model/simulation_runners/blade_internal_force_runner.py

from pathlib import Path
from blade_structural_model.containers.internal_force.blade_internal_force_row import BladeInternalForceRow
from blade_structural_model.containers.internal_force.blade_internal_force_bin import BladeInternalForceBin
from internal_force.blade_load_profile_reader import LoadProfileReader
from internal_force.blade_internal_force_processor import InternalForceProcessor


class BladeInternalForceRunner:
    def __init__(self, load_dir: Path, out_dir: Path, tsr_names: list, R: float, L: float):
        self.load_dir = load_dir
        self.out_dir = out_dir
        self.tsr_names = tsr_names
        self.R = R
        self.L = L
        self.result_bin = BladeInternalForceBin()
        self.rR_locii = None  # Set from first profile

    def run(self) -> BladeInternalForceBin:
        for i, tsr in enumerate(self.tsr_names):
            path = self.load_dir / f"{tsr}.csv"
            if not path.is_file():
                print(f"[WARNING] Skipping {tsr}: file not found → {path}")
                continue

            profile = LoadProfileReader(path, self.R, self.L)
            profile.read()

            if self.rR_locii is None:
                self.rR_locii = profile.rR_locii

            processor = InternalForceProcessor(profile.data, profile.dx)
            processor.compute()

            row = BladeInternalForceRow(
                tsr=tsr,
                r_over_R=profile.rR,
                f_y=profile.data["F_y"],
                f_z=profile.data["F_z"],
                m_x=profile.data["M_x"],
                V_y=processor.results["V_y"],
                M_z=processor.results["M_z"],
                V_z=processor.results["V_z"],
                M_y=processor.results["M_y"],
                T=processor.results["T"],
            )

            self.result_bin.rows.append(row)

        print(f"[INFO] Finished assembling internal force results for {len(self.result_bin.rows)} TSRs")
        return self.result_bin