from pathlib import Path
from simulation_runners.blade_internal_force_runner import BladeInternalForceRunner
from visualisation.internal_force.internal_force_visualisation_module import InternalForceVisualisationModule
 
R = 0.80
L = 0.70
TSR_NAMES = ["TSR4", "TSR5", "TSR6", "TSR7", "TSR8"]

SCRIPT_DIR = Path(__file__).resolve().parent
LOAD_DIR = SCRIPT_DIR / "load_profiles"
OUT_DIR = SCRIPT_DIR / "outputs"

runner = BladeInternalForceRunner(
    load_dir=LOAD_DIR,
    out_dir=OUT_DIR,
    tsr_names=TSR_NAMES,
    R=R,
    L=L,
)

runner.run()

# ─── Visualise after data is fully assembled ────
if runner.rR_locii:
    visualiser = InternalForceVisualisationModule(
        tsr_names=tsr_names,
        rR_locii=runner.rR_locii,
        output_dir=output_dir
    )
    visualiser.plot(result_bin)
