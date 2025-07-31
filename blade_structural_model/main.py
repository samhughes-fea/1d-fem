from pathlib import Path
from action.blade_internal_action_runner import BladeInternalActionRunner

R = 0.80
L = 0.70
TSR_NAMES = ["TSR4", "TSR5", "TSR6", "TSR7", "TSR8"]

COLORS = {
    "load":    "#4F81BD",
    "shear":   "#9BBB59",
    "bending": "#C0504D",
    "torsion": "#8064A2",
}

SCRIPT_DIR = Path(__file__).resolve().parent
LOAD_DIR = SCRIPT_DIR / "load_profiles"
OUT_DIR = SCRIPT_DIR / "outputs"

runner = BladeInternalActionRunner(
    load_dir=LOAD_DIR,
    out_dir=OUT_DIR,
    tsr_names=TSR_NAMES,
    R=R,
    L=L,
    colors=COLORS
)

runner.run()
