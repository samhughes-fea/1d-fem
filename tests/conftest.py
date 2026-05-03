# Silence deprecated [Modal] / legacy [Type] parser warnings unless a test clears this.
import os

os.environ.setdefault("FEM_SILENCE_LEGACY_SIMULATION_SETTINGS_WARNINGS", "1")
