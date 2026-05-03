# `processing.harmonic`

Frequency-domain structural response (**-ω² M + i ω C + K**) and related kernels.

## Staged operations

| Class | Role |
|-------|------|
| `AssembleHarmonicStructuralMatrices` | Global **K**, **M** |
| `AssembleHarmonicLoadVector` | Global **F** |
| `ModifyHarmonicStructuralMatrices` | BCs on **K**, **M**, **F** |
| `BuildHarmonicDampingMatrix` | **C** from ζ and optional Rayleigh terms |
| `SolveHarmonicFrequencySweep` | Direct sweep or modal superposition |

[`HarmonicSimulationRunner`](../../simulation_runner/harmonic/harmonic_simulation.py) wires these stages with telemetry.

## Entry-point convention

Stage classes expose **`run(...)`**; see [`processing/spectral/README.md`](../spectral/README.md) for the contrast with **`processing.static.operations`**.
