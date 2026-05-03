# `processing.dynamic`

Time integration for **M u'' + C u' + K u = F(t)** (Section 3).

## Staged operations

| Class | File | Role |
|-------|------|------|
| `AssembleDynamicGlobalSystem` | [`operations/assemble_dynamic_global.py`](operations/assemble_dynamic_global.py) | Global **K**, **M**, optional element **C** |
| `ModifyDynamicGlobalSystem` | [`operations/modify_dynamic_global.py`](operations/modify_dynamic_global.py) | BCs on **K**, **M**, **C** |
| `IntegrateTransientSystem` | [`operations/integrate_transient_system.py`](operations/integrate_transient_system.py) | Newmark wrapper |

[`TransientSimulationRunner`](../../simulation_runner/transient/dynamic_simulation.py) composes these with **`RuntimeMonitorTelemetry`** using the same staged-method pattern as linear static.

## Entry-point convention

Operation classes use **`run(...)`** as the main entry point (see [`processing/spectral/README.md`](../spectral/README.md) for contrast with static `processing.static.operations` naming).
