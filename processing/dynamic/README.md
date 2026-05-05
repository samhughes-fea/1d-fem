# `processing.dynamic`

Deprecated compatibility package for Section 3 transient dynamics.

Use [`processing.transient`](../transient/README.md) as the canonical package name.

## Staged operations

| Class | File | Role |
|-------|------|------|
| `AssembleTransientGlobalSystem` | [`../transient/operations/assemble_transient_global.py`](../transient/operations/assemble_transient_global.py) | Global **K**, **M**, optional element **C** |
| `ModifyTransientGlobalSystem` | [`../transient/operations/modify_transient_global.py`](../transient/operations/modify_transient_global.py) | BCs on **K**, **M**, **C** |
| `IntegrateTransientSystem` | [`operations/integrate_transient_system.py`](operations/integrate_transient_system.py) | Newmark wrapper |

[`TransientSimulationRunner`](../../simulation_runner/transient/dynamic_simulation.py) composes these with **`RuntimeMonitorTelemetry`** using the same staged-method pattern as linear static.

## Entry-point convention

Operation classes use **`run(...)`** as the main entry point (see [`processing/spectral/README.md`](../spectral/README.md) for contrast with static `processing.static.operations` naming).
