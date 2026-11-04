# `processing.transient`

Canonical Section 3 processing package for transient dynamics.

This package supersedes [`processing.dynamic`](../dynamic/README.md) as the canonical import path while preserving the same transient dynamics functionality.

## Staged operations

| Class | File | Role |
|-------|------|------|
| `AssembleTransientGlobalSystem` | [`operations/assemble_transient_global.py`](operations/assemble_transient_global.py) | Global **K**, **M**, optional element **C** |
| `ModifyTransientGlobalSystem` | [`operations/modify_transient_global.py`](operations/modify_transient_global.py) | BCs on **K**, **M**, **C** |
| `IntegrateTransientSystem` | [`operations/integrate_transient_system.py`](operations/integrate_transient_system.py) | Newmark wrapper |

Deprecated compatibility imports remain under [`processing.dynamic`](../dynamic/__init__.py).
