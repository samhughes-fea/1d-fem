# `processing.spectral`

Staged matrix operations for **Section 2** (undamped vibration pencil **K**, **M**) and **Section 5** linear buckling (**K**, **K_g** after prestress).

## Modules

| Module | Class | Role |
|--------|-------|------|
| [`operations/prepare_spectral_local.py`](operations/prepare_spectral_local.py) | `PrepareSpectralLocalMatrices` | COO-format element **K_e**, **M_e** |
| [`operations/assemble_spectral_global.py`](operations/assemble_spectral_global.py) | `AssembleSpectralGlobalSystem` | Global **K**, **M** |
| [`operations/modify_spectral_global.py`](operations/modify_spectral_global.py) | `ModifySpectralGlobalSystem` | Penalty / prescribed BCs |
| [`operations/solve_generalized_eigenproblem.py`](operations/solve_generalized_eigenproblem.py) | `SolveGeneralizedEigenproblem` | Smallest generalized eigenpairs |
| [`operations/buckling_stages.py`](operations/buckling_stages.py) | `AssembleBucklingGeometricStiffness`, `ModifyBucklingGlobalMatrices`, `SolveLinearBucklingEigenpairs` | Buckling pipeline |
| [`spectral_diagnostics.py`](spectral_diagnostics.py) | `log_spectral_diagnostics`, `log_spectral_constrained_dofs` | Matrix / BC logging |

Orchestration and the static-vs-spectral stage table: [`simulation_runner/spectral/README.md`](../../simulation_runner/spectral/README.md).

## Entry-point convention

Stage classes in this package expose a primary **`run(...)`** method (and sometimes a secondary helper). That differs from **`processing.static.operations`**, where linear static uses names such as **`assemble()`**, **`apply_boundary_conditions()`**, **`apply_condensation()`** on purpose-bound classes. New spectral/harmonic/dynamic stages should keep **`run`** for consistency across non-static families unless there is a strong readability reason not to.
