# Design note: nonlinear prestress for linear buckling

**Implemented:** `buckling.buckling_prestress = nonlinear_static` (legacy: `modal.buckling_prestress`; see [`BucklingSimulationRunner`](../../simulation_runner/buckling/buckling_simulation.py) / [`VibrationBucklingBackend`](../../simulation_runner/spectral/vibration_buckling_backend.py)).

This document records the **API and architecture** for using a converged **nonlinear static** state to build **K_σ** for the linearised buckling eigenproblem, in addition to the long-standing **`linear_static`** prestress.

## Current behaviour

- **`buckling_prestress = linear_static`** (default): shared buckling pipeline in [`BucklingSimulationRunner`](../../simulation_runner/buckling/buckling_simulation.py) / [`VibrationBucklingBackend`](../../simulation_runner/spectral/vibration_buckling_backend.py) runs [`LinearStaticSimulationRunner.solve_linear_system_only`](../../simulation_runner/static/linear_static_simulation.py) with load scale `buckling_load_factor`, uses **U_global** to evaluate element **K_σ**, then solves \((\mathbf{K} + \lambda \mathbf{K}_\sigma)\boldsymbol{\phi} = \mathbf{0}\).
- **`buckling_prestress = none`**: rejected with **`ValueError`** — no reference stresses → **K_σ** cannot be defined meaningfully.

## Intended nonlinear workflow (design)

1. Run [`NonlinearStaticSimulationRunner`](../../simulation_runner/static/nonlinear_static_simulation.py) with the same mesh, loads (scaled), and BCs **until convergence** per increment.
2. Extract **U_global** at the last converged state (or user-selected increment).
3. Evaluate **K_σ** from each element’s `linear_geometric_stiffness_matrix(U_e)` **using stresses consistent with that nonlinear configuration**. Note: TL/nonlinear elements already define stress from **Green–Lagrange** strains; linear buckling theory assumes **small** disturbances about a **prestressed** state — engineering acceptance often uses the **same** **K_σ** Gauss assembly as today but with **U** from nonlinear equilibrium (Koiter-type approximate buckling about a nonlinear path).

4. Clarify in settings: **load scaling** (does `buckling_load_factor` multiply both nonlinear prestress and eigenvalue interpretation?), **follow-up** increments (single vs arc-length branch), and **equilibrium tolerance** (nonlinear residual before accepting **U** for **K_σ**).

## Element stack

### Option A (default)

**Nonlinear beam types:** The mesh `element.txt` lists nonlinear beam elements (e.g. `NonlinearTimoshenkoBeamElement3D`) that implement `tangent_stiffness_matrix` and `internal_force_vector` for the prestress run.

### Option B — `buckling_nonlinear_prestress_twins = true` (legacy: `modal.buckling_nonlinear_prestress_twins`)

When **`buckling_prestress = nonlinear_static`** and **`buckling_nonlinear_prestress_twins`** is enabled, the buckling runner builds a **parallel** element array with registered nonlinear twins (`LinearTimoshenkoBeamElement3D` → `NonlinearTimoshenkoBeamElement3D`, `LinearEulerBernoulliBeamElement3D` → `NonlinearEulerBernoulliBeamElement3D`), runs `NonlinearStaticSimulationRunner` on those instances, then assembles **`K_σ`** with the **original** linear elements from `run_job` (same global DOF layout).

`assemble_global_geometric_stiffness` still calls `linear_geometric_stiffness_matrix(U_e)` on each **linear** element; nonlinear TL classes used only for prestress delegate `linear_geometric_stiffness_matrix` to a temporary linear twin when they appear in the mesh.

## Implementation checklist

- [x] `modal.buckling_prestress` value `nonlinear_static`, branch in `_run_linear_buckling`.
- [x] `NonlinearStaticSimulationRunner` for prestress; `buckling_load_factor` scales the reference load table for that solve; `nonlinear.num_increments` / `load_factors` may be forced for a single-load prestress as needed.
- [x] Optional twins: `modal.buckling_nonlinear_prestress_twins` with nonlinear twin instantiation ([`VibrationBucklingBackend`](../../simulation_runner/spectral/vibration_buckling_backend.py)).
- [x] Tests: `tests/test_modal_buckling_nonlinear_prestress.py`.
