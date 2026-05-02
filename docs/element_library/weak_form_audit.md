# Weak-form Gauss audit (`pre_processing/element_library`)

Reference: [FORMULATION_DOCSTRING_STANDARDS.md](../conventions/FORMULATION_DOCSTRING_STANDARDS.md) § *Weak-form assembly (required)*.

| Element module | `K_e` | `F_dist` | `F_point` | `M_e` | `F_int` | `K_mat` | `K_σ` | Notes |
|----------------|-------|----------|-----------|-------|---------|---------|-------|-------|
| `linear_euler_bernoulli_3D` | Gauss | Gauss | `N(ξ_p)ᵀ P` | Gauss | — | — | — | |
| `linear_timoshenko_3D` | Gauss (selective) | Gauss | same | Gauss | — | — | — | Two Gauss-style blocks for bend/shear |
| `linear_warping_euler_bernoulli_3D` (`euler_bernoulli_with_warp/`) | Gauss | Gauss | same | Gauss | — | — | — | |
| `linear_warping_timoshenko_3D` | Gauss | Gauss | same | Gauss | — | — | — | |
| `linear_levinson_3D` | Gauss (selective) | Gauss | same | Gauss | — | — | — | |
| `linear_reddy_3D` | Gauss (selective) | Gauss | same | Gauss | — | — | — | |
| `linear_curved_timoshenko_3D` | Gauss | Gauss | same | Gauss | — | — | — | |
| `linear_bar_3D` | Gauss | Gauss | same | Gauss | — | — | — | Constant `B` → quadrature still applied |
| `linear_truss_3D` | Gauss | Gauss | same | Gauss | — | — | — | Same |
| `nonlinear_euler_bernoulli_3D` | — | Gauss | same | — | Gauss | Gauss | **Gauss** | `K_σ`: sum `w_g |J| (N_g + M_z/L) h′h′ᵀ` (xy) + `(N_g + M_y/L)` (xz); axial `u_x` term per GP |
| `nonlinear_timoshenko_3D` | — | Gauss | same | Gauss ref. | Gauss | `assemble_timoshenko_K0` | **Gauss** | `K_0` selective quadrature via `TimoshenkoQuadratureOrders`; TL loops use `loop_order` |
| `geometrically_exact_shear_deformable_beam_3D` | — | — | — | — | — | — | — | Stub (NotImplementedError) |

**Resolved gap (historical):** Prior to the weak-form pass, `geometric_stiffness.py` (EB/Timoshenko) used a closed-form 4×4 beam-column template times `N/(30L)` and `M/(30L²)`. That matrix is **equivalent** to \(\sum_g w_g |J| (N_g + M/L) (\partial h/\partial x)(\partial h/\partial x)^\top\) for the EB Hermite bending DOFs (with the xz-plane sign convention on \(\theta_y\) slopes). Implementation now uses the **Gauss sum** explicitly.

**Utilities:** `B_matrix` / `D_matrix` / `shape_functions` / `interpolate_loads` should state which integral they support (`ε = B U_e`, `S = D ε`, `N` for loads).
