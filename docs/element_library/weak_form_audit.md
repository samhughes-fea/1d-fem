# Weak-form Gauss audit (`pre_processing/element_library`)

Reference: [FORMULATION_DOCSTRING_STANDARDS.md](../conventions/FORMULATION_DOCSTRING_STANDARDS.md) § *Weak-form assembly (required)*.

| Element module | `K_e` | `F_dist` | `F_point` | `M_e` | `F_int` | `K_mat` | `K_σ` | Notes |
|----------------|-------|----------|-----------|-------|---------|---------|-------|-------|
| `linear_euler_bernoulli_3D` | Gauss | Gauss | `N(ξ_p)ᵀ P` | Gauss | — | — | — | |
| `linear_timoshenko_3D` | Gauss (selective) | Gauss | same | Gauss | — | — | — | Two Gauss-style blocks for bend/shear |
| `linear_euler_bernoulli_3D` + warping (`euler_bernoulli_with_warp/utilities/`) | Gauss | Gauss | same | Gauss | — | — | — | 7 DOF when `[warping]` on; thin `LinearWarping*` alias removed |
| `linear_timoshenko_3D` + warping | Gauss | Gauss | same | Gauss | — | — | — | Same; use baseline type + `[warping]` |
| `linear_levinson_3D` | Gauss (selective) | Gauss | same | Gauss | — | — | — | |
| `linear_reddy_3D` | Gauss (selective) | Gauss | same | Gauss | — | — | — | |
| `linear_bar_3D` | Gauss | Gauss | same | Gauss | — | — | — | Constant `B` → quadrature still applied |
| `linear_truss_3D` | Gauss | Gauss | same | Gauss | — | — | — | Same |
| `nonlinear_euler_bernoulli_3D` | — | Gauss | same | — | Gauss | Gauss | **Gauss** | `K_σ`: sum `w_g |J| (N_g + M_z/L) h′h′ᵀ` (xy) + `(N_g + M_y/L)` (xz); axial `u_x` term per GP |
| `nonlinear_timoshenko_3D` | — | Gauss | same | Gauss ref. | Gauss | `assemble_timoshenko_K0` | **Gauss** | `K_0` selective quadrature via `TimoshenkoQuadratureOrders`; TL loops use `loop_order` |
| `nonlinear_timoshenko_3D` + warping | — | Gauss | first 12 only | Gauss | Gauss | full Gauss on `B` (7×14) | **Gauss** (12×12 block) | Same class + `[warping]`: `K_0` from full `B_7x14` (not selective); `D[6,6]=E·Γ_eff` via `beam_warping_policy`; at `U=0`, `K_T` equals linear Timoshenko `K_e` (see `tests/test_warping_nl_linear_material_agreement.py`) |
| `nonlinear_euler_bernoulli_3D` + warping | — | Gauss | first 12 only | Gauss | Gauss | Gauss on `B_7x14` | **Gauss** (12×12 block) | At `U=0`, material tangent matches linear EB warping `K_e`; same Γ policy as Timoshenko |
| `GeometricallyExactShearDeformableBeam3D` | — | Gauss | same as parent | Gauss | Gauss | same as NL Timoshenko | same | Subclasses TL Timoshenko; warping policy identical |

**Resolved gap (historical):** Prior to the weak-form pass, `geometric_stiffness.py` (EB/Timoshenko) used a closed-form 4×4 beam-column template times `N/(30L)` and `M/(30L²)`. That matrix is **equivalent** to \(\sum_g w_g |J| (N_g + M/L) (\partial h/\partial x)(\partial h/\partial x)^\top\) for the EB Hermite bending DOFs (with the xz-plane sign convention on \(\theta_y\) slopes). Implementation now uses the **Gauss sum** explicitly.

**Utilities:** `B_matrix` / `D_matrix` / `shape_functions` / `interpolate_loads` should state which integral they support (`ε = B U_e`, `S = D ε`, `N` for loads).

**Modal linear buckling (global assembly):** Not an element Gauss row — [`assemble_global_geometric_stiffness`](../../processing/modal/buckling.py) builds global **K_σ** from each beam element’s `linear_geometric_stiffness_matrix(U_e)` after a linear static prestress solve. EB/Timoshenko smoke: [`tests/test_modal_buckling_euler_column.py`](../../tests/test_modal_buckling_euler_column.py).
