# Timoshenko 3D beam element – formulation review

**Date:** 2026-02-22  
**Scope:** `pre_processing/element_library/timoshenko/` (B matrix, D matrix, shape functions, stiffness/force assembly).

---

## 1. Formulation summary

- **Strain vector:** ε = [ε_x, κ_y, κ_z, γ_xy, γ_xz, φ_x]  
  - κ_y = ∂θ_y/∂x (DOFs 4, 10), κ_z = ∂θ_z/∂x (DOFs 5, 11)  
  - γ_xy = ∂u_y/∂x − θ_z (DOFs 1, 7, 5, 11), γ_xz = ∂u_z/∂x − θ_y (DOFs 2, 8, 4, 10)
- **D matrix:** Diagonal (EA, EI_y, EI_z, κGA, κGA, GJ_t); ordering matches B.
- **Stiffness:** K = K_axial + K_bend + K_shear + K_torsion; bending uses selective integration (from element_array); shear block uses **1-point** reduced integration to avoid locking.
- **Force vector:** F_e from distributed loads (∫ N^T q dx) and point loads (N(x_p)^T P). Translation DOFs [0,1,2,6,7,8] receive F_x,F_y,F_z; rotation DOFs [3,4,5,9,10,11] receive M_x,M_y,M_z. Consistent with 6-DOF/node ordering.

---

## 2. Root cause of wrong results (historical): shear locking

- **Shape functions (before fix):** The Timoshenko element previously used the same shape functions as Euler–Bernoulli (Hermite: cubic u_y, u_z and θ_z, θ_y with **du_y/dx = θ_z** and **du_z/dx = θ_y** at the nodes), so γ_xy, γ_xz ≈ 0 and the element was over-stiff (shear locking).
- **Evidence:** Job 3 (Timoshenko) vs job 0 (EB): same load; Timoshenko u_y was much smaller than EB and θ_z at tip had wrong sign and was noisy.

---

## 3. Implementation (2026-02-22)

1. **Timoshenko-specific shape functions**  
   - **Implemented:** `timoshenko/utilities/shape_functions.py` now uses **linear Lagrange** for u_y, u_z and for θ_z, θ_y so that u and θ are independent and γ = du/dx − θ can be non-zero.

2. **Shear quadrature**  
   - **Implemented:** In `timoshenko_3D.py`, the shear block K_shear is integrated with **1-point** Gauss–Legendre (reduced integration) to avoid locking.

3. **B and D**  
   - Unchanged; strain and material ordering remain correct.

---

## 4. Force vector

- Point and distributed loads are applied via N^T; DOF ordering (translation vs rotation, node 1 vs node 2) is consistent with the 12-DOF element and with the B matrix. No change required for load application once the element uses proper Timoshenko shape functions and quadrature.

---

## 5. References

- B matrix: `timoshenko/utilities/B_matrix.py`
- D matrix: `timoshenko/utilities/D_matrix.py`
- Shape functions: `timoshenko/utilities/shape_functions.py` (Timoshenko linear u_y, u_z, θ_y, θ_z; axial/torsion linear)
- Stiffness/force assembly: `timoshenko/timoshenko_3D.py`
