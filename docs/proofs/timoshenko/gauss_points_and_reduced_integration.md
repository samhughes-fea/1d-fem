# Timoshenko: Gauss points and reduced integration in post-processing

## Stiffness vs post-processing

- **Stiffness assembly** uses **selective/reduced integration**: the shear block (rows 3,4 of B) is integrated with **1-point** quadrature to avoid shear locking; bending and the rest use full order (e.g. 3 points).
- **Formulation cache** (and thus strain, stress, section forces) is built from the **full-order** loop over `xi_full = leggauss(max_order)`, so there are **3 Gauss points per element** (when `max_order = 3`).
- Post-processing (secondary/tertiary results) therefore evaluates **B @ U_e** and **D @ ε** at those **same 3 points**. Reduced integration is only used when assembling K_e; it does not change the number or positions of the cached GPs.

## Why three GPs at (nearly) constant value?

For a 2-node Timoshenko element under **constant shear** (e.g. cantilever with tip load), the shear force Vy is constant along the element. So at all three Gauss points we get the same Vy. That is correct: the "three Gauss points per element at a constant value" for Vy is expected.

## Gauss point placement in the plot

GP positions are **3-point Gauss–Legendre** in natural coordinates: ξ ≈ −0.775, 0, 0.775. In physical space they are **not** equally spaced (they cluster slightly toward the element ends). The section forces visualiser infers xi from the number of rows in the CSV (`n_gp`) and uses `leggauss(n_gp)[0]`, so the plotted positions match the formulation cache order (ascending ξ).

## Expected magnitude (e.g. job_0003)

For a cantilever with tip load **F_y = −500 N**, internal shear should be **|Vy| ≈ 500 N** (constant along the beam). If the section forces plot shows a very different scale (e.g. ~20 000 N), check:

1. **Re-run the simulation** – Ensure tertiary results (and section force CSVs) are regenerated so that files start with `# column_order=resultant` and columns are [N, Vy, Vz, T, My, Mz]. Legacy CSVs are reordered on read, but re-running is the reliable fix.
2. **Point load magnitude** – In `jobs/job_0003/point_load.txt` confirm F_y is as intended (e.g. −500 N).
3. **Diagnostic** – Run the section forces visualiser with `DEBUG_SECTION_FORCES=1` (e.g. `set DEBUG_SECTION_FORCES=1` then run the script); it will print min/max per component so you can verify Vy is on the order of 500 N.

## References in code

- Timoshenko stiffness and cache: `pre_processing/element_library/timoshenko/timoshenko_3D.py` (element_stiffness_matrix: gauss_cache from xi_full loop; shear_order = 1 for Ke only).
- Section forces visualiser GP positions: `post_processing/.../section_forces/section_forces_visualisation.py` (GAUSS_3PT_XI, _gauss_point_x_for_element).
