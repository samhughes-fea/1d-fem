# Problem: Section forces nodal value spike at support

## Summary

In the section forces plot (Gaussian resolution, with nodal markers), the **nodal** value of Vy at the **support** (root) can show a large spike (e.g. −20 000 N) while:

- Vy along the beam (Gauss point markers and interpolant) is approximately −500 N.
- The **integrated** section forces plot (element-wise mean) shows Vy constant at −500 N.

So the underlying GP data and element means are consistent with the expected shear (~500 N); the error is confined to **how nodal values are obtained** for the section-forces visualisation at boundary nodes (e.g. the fixed support).

## Observed behaviour

- **Job 0003** (Timoshenko, cantilever, tip load F_y = −500 N): Vy should be |Vy| ≈ 500 N everywhere.
- **Integrated section forces plot**: Vy = −500 N constant (correct).
- **Section forces plot (Gauss + nodal)**: Gauss markers and profile show ~−500 N; the **nodal marker at the support** shows ~−20 000 N.

## What was tried (not root cause fixes)

1. **Clipping** nodal values to the min/max of contributing GP values — removes the spike in the plot but hides the symptom; does not fix why the extrapolation is wrong.
2. **Replacing least-squares extrapolation** with the mean of element-wise mean of GP values at each node — gives stable, in-range nodal values but changes the definition of “nodal” and is still a workaround rather than a correction of the extrapolation.

## Shape functions: interpolation vs recovering nodal values

**Interpolation** (what you’re thinking of): we have **nodal** values \(f_1, f_2\). We use the shape functions to get the value at any \(\xi\):
\[
f(\xi) = N_1(\xi)\,f_1 + N_2(\xi)\,f_2.
\]
So we interpolate *from* nodes *to* any point (e.g. Gauss points or along the element). The shape functions are used in this **forward** direction.

**What we have for section forces:** The solution gives us values **at the Gauss points** (e.g. 3 points per element), not at the nodes. So we have \(f_{\mathrm{gp}} = [f(\xi_1), f(\xi_2), f(\xi_3)]\). We still assume the same model: the field is the shape-function interpolant of *some* nodal values \(f_1, f_2\). So
\[
f(\xi_i) = N_1(\xi_i)\,f_1 + N_2(\xi_i)\,f_2.
\]
Stacking the three Gauss points gives \(\boldsymbol{f}_{\mathrm{gp}} = \mathbf{N}_{\mathrm{mat}} \boldsymbol{f}_{\mathrm{nodal}}\), with \(\mathbf{N}_{\mathrm{mat}}\) of size \(3\times2\). That is **3 equations in 2 unknowns** — overdetermined — so we **solve for** \(\boldsymbol{f}_{\mathrm{nodal}}\) in the least-squares sense: find nodal values such that the shape-function interpolant best fits the GP data. So the shape functions **are** used; we’re just solving the **inverse** problem (GP values → nodal values) instead of the forward one (nodal → anywhere).

The word “extrapolation” in the plan only means that the **nodes** are at \(\xi = \pm 1\), while the 3 GPs are at \(\xi \approx -0.77, 0, 0.77\). So the nodal positions lie *outside* the range of the GP positions — we’re inferring values at the ends from values in the middle. The *method* is still “fit the shape-function model to the GP data” (least-squares), not a different kind of maths.

## Root cause to be established

The section forces visualiser currently derives **nodal** values from GP data by either:

- Least-squares extrapolation: fit nodal values so that the shape-function interpolant matches the GP values (overdetermined when 3 GPs, 2 nodes), then average over elements sharing the node; or  
- (After workaround) element-mean averaging at each node.

The spike at the support suggests one or more of the following:

1. **Extrapolation method**
   - Wrong choice or implementation of shape functions (e.g. linear vs element formulation) for the section-force field.
   - Wrong assignment of extrapolated values to node IDs (e.g. connectivity order vs physical left/right, or vs the order of xi used in the shape matrix).
   - Incorrect or inconsistent order of GP rows in the CSV vs the xi order used when building the shape matrix for the least-squares system (so the “left” GP value is paired with the wrong xi).

2. **Data order / convention**
   - Mismatch between the order in which GPs are written (formulation cache / save) and the order assumed when reading and building the extrapolation matrix (e.g. ascending xi in one place, different order in another).
   - For 2-node elements, whether “first node” in connectivity is always the physical left (smaller x) and how that aligns with xi = −1 vs xi = +1 in the visualiser.

3. **Numerical / conditioning**
   - Least-squares system (e.g. 3×2) ill-conditioned or sensitive at the first element so that the “left” nodal value is strongly amplified (would need checking with actual GP values and N_mat for element 0).

4. **Actual GP value at the support**
   - Only if the GP value in the CSV for the first GP of the first element (nearest the support) were ~−20 000 could the spike be “correct” from data; the integrated plot (constant −500) contradicts that unless the other two GPs of that element compensate (e.g. average still −500). That scenario would point to the **computation** of stress/section force at that GP (formulation or boundary effect), not the visualisation.

## Required outcome

- **Root cause** identified (extrapolation logic, data/convention order, or upstream GP value at the support).
- **Fix** in the appropriate layer (visualisation extrapolation and/or ordering, or formulation/post-processing that writes section forces), not a workaround (no clipping or ad-hoc redefinition of nodal values solely to remove the spike).
- Nodal markers in the section forces plot consistent with GP data and with the integrated section forces (e.g. Vy ≈ −500 N at the support for job 0003).

## Recommended next steps (robust for any load case)

Implement in this order:

1. **Confirm data (5 min)**  
   Run the section forces visualiser with `DEBUG_SECTION_FORCES=1` for job_0003 (and optionally one other job). Check that Vy min/max over all GPs are on the order of 500 N, not 20 000. If a single GP is wrong, fix the upstream section-force computation first; otherwise the bug is in the visualisation.

2. **Lock GP ↔ CSV convention in code**  
   - In the visualiser, add a short comment above where CSV rows are read and `xi_used` is set: “CSV row i = GP at natural coordinate xi_used[i] (ascending).”  
   - Optionally in the saver (`save_tertiary_container._save_section_forces`), write one comment line with the xi values used (e.g. `# xi_per_row=...`) so reader and writer stay in sync.  
   - Keep using the same `xi_used` for both GP positions and any future extrapolation matrix (no second `leggauss` with a different order).

3. **Use boundary-aware nodal values (main fix)**  
   - Keep computing an **element-wise mean** per element (mean of the 3 GP values).  
   - For each node, take the **mean of the element-wise means** of all elements that share that node (current behaviour).  
   - **Exception:** For nodes that have **only one** connected element (boundary nodes), define the nodal value as that single element’s mean. (You already do this implicitly when there’s only one contributor; the important part is to **not** use a different rule at boundaries, e.g. no least-squares extrapolation only at boundaries.)  
   - Result: nodal values are always in the range of the GP data, match the integrated section forces at boundaries, and are stable for any load case (point load, distributed, mixed). No extrapolation at the ends.

4. **Optional: least-squares for interior nodes only**  
   If you want “nodal” to mean “shape-function fit” at **interior** nodes (where two elements meet), you can:  
   - For each element, solve the 3×2 least-squares system (shape matrix at `xi_used` vs GP values) to get two nodal values; assign to the element’s left/right nodes using connectivity (column 0 → node at ξ = −1, column 1 → node at ξ = +1).  
   - For **interior** nodes, average the two extrapolated values from the two adjacent elements.  
   - For **boundary** nodes (one element), do **not** use that extrapolate; use the element mean instead (step 3).  
   - This gives smooth shape-function-based values in the interior and avoids boundary spikes.

5. **Add a simple test**  
   - In a test, build a minimal case: one 2-node element, 3 GP values all equal (e.g. −500 for Vy).  
   - Compute nodal values with your chosen rule (element-mean, or boundary-aware with optional LS for interior).  
   - Assert: both nodal values are approximately −500 (e.g. within 1% or exact for pure mean).  
   - If you add least-squares, add a second test: two elements, constant GP value on both; interior node should be ≈ that value, boundary nodes ≈ that value (no spike).

6. **Remove dead code**  
   - Drop any unused clipping or old extrapolation paths; keep one clear rule (e.g. “nodal = mean of element means, with boundary = single-element mean”) so behaviour is predictable for any load case.

Summary: the robust solution for **any** load case is to **never** extrapolate at boundary nodes — use the element mean there — and to lock the GP order convention. Optionally use least-squares only at interior nodes and element mean at boundaries.

**How stress and strain do nodal projection (shape functions)**  
Stress and strain nodal values are **not** computed in the visualisers. They are computed in the pipeline by **NodalResultProjector** (`processing/static/results/compute_secondary/nodal_result_projector.py`), which writes `nodal_strain.csv` and `nodal_stress.csv`; the visualisers only read those files.

1. **Build the shape matrix from the formulation cache**  
   For each element, `N_mat` has shape `(n_gauss, n_nodes_elem)`. Each row is the **element’s** shape functions at one Gauss point: the projector iterates over `elem_obj.gauss_data` and takes `gp.shape_functions` for each GP. So row order matches the formulation’s GP order (ascending ξ for Timoshenko), and the functions are the element’s own (e.g. Timoshenko/EB), not a separate linear assumption.

2. **Solve for nodal values**  
   It assumes values at GPs = N_mat @ values_nodal. For each component (each column of strain or stress):  
   - If **n_gauss == n_nodes_elem** (e.g. 2×2): `np.linalg.solve(N_mat, values_g[:, j])`.  
   - If **n_gauss > n_nodes_elem** (e.g. 3×2): `np.linalg.lstsq(N_mat, values_g[:, j], rcond=None)[0]`.

3. **Average at shared nodes**  
   For each node, sum the extrapolated values from all elements sharing that node and divide by the number of elements.

So stress and strain use **least-squares (or exact solve) with the element’s own shape functions** from the formulation cache. The section-forces visualiser only has CSV data (no cache), so to match that behaviour you would either (a) add a nodal section-forces projection in the same pipeline as NodalResultProjector (using the same N_mat) and save e.g. `nodal_section_forces.csv`, or (b) in the visualiser build the same N_mat using the same ξ convention and linear shape functions and do the same least-squares solve, then average at shared nodes (with the boundary rule: element mean at boundary nodes).

---

## Suggested solutions (detail)

### 1. Verify and fix GP ↔ CSV row ordering (data convention)

**Goal:** Eliminate any mismatch between the order in which GPs are written (formulation) and the order assumed when reading (visualiser).

- **Writer:** In `timoshenko_3D.py`, `gauss_cache` is built in the loop over `xi_full, w_full = np.polynomial.legendre.leggauss(max_order)` (ascending ξ). Stress/section forces are computed in the same GP order, and `save_tertiary_container` writes `elem_section_forces` row-by-row, so CSV row `i` = GP at `xi_full[i]` (ascending).
- **Reader:** The visualiser uses `xi_used = GAUSS_3PT_XI` (i.e. `leggauss(3)[0]`, also ascending) and pairs `data[i]` with `xi_used[i]` for positions and for any extrapolation matrix.
- **Action:** Add an explicit contract and optional check:
  - Document in code (and in this plan): “CSV row index `i` corresponds to Gauss point at natural coordinate `xi_used[i]` in ascending order.”
  - Optionally: in the visualiser, when `DEBUG_SECTION_FORCES` is set, compute `x_gp` from `xi_used` and from element geometry and assert or log that the first GP of element 0 is at the left (smaller x) and the last at the right (larger x), to catch any future reordering.

### 2. Re-introduce least-squares extrapolation with correct conventions

**Goal:** If you want “nodal” to mean “shape-function extrapolation from GP data” (rather than element-mean), implement it in a way that cannot produce the support spike.

- **Convention:** For each element, build the shape matrix so that row `r` corresponds to `data[r]` at `xi_used[r]`, and column 0 = node at ξ = −1 (left), column 1 = node at ξ = +1 (right). Use `_nodal_shape_matrix_at_xi(xi_used, n_nodes_elem)` (already present) and ensure `node_ids[0]` is the node at the left (smaller x) — which is the case when connectivity is `[n_left, n_right]` and `get_element_node_coords` returns coords in connectivity order.
- **Assignment:** Solve the 3×2 least-squares system per component; assign the solution slot for ξ = −1 to `node_ids[0]` and for ξ = +1 to `node_ids[1]`. Then average over elements sharing each node as now.
- **Safeguards:** (a) Use the same `xi_used` as for GP positions (no separate leggauss call with different ordering). (b) Optionally regularise the 3×2 system (e.g. small Tikhonov) to avoid amplification if the system is ill-conditioned. (c) Add a unit test: constant GP values (e.g. all −500) must yield nodal values ≈ −500 for both nodes.

### 3. Boundary-aware nodal value (recommended if spike persists)

**Goal:** Treat boundary nodes in a physically consistent way so that a single-element contribution does not produce an unreasonable extrapolation.

- **Idea:** For nodes that belong to only one element, least-squares extrapolation from 3 GPs to 2 nodes is well-defined but can amplify small errors near the boundary (extrapolation “outside” the GP cloud). A robust approach is: **at nodes that have only one connected element, use the element-wise mean of the GP values** (same as integrated section forces) instead of the extrapolated nodal value.
- **Rationale:** This is not an ad-hoc clip; it matches the fact that at the boundary we do not have a second element to average with, and the element mean is the quantity that matches the integrated section forces plot. So “nodal value = element mean at boundary” is a consistent definition.
- **Implementation:** When building nodal values, detect nodes with `weight[nid] == 1` (or count of contributing elements == 1). For those nodes, set the value to the mean of the GP values of that single element (per component), instead of using the least-squares extrapolate for that element. Interior nodes can keep least-squares (or current element-mean) as desired.

### 4. Confirm GP values at the support (upstream sanity check)

**Goal:** Rule out that the spike comes from an incorrect section force value at the first GP of the first element.

- **Check:** For job_0003, with `DEBUG_SECTION_FORCES=1`, inspect the printed min/max for Vy. If the GP values (and thus CSV rows) are all on the order of −500 N, the bug is in the visualisation (extrapolation or assignment). If one of the GP values (e.g. first row of element 0) is ~−20 000, the bug is upstream (formulation or post-processing that fills section forces at that GP).
- **Reference:** Integrated section forces are element-wise means; if they show −500 N constant, the mean of the 3 GPs per element is correct, so a single GP could still be wrong only if the other two compensate (e.g. two at +9 500 and one at −20 000). Checking the actual CSV for element 0 (all three Vy values) quickly confirms this.

### 5. Optional: write ξ (or GP index) in the CSV for traceability

**Goal:** Make GP order unambiguous for any future reader or tool.

- **Action:** In `save_tertiary_container._save_section_forces`, optionally write a comment line listing the natural coordinates used for each row, e.g. `# xi_per_row=-0.7745966692414834,0.0,0.7745966692414834`, or add a small header that documents “rows in ascending xi”. Then the visualiser (or a test) can assert that its `xi_used` matches.

---

## References

- Section forces visualisation: `post_processing/graphical_visualisers/tertiary_visualisers/section_forces/section_forces_visualisation.py`
- GP positions and CSV order: formulation cache writes GPs in ascending xi order; visualiser uses `GAUSS_3PT_XI` / `leggauss(n_gp)[0]` for positions and for building the extrapolation matrix.
- Integrated section forces (element mean): `post_processing/.../integrated_section_forces/` (shows correct constant Vy).
- Job 0003: Timoshenko, tip load F_y = −500 N; `jobs/job_0003/`.

## Option A implemented

Nodal section forces are now projected in the tertiary pipeline using the same shape-function extrapolation as stress/strain (formulation cache N_mat), with a boundary rule (single-element nodes use element mean). Implemented components:

- **Projector:** `processing/static/results/compute_tertiary/nodal_section_forces_projector.py` — projects GP section forces to nodes via N_mat (solve/lstsq), averages at shared nodes, and uses element mean at boundary nodes.
- **Storage:** `tertiary_results/nodal/nodal_section_forces.csv` (one row per node, columns N,Vy,Vz,T,My,Mz).
- **Visualiser:** Section-forces visualisation prefers this file when present; otherwise falls back to element-mean from GP CSVs. GP **ξ** for plot positions uses **`# xi_per_row=`** from each element CSV when present (`_read_elem_section_forces_csv` in `section_forces_visualisation.py`).
- **Tests:** `tests/test_nodal_section_forces_projector.py` — constant GP section forces yield matching nodal values (no boundary spike).
