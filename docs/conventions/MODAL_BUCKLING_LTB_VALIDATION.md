# Modal buckling: lateral–torsional buckling (LTB) and warping χ — validation playbook

This note complements [`JOB_INPUT_BEAM_WARPING.md`](JOB_INPUT_BEAM_WARPING.md) (buckling subsection): the code path embeds a **12×12** beam-column geometric stiffness **K_σ** on translation/rotation DOFs and leaves **χ** rows/cols zero in **K_σ**. That is appropriate for **in-plane column buckling** about a prestressed line model; it does **not** fully represent thin-walled **lateral–torsional buckling** where warping and bending–torsion coupling in **K_σ** matter.

## When to treat LTB as out of model scope

- Open or mono-symmetric sections with significant **Γ** and loads that drive **minor-axis** bending or **torsional** response.
- Any design limit state explicitly framed as **LTB** or **distortional** modes on the thin-walled cross-section.

In those cases, treat this 1D beam-line buckling eigenvalue as **non-conservative** unless you validate against a richer model.

## Practical checklist

1. **Clarify the limit state**  
   Column flexural buckling about a principal axis → current **K + λK_σ** workflow may suffice after mesh refinement. True LTB → require external validation.

2. **Conservative screening run**  
   Run buckling with **`[warping]` off** (or a 6 DOF/node mesh) for a simple Euler-type check when the dominant mode is column-like; see [`JOB_INPUT_BEAM_WARPING.md`](JOB_INPUT_BEAM_WARPING.md).

3. **Thin-walled section modelling**  
   Review [`BEAM_SHEAR_CORRECTION_AND_THINWALL.md`](BEAM_SHEAR_CORRECTION_AND_THINWALL.md): scalar **κGA** and line resultants do not replace plate/shell shear flow. Non-zero **y_sc**, **z_sc**, and **Γ** improve fidelity but do not, by themselves, fix missing **K_σ** coupling terms for LTB.

4. **External reference**  
   - Handbook or code formula (e.g. idealised simply supported I-beam **M_cr**) with **matching boundary conditions**.  
   - Higher-fidelity FE (shell/solid or beam with dedicated LTB DOFs in another code).  
   - Optional workflow: [`post_processing/validation_visualisers/README_VALIDATION_ABAQUS.md`](../../post_processing/validation_visualisers/README_VALIDATION_ABAQUS.md) for comparing critical loads or mode shapes when exporting reference solutions.

5. **In-repo regressions**  
   [`tests/test_modal_buckling_warping_benchmark.py`](../../tests/test_modal_buckling_warping_benchmark.py) checks embedding of **K_σ** on the first 12 DOFs and a warping-off Euler band — **not** a thin-walled LTB benchmark.

## Closed-form or numerical benchmark test (future)

A single regression test against a published **P_cr** or **M_cr** is reasonable **only after** the reference geometry, material, BCs, and sign conventions are fixed. Until then, document the comparison in job-specific validation notes rather than hard-coding uncertain tolerances.

### Repository policy (deferral)

No in-repository closed-form **P_cr** / **M_cr** / eigenvalue benchmark is added until a reference solution (handbook, code formula, or independent FE) is agreed with matching BCs and load definition. Until then, track comparisons in job-specific validation notes or external reports rather than committing tolerance-bound assertions that may drift with undocumented assumptions.

### In-repo regression template (when reference is fixed)

When promoting a comparison to a permanent test (e.g. under `tests/test_modal_buckling_ltb_reference.py`):

1. **Pin the job** — commit or document a `jobs/job_*` directory (mesh, section, material, loads, buckling flags) that reproduces the reference case; cite the handbook equation or FE report in the test docstring.
2. **State the scalar** — assert **λ**, **P_cr**, or **M_cr** with units and whether the reference is elastic eigenvalue or another convention.
3. **Tolerance** — use a tight relative tolerance only if the reference digits are known; otherwise document sensitivity (mesh **n**, warping on/off) in the test module docstring.
4. **Scope** — keep the test focused on the agreed BC class; do not generalise to arbitrary thin-walled LTB without extending **K_σ** modelling.

Existing smoke coverage remains [`tests/test_modal_buckling_warping_benchmark.py`](../../tests/test_modal_buckling_warping_benchmark.py) (embedding / Euler band), not LTB closure.
