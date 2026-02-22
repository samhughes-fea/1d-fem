# Job layout

Summary of what each job is for. Use this when interpreting `post_processing/results` and when running validation visualisers (e.g. Roark vs FEA, deformation convergence).

## Point load – Euler–Bernoulli, three positions (job_0000, job_0001, job_0002)

Same cantilever, Euler–Bernoulli element; point load position varies (end, midspan, quarter).

| Job       | Load position   | Formula / x   |
|-----------|------------------|---------------|
| job_0000  | End (tip)        | P(x=L)        |
| job_0001  | Midspan          | P(x=L/2)      |
| job_0002  | Quarter-point    | P(x=L/4)      |

Use these to validate Euler–Bernoulli tip/midspan/quarter deflection vs analytical (e.g. Roark point-load formulas with different `a`).

## End load – different element types (job_0003, job_0004)

Same cantilever with tip (end) point load; only the beam element type changes.

| Job       | Element type              | Load        |
|-----------|---------------------------|------------|
| job_0003  | TimoshenkoBeamElement3D   | End (tip) point load |
| job_0004  | LevinsonBeamElement3D     | End (tip) point load |

Use these to compare element formulations (Timoshenko vs Levinson) under the same end load.

## Truss and Bar elements

| Element type        | Carries                    | Use case                          |
|---------------------|----------------------------|-----------------------------------|
| TrussElement3D      | Axial, transverse, torsion | Pin-jointed + shear + torsion     |
| BarElement3D        | Axial, torsion only        | Axial + torsion, no transverse    |

Both use 6 DOF per node and return 12×12 K_e, 12×1 F_e. Register in the [Element] file as `TrussElement3D` or `BarElement3D`.

**Optional job**: `job_bar_single` – one Bar element, fixed at node 0, axial point load at node 1; for pipeline check.

## Distributed loads – Euler–Bernoulli (job_0005, job_0006, job_0007)

Same cantilever, Euler–Bernoulli element; load type varies (UDL, triangular, parabolic).

| Job       | Load type   | Formula / profile        |
|-----------|-------------|---------------------------|
| job_0005  | UDL         | q(x) = w (constant)       |
| job_0006  | Triangular  | q(x) = w·(x/L)            |
| job_0007  | Parabolic   | q(x) = w·(x/L)²           |

Use these for Roark vs FEA validation of distributed-load formulas (e.g. `post_processing/validation_visualisers/roarks_formulas`).

## jobs/bin

Duplicate/canonical definitions for some of the above (e.g. UDL, triangular, parabolic) also live under `jobs/bin` (job_0004, job_0005, job_0006 there map to UDL, Triangular, Parabolic). Validation scripts may reference either `jobs/` or `jobs/bin` for job input files.
