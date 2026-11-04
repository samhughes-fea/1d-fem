# Abaqus results: rotation_source and run_log summary

Checked all job dirs under `abaqus_results/`.

## rotation_source.txt

| Job              | rotation_source |
|------------------|-----------------|
| All 41 jobs      | **ODB**         |

So for every job, rotation (θ_z) on the deformation plots comes from **Abaqus UR in the ODB** (with our sign flip applied), not from our derived θ_z = −du_y/dx.

## run_log.txt (consistent across jobs)

- **Field output request:** Fails (Model / StaticStep has no `FieldOutputRequest` on this Abaqus version). Messages refer to (U, UR, SF) or (U, UR, SF, SM) depending on when the script was generated.
- **ODB field outputs:** Always include `'U'` and `'UR'` (Abaqus default output). Some jobs also show `'CF', 'CM'`; job_0005_* sometimes show `['E', 'RF', 'RM', 'S', 'U', 'UR']` (no CF, CM).
- **ODB has UR:** **True** for every job.

So although the explicit field output request fails, the ODB still contains U and UR, and the extract script correctly reads them and writes `rotation_source.txt` = ODB.

## Jobs with multiple run_log blocks

Some dirs (e.g. job_0000_n10, job_0001_n10, job_0007_n8) have run_log content from several runs (older "U, UR, SF" failures and newer "U, UR, SF, SM" plus the note about default U/UR). The last run’s ODB state is what the current U_global.csv and rotation_source reflect.
