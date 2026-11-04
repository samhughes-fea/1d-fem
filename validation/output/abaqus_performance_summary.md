# Abaqus results review summary

- **Result directories scanned:** 40
- **Completed successfully:** 40
- **Aborted / not completed:** 0
- **No .sta file:** 0
- **Total solver time (sec):** 40.00

Validation performance (FEM vs Abaqus agreement) is produced by existing scripts:
- `deflection_tables/deformation_comparison.py`, `section_forces/section_forces_comparison.py`
- `deflection_tables/gci_richardson_abaqus_report.py` → `output/gci_richardson_abaqus_deflection_rotation.csv`
- See `output/` for overlay plots and error CSVs.

## Expected jobs with no result directory

The following jobs exist under `jobs/` but have no `abaqus_results/job_XXXX_nN/` dir:

- job_0004_n10
