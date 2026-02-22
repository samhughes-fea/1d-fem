# Tertiary result visualisers

These scripts read saved tertiary results from `post_processing/results/job_*/tertiary_results/` and produce plots. They use the same job-discovery pattern as the primary and secondary visualisers (grid and element data from `jobs/job_{id}/`).

Layout is **analogous to primary and secondary**: one folder per quantity (one script + its `*_plots/` output dir per folder).

- **total_strain_energy/** – `total_strain_energy_visualisation.py` → `total_strain_energy_plots/` (elemental)
- **integrated_section_forces/** – `integrated_section_forces_visualisation.py` → `integrated_section_forces_plots/` (elemental)
- **section_forces/** – `section_forces_visualisation.py` → `section_forces_plots/` (Gaussian)
- **principal_stress/** – `principal_stress_visualisation.py` → `principal_stress_plots/` (Gaussian)
- **mohrs_circle_3d/** – `mohrs_circle_3d_visualisation.py` → `mohrs_circle_3d_plots/` (Gaussian; uses principal_stress data)

## Input paths

- **Elemental:** `job_*/tertiary_results/elemental/total_strain_energy.csv`, `integrated_section_forces.csv`
- **Gaussian:** `job_*/tertiary_results/gaussian/section_forces/section_forces_elem_*.csv`, `principal_stress/principal_stress_elem_*.csv`. Mohr's circle 3D reads the same principal_stress CSVs and plots the three circles for the stress state with max shear.

## GP position inference (Gaussian visualisers)

The saved CSVs do not store Gauss point coordinates. The section-forces and principal-stress visualisers infer GP physical positions from:

- Element geometry (node coordinates from `jobs/job_{id}/grid.txt` and `element.txt`)
- A **fixed quadrature rule**: 3-point Gauss-Legendre (ξ = ±√(3/5), 0), matching the default Euler-Bernoulli 3D integration order.

If the pipeline uses a different element type or quadrature order, GP positions may not align exactly; a future option is to save GP ξ or x in the tertiary export (see plan Option B).

## Plotting convention

Markers follow the B2 convention: small solid circle at Gauss points or element midpoints. See `docs/plans/b2_shape_function_coefficients_eb_implementation_plan.md` and `resolution_plotting_utils.py`.

## Running

From project root (or with `sys.path` set):

```text
python post_processing/graphical_visualisers/tertiary_visualisers/total_strain_energy/total_strain_energy_visualisation.py
python post_processing/graphical_visualisers/tertiary_visualisers/integrated_section_forces/integrated_section_forces_visualisation.py
python post_processing/graphical_visualisers/tertiary_visualisers/section_forces/section_forces_visualisation.py
python post_processing/graphical_visualisers/tertiary_visualisers/principal_stress/principal_stress_visualisation.py
python post_processing/graphical_visualisers/tertiary_visualisers/mohrs_circle_3d/mohrs_circle_3d_visualisation.py
```

Figures are written into each script’s `*_plots/` directory next to the script.

## Deprecation

The following are **deprecated and removed**; do not use:

- **Old layout:** Scripts at the root of `tertiary_visualisers/` or under `elemental/` and `gaussian/` subfolders. Use the quantity-named folders above (e.g. `total_strain_energy/`, `section_forces/`).
- **Old run paths:** e.g. `tertiary_visualisers/total_strain_energy_visualisation.py` or `tertiary_visualisers/elemental/total_strain_energy/...`. Use the paths in the Running section.
