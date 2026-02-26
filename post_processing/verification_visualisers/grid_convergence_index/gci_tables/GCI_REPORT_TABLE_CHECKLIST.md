# GCI report table checklist

For a complete, auditable GCI (Grid Convergence Index) report table, the following are typically required or recommended (Roache 1994; Celik et al.; ASME V&V 20–style reporting).

## Essential (already in this report)

| Item | Description | In this table |
|------|-------------|---------------|
| **Grid sizes** | \(N_1\), \(N_2\), \(N_3\) (or \(h_1,h_2,h_3\)) | In table note: \(n = 16,\,8,\,4\) |
| **Refinement ratio** \(r\) | Ratio between successive grids (e.g. \(r=2\)) | In table note |
| **GCI\_12, GCI\_23** | GCI on fine and medium grids (%) | ✓ Columns |
| **Apparent order** \(p\) | Observed order of convergence | ✓ Column \(p\) |
| **Richardson extrapolate** | \(\phi_\mathrm{ext}\) (Rich) | ✓ Column Rich |
| **Reference solution** | Analytical or benchmark (here: Roark) | ✓ Column Roark |
| **Error vs reference** | e.g. \(\Delta = 100(\phi_\mathrm{ext}-\mathrm{Roark})/\mathrm{Roark}\) (%) | ✓ Column \(\Delta\) (%) |
| **Quantity of interest** | Name and units (e.g. tip def mm, tip rot deg) | ✓ Job, Load type, QoI |

## Recommended (methodology and interpretation)

| Item | Description | In this table |
|------|-------------|---------------|
| **Safety factor** \(F_s\) | Used in GCI (e.g. 1.25 for three-grid) | In table note |
| **Theoretical order** \(p_\mathrm{th}\) | Expected order (e.g. 2 for this scheme); allows asymptotic check \(p \approx p_\mathrm{th}\) | In table note: \(p_\mathrm{th}=2\) |
| **Table note / caption** | States \(r\), \(F_s\), grid levels, and that Roark is reference | ✓ See .tex |

## Optional (for deeper verification) — included in this table

| Item | Description | In this table |
|------|-------------|---------------|
| **Fine-grid solution** \(\phi_1\) | So readers see the discrete solution, not only \(\phi_\mathrm{ext}\). | ✓ Column \(\phi_1\) |
| **Relative error fine vs Rich** | \(100|\phi_1 - \phi_\mathrm{ext}|/|\phi_\mathrm{ext}|\) (%) — shows benefit of extrapolation. | ✓ Column Fine→Rich (%) |
| **GCI consistency ratio** | \(\mathrm{GCI}_{23}/\mathrm{GCI}_{12}\); in asymptotic range should be \(\approx r^p\) (e.g. \(\approx 4\) for \(r=2,\,p=2\)). | ✓ Column GCI₂₃/GCI₁₂ |
| **Solution on all three grids** | \(\phi_1,\phi_2,\phi_3\) — full triple for reproducibility. | \(\phi_1\) only; \(\phi_2,\phi_3\) in CSV |

## Summary

- **Minimum**: Grid info (in note), GCI\_12, GCI\_23, \(p\), Rich, reference, error vs reference, QoI + units.
- **Good practice**: Add note with \(r\), \(F_s\), \(p_\mathrm{th}\), and short definition of reference.
- **Strict / journal**: Consider adding \(\phi_1\), fine-vs-Rich error, and/or GCI ratio for asymptotic check.
