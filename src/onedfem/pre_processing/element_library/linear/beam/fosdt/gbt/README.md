# Linear GBT beam (fosdt)

`fosdt` = first-order shear deformation theory. This is a **spanwise 1D** beam: Gauss
integrals run along the element axis. Extra nodal DOFs (`b` warping amplitudes,
Hermite distortional modes `a`) come from **section-derived** GBT enrichment, not
from a 2D section-plane finite element.

Operators live in `utilities/` (`B_matrix`, `D_matrix`, `shape_functions`).
The orchestrator is `linear_gbt_3D.py`.
