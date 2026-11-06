# Linear Benscoter beam (fosdt)

`fosdt` = first-order shear deformation theory. This is a **spanwise 1D** beam with
an independent warping amplitude `psi` from thin-walled section analysis.

Operators live in `utilities/` (`B_matrix`, `D_matrix`). Shape functions and
distributed-load interpolation are inherited from the Timoshenko parent.
The orchestrator is `linear_benscoter_3D.py`.
