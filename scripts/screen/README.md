# Repository screen

Scripts that check 1d-fem against its aims: a general 1D FE code with a research-grade
element registry, whose results can be probed at every level of resolution.

- `capability_matrix.py` runs every registered element through every analysis (static,
  nonlinear static, eigen, buckling, harmonic, transient, layup sensitivity) on three
  section kinds (closed form, thin-walled, 2D) and three geometries (straight, rotated,
  curved) on one composite box cantilever. Each case is `ok`, `refused` (a clear
  `NotImplementedError` or `ValueError`), `crash` or `differs` (more than 3 % from the
  linear Timoshenko beam on the same section and geometry). Failures that `process_job`
  logs without raising are read from `logs/traceback.log`. About 1450 cases, roughly an
  hour on four workers.
- `static_probes.py` scans the source: capabilities decided by element names, private
  attribute reads across modules, blade words in the core, the columns of every result
  file, smoothing calls and figures saved outside `figure_style`, and the theory pages
  and tests of each element.

```
python scripts/screen/capability_matrix.py --out screen_out --workers 4
python scripts/screen/static_probes.py --out screen_out
```
