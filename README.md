# fem_model

## Testing

Install test dependencies and run the suite (matches CI):

```bash
pip install -r requirements-ci.txt
python -m pytest tests -q
```

CI runs with **`FEM_FORMULATION_CACHE_STRICT_SHAPE=1`** so missing Gauss **`shape_functions`** fail fast. Match that locally:

```bash
export FEM_FORMULATION_CACHE_STRICT_SHAPE=1   # POSIX
python -m pytest tests -q
```

```powershell
$env:FEM_FORMULATION_CACHE_STRICT_SHAPE = "1"
python -m pytest tests -q
```

CI: [.github/workflows/pytest.yml](.github/workflows/pytest.yml) (Python **3.10** and **3.11**).

## Smoke jobs (warping)

Example nonlinear static jobs under [`jobs/`](jobs/):

- [`jobs/job_smoke_nl_eb_warp`](jobs/job_smoke_nl_eb_warp) — TL Euler–Bernoulli with `[warping]`
- [`jobs/job_smoke_nl_ts_warp`](jobs/job_smoke_nl_ts_warp) — TL Timoshenko with `[warping]`
