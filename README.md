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

After features land on **`main`**, update other clones with:

```bash
git checkout main && git pull origin main
```

Re-install **`requirements-ci.txt`** if it changed (new imports for CI parity).

## Smoke jobs (warping)

Example nonlinear static jobs under [`jobs/`](jobs/):

- [`jobs/job_smoke_nl_eb_warp`](jobs/job_smoke_nl_eb_warp) — TL Euler–Bernoulli with `[warping]`
- [`jobs/job_smoke_nl_ts_warp`](jobs/job_smoke_nl_ts_warp) — TL Timoshenko with `[warping]`
