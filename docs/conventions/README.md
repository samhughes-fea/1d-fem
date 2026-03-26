# 1D FEM project conventions

Standards for this repository: **1D finite-element beam/frame models in 3D space** (six DOF per node by default), text-driven jobs, and research-oriented validation.

| Document | Scope |
|----------|--------|
| [API_STANDARDS.md](API_STANDARDS.md) | Weak-form assembly (**B**, **D**, **J**), equivalent nodal loads (**N**, **f**), linear `utilities/` layout, nonlinear TL extension, runners and parsers. |
| [COMMIT_STANDARDS.md](COMMIT_STANDARDS.md) | [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/). |
| [TESTING_STANDARDS.md](TESTING_STANDARDS.md) | Pytest layout under [tests/](../../tests), naming, fixtures. |
| [FORMULATION_DOCSTRING_STANDARDS.md](FORMULATION_DOCSTRING_STANDARDS.md) | Checklist for element and utility docstrings (shapes, DOF order, limits). |

**Running tests:** From the repository root, ensure `1d-fem` is on `PYTHONPATH` (or run from root with `python -m pytest`). See root [pytest.ini](../../pytest.ini).

**Jobs:** Input layout and benchmarks are described under [jobs/README_JOBS.md](../../jobs/README_JOBS.md).
