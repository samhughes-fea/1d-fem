# Contributing

The rules are in [AGENTS.md](AGENTS.md) and the shared standard in
[docs/conventions/standard/](docs/conventions/standard/README.md): work through a
pull request with a Conventional Commit title, one concern per pull request, and a
line in [docs/CHANGELOG.md](docs/CHANGELOG.md) under `## [Unreleased]` ending with
the PR number.

## Checks

CI runs the same checks (`pytest.yml`). System package: `libgeos-dev`, for the
optional geometry kernel.

```bash
pip install -r requirements-ci.txt && pip install -e .
pip install -r requirements-sections.txt -r requirements-geometry.txt
FEM_FORMULATION_CACHE_STRICT_SHAPE=1 python -m pytest tests -q -n auto
```

This repository's own conventions are in
[docs/conventions/README.md](docs/conventions/README.md).
