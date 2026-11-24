# Releasing

The rules are the shared standard's
[releasing](standard/releasing.md) and [versioning](standard/versioning.md): a
release is a pull request titled `chore(release): prepare X.Y.Z`, and the synced
`release.yml` tags and releases its squash commit when it merges. Never push or
move a tag by hand.

This repository's specifics:

- `__version__` lives in `src/onedfem/__init__.py`; `pyproject.toml` reads it
  from there. On `main` it is the next release.
- The changelog is [docs/CHANGELOG.md](../CHANGELOG.md). The release PR renames
  its `## [Unreleased]` section `## [X.Y.Z] - YYYY-MM-DD`.
- 1d-fem pins anisotropic-section-model (the default section provider) and, for
  the geometry extra, the kernel by tag in `requirements-sections.txt` and
  `requirements-geometry.txt` (and `pyproject.toml`): release those first and bump the pins.
- Releases v0.1.0 to v0.4.0 were tagged retrospectively by the old workflow's
  backfill mode; their code reports `__version__` 0.1.0. The shared workflow has
  no backfill mode, and none is needed now that every release is tagged.
- PyPI: `pypi: false` in `.github/repo.yml`. To publish, add a PyPI trusted
  publisher for `onedfem` (workflow `release.yml`, environment `pypi`), create the
  `pypi` environment in the repository settings, and set `pypi: true`.
