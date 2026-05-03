# Commit standards

This repository follows **[Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/)** for commit messages.

## Format

```text
<type>[optional scope]: <description>

[optional body]

[optional footer(s)]
```

- **type** (required): e.g. `feat`, `fix`, `docs`, `test`, `refactor`, `chore`, `ci`, `build`, `style`, `perf`.
- **scope** (optional): area of the codebase, e.g. `elements`, `processing`, `simulation`, `post_processing`, `workflow`, `conventions`.
- **description** (required): short summary after the colon and space.
- **body** (optional): extra context, one blank line after the description.
- **footer** (optional): e.g. `BREAKING CHANGE:`, `Refs: #123`.

## Examples

- `docs(conventions): tailor API standards for 1D FEM weak form`
- `feat(elements): add Timoshenko consistent mass for modal analysis`
- `fix(processing): correct condensation index map for prescribed DOF`
- `test(elements): parametrize geometric stiffness buckling check`

For full rules, see **[Conventional Commits v1.0.0](https://www.conventionalcommits.org/en/v1.0.0/)`.

**Note:** Present proposed commit messages for review before pushing shared branches when working in a team.

## Pull request size (optional)

For easier review, split unrelated work into separate PRs by theme (e.g. `docs:` only, `chore(jobs):` for `jobs/**` `simulation_settings.txt` churn, benchmarks isolated from library changes). Keep each PR focused and conventionally titled.

## Releases and changelog

When tagging a release, move **`## [Unreleased]`** in **[`CHANGELOG.md`](../CHANGELOG.md)** into a dated section **`## [x.y.z] - YYYY-MM-DD`** per [Keep a Changelog](https://keepachangelog.com/).
