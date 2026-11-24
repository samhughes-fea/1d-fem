<!-- guidelines:start (shared development standard 0.1.0; edit it there, not here) -->
## Shared rules

These rules hold in every tool repository; the full standard is in
[docs/conventions/standard/](docs/conventions/standard/README.md). The rules
under "Hard rules" below this block are this repository's own.

1. **Read before changing.** This file, `SCOPE.md`, then the plan doc named in
   `.github/repo.yml`.
2. **Work through pull requests.** `main` is protected. One concern per pull
   request; the title is a Conventional Commit
   ([commits](docs/conventions/standard/commits.md)) and becomes the squash commit.
3. **Keep the plan doc current** in the same pull request as the work.
4. **No private data.** Never commit confidential inputs or name a private
   project ([data and security](docs/conventions/standard/data-and-security.md)).
5. **Generic source.** `src/` and `tests/` never name or special-case a
   downstream project.
6. **Deprecate before removing** ([API stability](docs/conventions/standard/api-stability.md)).
7. **Golden numbers move only on purpose**, with the `golden-change` label and
   the reason in the changelog.
8. **Never remove `needs-decision`**; record a lasting decision in `docs/decisions/`.
9. **Releases follow the workflow.** Prepare the release pull request; never push
   or move a tag ([releasing](docs/conventions/standard/releasing.md)).
10. **Label everything**: one type, its areas, any flag
    ([labels](docs/conventions/standard/labels.md)).
11. **Example figures are committed** and embedded in the example README; commit a
    regenerated PNG only when the example's subject changed.
12. **Re-run only what failed** while fixing a failing check.
<!-- guidelines:end -->

Read this before changing the repository. It applies to people and to coding
agents alike. [SCOPE.md](SCOPE.md) says what belongs here;
[docs/conventions/README.md](docs/conventions/README.md) lists this repository's
own conventions; [CONTRIBUTING.md](CONTRIBUTING.md) lists the checks.

## Hard rules

These are this repository's own rules, numbered on from the shared ones above.

13. **Along the span only.** Section properties come from the section provider
    (`onedfem.sections`, anisotropic-section-model by default) and laminates from
    laminate-constitutive-model; 1d-fem does not compute a cross-section's
    stiffness itself.
14. **Imports point down.** 1d-fem may import ASM, LCM and (optionally) the
    geometry kernel; none of them imports 1d-fem.
15. **The section contract** (`SECTION_API`) is agreed with ASM; a change to it is
    a new contract version in both repositories.
16. **Axis and strain order.** 1d-fem keeps its own axis and strain order
    internally ([NODAL_DOF_LAYOUT](docs/conventions/NODAL_DOF_LAYOUT.md),
    [SECTION_STIFFNESS_ORDER](docs/conventions/SECTION_STIFFNESS_ORDER.md)) and
    converts at its boundaries, for example to BeamDyn order.
17. **Weak-form assembly.** Element matrices follow
    [API_STANDARDS](docs/conventions/API_STANDARDS.md) and their docstrings
    [FORMULATION_DOCSTRING_STANDARDS](docs/conventions/FORMULATION_DOCSTRING_STANDARDS.md).

## Working conventions

- Commit scopes and `area:` labels are in [.github/repo.yml](.github/repo.yml) and
  [.github/labels.yml](.github/labels.yml).
- Run the checks in [CONTRIBUTING.md](CONTRIBUTING.md) before pushing.
- Tests follow [TESTING_STANDARDS](docs/conventions/TESTING_STANDARDS.md) on top
  of the shared [testing](docs/conventions/standard/testing.md) standard.
- The changelog is [docs/CHANGELOG.md](docs/CHANGELOG.md).
