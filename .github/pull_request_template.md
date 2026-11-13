## What this changes

<!-- One or two sentences. Link the issue if there is one. -->

## Why

<!-- What was wrong or missing. For a layout or rendering change, say what the
     figure looked like before and after. -->

## Checks

- [ ] `uv run pytest` passes
- [ ] `uv run ruff check src/ tests/` and `uv run ruff format --check src/ tests/` pass
- [ ] `uv run mypy` passes
- [ ] Tests added for the change (coordinate assertions for layout, golden image for rendering)
- [ ] `CHANGELOG.md` updated under `[Unreleased]`
- [ ] Docs updated if user-visible (`README.md`, `docs/api.md`, `docs/pypedigree_documentation.html`)

## Regenerated golden images

<!-- If tests/visual/golden/ changed, list each image and say why it changed.
     A golden regenerated without explanation can hide a regression.
     Write "none" if you did not touch them. -->

## Test data

- [ ] Any new fixture is synthetic or fully de-identified — no cohort accession
      numbers or other identifiers linkable to a registry
