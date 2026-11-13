# Releasing pypedigree

Publishing is automated by `.github/workflows/release.yml`, which runs on a
version tag. Authentication uses PyPI **Trusted Publishing** (OIDC), so there
is no API token stored anywhere in the repository or in GitHub secrets.

## One-time setup on PyPI

This must be done once, by hand, before the first release. It cannot be
scripted from the repository.

1. Sign in to <https://pypi.org> and go to **Your projects → Publishing**, or
   for a project that does not exist yet:
   <https://pypi.org/manage/account/publishing/>.
2. Add a **pending publisher** (GitHub) with exactly these values:

   | Field | Value |
   |-------|-------|
   | PyPI project name | `pypedigree` |
   | Owner | `darkxander` |
   | Repository name | `Pypedigree` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

   The environment name must match the `environment: name: pypi` in the
   workflow, or PyPI will reject the token.

3. In the GitHub repository, go to **Settings → Environments → New
   environment** and create one named `pypi`. Adding a required reviewer
   there is recommended: it makes every publish pause for an explicit
   approval.

The project name `pypedigree` was free on PyPI at the time of writing. Verify
it still is before setting up the publisher:

```bash
curl -s -o /dev/null -w '%{http_code}\n' https://pypi.org/pypi/pypedigree/json
```

`404` means the name is available; `200` means it is taken.

## Cutting a release

1. Make sure `main` is up to date and CI is green.

2. Move the `[Unreleased]` section of `CHANGELOG.md` under a new version
   heading with today's date:

   ```markdown
   ## [0.3.1] - 2026-09-02
   ```

3. Bump the version in **one** place — `pyproject.toml`:

   ```toml
   version = "0.3.1"
   ```

   `pypedigree.__version__` reads the installed package metadata, so nothing
   else needs editing. Update `version` and `date-released` in
   `CITATION.cff` to match.

4. Verify locally:

   ```bash
   uv sync --extra dev --extra pandas
   uv run pytest -q
   uv run ruff check src/ tests/ && uv run ruff format --check src/ tests/
   uv build && uvx twine check --strict dist/*
   ```

5. Commit, tag and push. The tag must be the version prefixed with `v` — the
   workflow fails the build if it disagrees with `pyproject.toml`:

   ```bash
   git commit -am "Release 0.3.1"
   git tag -a v0.3.1 -m "Release 0.3.1"
   git push origin main --follow-tags
   ```

6. The workflow then runs the test suite, builds the sdist and wheel, checks
   the metadata with `twine check --strict`, re-runs the tests from inside the
   extracted sdist, publishes to PyPI, and creates the GitHub release with
   the artifacts attached.

## Rehearsing without publishing

Run the workflow manually from the **Actions** tab with `dry_run` left at
`true`. Everything is built and checked, and the publish step is skipped.

## Making the release citable

Once the GitHub release exists, connect the repository to
[Zenodo](https://zenodo.org/account/settings/github/) to mint a DOI for each
release. Add that DOI to `CITATION.cff` as:

```yaml
identifiers:
  - type: doi
    value: 10.5281/zenodo.XXXXXXX
    description: Archived snapshot of this release
```

GitHub reads `CITATION.cff` and shows a **Cite this repository** button on the
repository page.
