# Contributing to pypedigree

Thanks for considering a contribution. Bug reports with a `.ped` file that
reproduces the problem are especially useful — pedigree topologies are varied
enough that most layout bugs only show up on real data.

## Getting set up

```bash
git clone https://github.com/darkxander/Pypedigree.git
cd Pypedigree
uv sync --extra dev --extra pandas
```

Everything runs through [uv](https://docs.astral.sh/uv/). Python 3.10 or newer.

## Before opening a pull request

```bash
uv run pytest                                   # tests
uv run ruff check src/ tests/                   # lint
uv run ruff format src/ tests/                  # format
uv run mypy                                     # type check
```

CI runs the same four, plus the test suite on Python 3.10–3.14 across Linux
and macOS, and once more with every dependency at its declared minimum
version. All must pass.

## The architecture matters

The package is a four-stage pipeline with no circular dependencies:

```
parser.py -> models.py -> layout.py -> renderer.py / svg_renderer.py
```

Two rules keep it testable, and a change that breaks either will be asked to
change:

- **`layout.py` never imports matplotlib.** It produces pure numeric
  coordinates in `FamilyLayout` dataclasses. This is what makes layout
  testable by asserting on numbers rather than on pixels.
- **`renderer.py` never computes positions.** It only draws from a
  pre-computed `FamilyLayout`.

`CLAUDE.md` documents the full set of design rules, including why phenotype
display decisions are resolved once per figure and why the renderer holds no
module-level mutable state.

## Adding a feature

Work down the pipeline, adding tests at each stage you touch:

1. Data representation → `models.py`
2. File parsing → `parser.py`
3. Positioning → `layout.py`, with tests asserting coordinates
4. Drawing → `renderer.py` and `svg_renderer.py`. A decision about what to
   draw (a color, a text, a legend row, a symbol's shape) goes in
   `phenotypes.py`, `carriers.py` or `symbols.py`, which both renderers read, so the PNG and the interactive
   SVG stay identical
5. New style parameters → `config.py` (all visual parameters live there; no
   magic numbers in layout or renderer)
6. New CLI flags → `cli.py`
7. A user-visible feature → add an example `.ped` under `tests/data/`

Then update the docs: `README.md`, `CHANGELOG.md`, `docs/api.md` and
`docs/pypedigree_documentation.html`. The PDF is regenerated from the HTML.

## Tests

- Unit tests mirror the `src/` layout: `tests/test_<module>.py`
- Fixtures live in `tests/data/`, grouped by what they exercise
- Layout changes want coordinate assertions, not just smoke tests
- Visual changes want a golden image:

  ```bash
  uv run pytest tests/visual/ --regen-golden
  ```

  Inspect the regenerated images before committing them, and say in the pull
  request why each one changed. A golden image regenerated without explanation
  hides a regression. The comparison counts changed pixels and fails above
  0.01% of them: when a test fails, `tests/visual/golden/<name>.diff.png`
  shows the changes in red. The golden images hold for matplotlib 3.7 to
  3.10 (`GOLDEN_MATPLOTLIB`); with another version the comparison is skipped,
  so regenerate them with the version pinned in `uv.lock` (`uv run`).

**Do not commit test data derived from real individuals.** Every fixture in
`tests/data/` is synthetic — see `tests/data/README.md`. Accession numbers and
sample barcodes are pseudonyms, not anonymisation, and a large pedigree
carrying a rare phenotype can be identifying on its own. `tests/test_fixtures.py`
enforces this: it fails on any accession-shaped identifier or populated
date-of-birth column.

Generated fixtures come from `tests/data/generate_fixtures.py` with a fixed
seed. If you edit one by hand the reproducibility test fails — change the
generator instead and re-run it.

## Reporting a bug

Include:

- A minimal `.ped` file that reproduces it (synthetic — see above)
- The exact command or code you ran
- What you expected and what you got; attach the rendered image if the problem
  is visual
- `pypedigree --version`, your Python version and your OS

## Style

- Google-style docstrings on public functions, with `Args:` and `Returns:`
- Type hints on everything; `mypy --strict` must pass
- Comments explain *why*, not *what*
- Line length 100, enforced by ruff

## Releases

Maintainers only: see [docs/RELEASING.md](docs/RELEASING.md).
