# CLAUDE.md - PyPedigree Development Guide

## Project Overview

pypedigree is a Python package that generates publication-quality pedigree (family tree) diagrams
from standard pedigree files used in genetics research.

Architecture follows a strict four-stage pipeline:

```
Parser -> Data Model -> Layout Engine -> Renderer
```

Each stage is independently testable with no circular dependencies.

## Quick Reference

```bash
# Install (development mode with all dev dependencies)
uv sync --extra dev

# Install with pandas support
uv sync --extra dev --extra pandas

# Run the CLI
uv run pypedigree tests/data/basic/simple_trio.ped -o output.png

# Run as module
uv run python -m pypedigree tests/data/basic/simple_trio.ped --show

# Run all tests
uv run pytest

# Run tests with coverage
uv run pytest --cov=pypedigree --cov-report=term-missing

# Lint
uv run ruff check src/ tests/

# Format
uv run ruff format src/ tests/
```

## Architecture

```
src/pypedigree/
  errors.py    -> Custom exceptions + warning categories (no dependencies)
  models.py    -> Data structures (depends on: errors)
  config.py    -> PedigreeStyle dataclass + style resolution (depends on: errors)
  parser.py    -> File I/O + validation + pandas (depends on: models, errors; optional: pandas)
  layout.py    -> Coordinate math, incl. circular layout (depends on: models, config, symbols, errors)
  carriers.py  -> Carrier notation decisions shared by both renderers, no matplotlib (depends on: models, errors)
  symbols.py   -> Standard notation (pregnancies, adoption, childlessness, donors, gender, consultand):
                  decisions and geometry shared by layout and both renderers, no matplotlib (depends on: models)
  phenotypes.py -> Phenotype display decisions shared by both renderers: column modes, band colors/text,
                   legend rows and geometry; uses only matplotlib.colors (depends on: models, config, carriers, errors)
  renderer.py      -> Matplotlib drawing (depends on: models, config, layout, carriers, phenotypes, symbols)
  svg_renderer.py  -> Interactive SVG rendering (depends on: models, config, layout, carriers, phenotypes, symbols)
  api.py           -> PedigreeBuilder fluent API (depends on: models, config, renderer, svg_renderer)
  cli.py       -> Click CLI (depends on: parser, renderer, config)
  __init__.py  -> Public re-exports (includes compute_layout_circular)
  __main__.py  -> python -m entry point
```

## Key Design Rules

1. **layout.py never imports matplotlib.** It produces pure numeric
   coordinates in FamilyLayout dataclasses.

2. **Hierarchical edges are orthogonal.** Every edge `compute_layout` emits is
   horizontal or vertical, as pedigree nomenclature requires; an offset single
   child hangs from a short sibship-level segment. `TestOrthogonalDescent`
   checks this on every fixture and on the large synthetic pedigrees. (The
   circular layout uses straight and curved lines by design.)

3. **Sibship lines never overlap.** `_HierarchicalPlacer` packs each union's
   descendants as a block against its neighbours' contours, so two couples'
   children never interleave. Loops are drawn side by side unless that
   stretches a sibship line over another one; then `_place_hierarchically`
   duplicates the partner and lays out again. To make that rare,
   `_children_group` orders each sibship so a loop's branches meet: a child
   whose line leads to a partner in a family placed later goes last, one whose
   line leads to a partner already drawn goes first. `loops_to_break` breaks
   one loop per tangle (group of overlapping lines) per pass, the one
   stretching the widest line, and every tangled loop only after half of
   `_MAX_LOOP_BREAKING_PASSES`. Don't add passes that move nodes individually
   after placement (the old sweep / barycenter / refinement passes broke
   this). `TestSibshipSeparation` checks it; `TestLoopDuplicates` guards the
   duplicate count. A partner normally takes the side given by the
   partnership's order, but the end twins of a twin group take theirs on the
   group's outer side (`_twin_outer_sides`), so no spouse sits between twins
   under the twin connector; loop partners keep their forced right side.
   `TestTwinSpouses` checks it.

4. **renderer.py never computes positions.** It only draws from
   pre-computed FamilyLayout objects. It supports waypoint-based
   polyline edges for curved/non-orthogonal connections (circular layout).

5. **models.py has no I/O.** It is pure data structures.

6. **All visual parameters live in PedigreeStyle.** No magic numbers
   anywhere in layout or renderer code.

7. **Type hints on all public functions.** Google-style docstrings.

8. **Style resolution is centralized in `config.py`.** `resolve_style()` fills
   in pedigree-dependent style defaults (palette colors) and is pure and
   idempotent — it never mutates its argument. Every entry point (CLI,
   `PedigreeBuilder`, `render_pedigree`, `render_interactive_svg`) calls it, so
   the API and the CLI always produce the same figure. `apply_pgs_spacing()` is
   *not* idempotent: each renderer's entry point (`render_pedigree()`,
   `render_interactive_svg()`) applies it exactly once, and nothing else does.

9. **Phenotype display decisions are resolved once per figure, in
   `phenotypes.py`, and shared by both renderers.** `PhenotypeRender` maps each
   column to its effective mode (`binary`/`gradient`/`categorical`/`outline`/
   `carrier`/`pgs`) and, for categorical columns, to its value -> color mapping.
   `band_cell()` gives each band's color and text, `legend_rows()` the legend's
   rows. The matplotlib renderer and the SVG renderer only translate these into
   artists or elements, so the PNG, its legend and the interactive SVG can never
   disagree. Put any new display decision there, not in a renderer. Nothing is
   cached in module-level state between renders.

10. **No module-level mutable state in the renderer.** Anything derived from a
   pedigree or a style is computed per call and passed down. Caching it at
   module level leaks between renders in one process (Jupyter, batch loops,
   web services) and makes output depend on render order.

11. **Diagnostics use the `Pedigree*Warning` categories** from `errors.py`, never
   a bare `UserWarning`, so callers can filter or escalate them. The CLI
   captures them with `warnings.catch_warnings` and re-emits them as
   `Warning: <message>` on stderr.

12. **Standard notation is decided in `symbols.py`.** `models.Notation` holds an
    individual's pregnancy, adoption, childlessness, assisted-reproduction role,
    gender identity and consultand flag. `symbols.py` decides what they look like
    (shape sex from gender, AMAB/AFAB and "SB"/"ECT" notes, triangle and bracket
    geometry, the letter inside) with no matplotlib, and functions whose result
    is not symmetric top to bottom take an `up` sign so the SVG (y down) uses
    them too. The layout copies the decisions into `NodeLayout` (`sex`,
    `is_deceased`, `notes`, `marker_depth`) and emits the childlessness marker as
    `EdgeType.CHILDLESS` edges and an adopted-in child's line as dashed, so both
    renderers draw them without knowing about notation.

13. **The legend panel derives its own geometry, from the font.** `legend_metrics()`
    (in `phenotypes.py`) returns a `LegendMetrics` from the row count and
    `style.font_size`; one legend unit is a fixed multiple of the font size
    (`LEG_UNIT_X_EM`, `LEG_UNIT_Y_EM`), so text and swatches scale together and
    never with the figure. Both renderers use it: `_add_legend_axes()` grows the
    matplotlib figure downward by the panel height rather than taking a share of
    the pedigree's, and `_svg_legend()` places the panel below the pedigree.
    Each row type has its own drawing function. Legend geometry constants
    (`LEG_*`) are not user-tunable style. `TestLegendPhysicalSize` checks it.

## File Format

See `docs/pedigree_standard.md` for the complete specification.

Standard 6-column format (TSV/CSV/space-separated):

| Col | Name | Values |
|-----|------|--------|
| 1 | family_id | Family identifier |
| 2 | sample_id | Individual identifier |
| 3 | paternal_id | Father's sample_id (0/NaN = unknown) |
| 4 | maternal_id | Mother's sample_id (0/NaN = unknown) |
| 5 | sex | 1=male, 2=female, 0/-9=unknown |
| 6 | phenotype | 1=unaffected, 2=affected, 0/-9=unknown |

Extended format adds:
- Column 7: date_of_birth (YYYYMMDD, YYYYMM, or YYYY)
- Column 8: twins (comma-separated TYPE:sample_id entries)
- Column 9: deceased (1=deceased, 0=alive)
- Column 10: proband (1=proband, 0=not proband)
- Columns with `a_` prefix: annotation columns (text displayed below individual's ID label, not phenotypes). Missing values (0, -9, NA, empty) produce no annotation line. Generation gap auto-increases when present.
- Standard notation columns, by header name (`_NOTATION_ROLES` in parser.py; codes
  case-insensitive, unknown code -> `PedigreeParseError` with the line):
  `pregnancy` (P/SAB/TOP/ECT/SB), `adopted` (in/out), `childless` (choice/infertility),
  `art` (D/S), `gender_identity` (man/woman/nonbinary; `gender` means this only when
  a `sex` column is also present, otherwise it is the sex column), `consultand` (1).
  `build_notation()` accepts the same codes (used by `PedigreeBuilder`).
- Columns 11+: extra phenotypes (requires header row with column names)

Smart header-based column mapping:
- When a header is present, columns are matched by name (not position)
- Any column after the core 5 (FID, ID, Father, Mother, Sex) that doesn't match
  a known structural name is treated as a phenotype column
- Supports aliases: Pat/Father/Dad, Pheno/Phenotype/Affected, etc.
- `ColumnMapping` dataclass in parser.py handles the name→role mapping
- `_parse_individual_mapped()` uses ColumnMapping instead of fixed indices.
  `_parse_text` uses it whenever the header has extra phenotype, annotation or
  notation columns, or a structural column outside the standard order
  (`ColumnMapping.is_positional`); only a header in the standard order goes
  through the positional `_parse_individual()`

Multi-phenotype support:
- Inline phenotypes parsed from columns 11+ when header is present
- Separate .pheno TSV file via `--pheno` CLI option
- Binary phenotypes: 1=unaffected (white), 2=affected (color), 0=unknown (grey), -9=no data (white with "NA" text)
- When all bands are NA for an individual, a single centered "NA" is shown instead of repeating per band
- Continuous phenotypes → color gradient bands
- PGS columns: any column starting with `PGS` (case-insensitive) → gaussian distribution curve. Multiple PGS columns (e.g. `PGS_autism`, `PGS_IQ`) show colored markers on a single curve. Marker colors via `pgs_colors` style option. `couple_gap` and `sibling_gap` auto-increase when PGS present.
- Colors from `phenotype_colors` in style JSON, or auto-assigned from palette
- Comma decimal separator: parser accepts `0,712` as `0.712` (European locale)
- Gradient band values always displayed regardless of intensity
- AUTO mode (no prefix) is resolved **once per column** from all of its values:
  any string value -> categorical, any numeric value outside 0/1/2 -> gradient,
  otherwise binary. Never per individual, so a column cannot mix band types.

Additional PedigreeStyle fields:
- `pgs_colors: dict[str, str]` — PGS column name → hex color for marker line (default: `#e74c3c`)
- `gradient_range: dict[str, list[float]]` — custom min/max for gradient normalization (e.g. `{"SRS": [40, 90]}`)
- `gradient_diverging: dict[str, dict]` — tricolor diverging gradient config (e.g. `{"IQ": {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}}`)
- `zero_pad_ids: int` — zero-pad IDs to this width (0 = no padding)
Warnings and errors:
- `PedigreeWarning` (base), `PedigreeStyleWarning`, `PedigreeDataWarning`
- Unknown keys in a style JSON warn with a suggested spelling (they are still ignored)
- Invalid colors raise `PedigreeError` naming the option, checked before any drawing

- Prefix-based display modes in column names:
  - `b_` → binary (color/white/grey/NA, no text by default)
  - `g_` → gradient (white→color continuous scale)
  - `gi_` → gradient inverted (color→white scale; high value = white, low = full color; for metrics where low = bad)
  - `o_` → outline (colored shape border, not a band)
  - `c_` → categorical (distinct color per unique string value)
  - `cr_` → carrier status of one gene or condition (0/1/2), drawn as a patterned section, not a band
  - No prefix → auto-detect from data

## CLI Options

```bash
pypedigree INPUT.ped [OPTIONS]

-o, --output PATH        Output file path
--format [png|svg|pdf]   Output format (auto-detected from extension)
--family ID              Render only this family ID
--dpi INT                Output resolution (default: 150)
--show                   Display interactively instead of saving
--style PATH             JSON file with style overrides
--separator [tab|comma|space]  Force field separator
--pheno PATH             Separate .pheno TSV file with extra phenotype data
--annotations PATH       Separate TSV file with annotation text (displayed below IDs)
--outline-color COLOR    Shape/deceased outline color (e.g. 'grey', '#888888')
--carrier-notation [nsgc2022|legacy]  Carrier notation (default: nsgc2022)
--no-labels              Hide sample ID labels
--no-legend              Hide the phenotype legend
--show-phenotype-labels  Show YES/NO text inside binary phenotype bands
--interactive            Generate interactive SVG with tooltips and click highlighting
--circular               Use circular/radial layout (founders at center, descendants on concentric rings)
--plink                  Parse as PLINK format (ignore genotype columns after col 6)
```

Circular layout: `compute_layout_circular` reuses `_place_hierarchically` and
wraps its x order around the circle (`_circular_positions`), so it inherits the
hierarchical guarantees (no interleaved sibships, loops broken with duplicates).
Shapes stay upright; `NodeLayout.label_angle` sends the ID label radially
outward (`_draw_radial_label`). Twin connectors hang from the sibship arc, and
duplicate arcs are emitted as in the hierarchical layout. The renderer still
supports `NodeLayout.rotation` for API callers: rotated labels sit below the
rotated outline, band text is never upside down, and the proband marker stays
upright.

Interactive SVG: `svg_renderer.py` draws the same figure as `renderer.py`
(bands, outline colors, PGS curves, annotations, carrier notation, legend) from
the same `phenotypes.py` / `carriers.py` / `symbols.py` decisions, plus per-individual tooltips,
click highlighting and zoom/pan. Families are laid side by side; element ids are
generated (`f<family>n<node>`), never derived from sample IDs. Text is placed
with explicit baseline offsets (`_BASELINE_SHIFT`), not `dominant-baseline`,
which librsvg and older Inkscape ignore. Font sizes and line widths in points are
scaled by `_PT` layout units per point. No circular layout in the SVG.

Validation warnings are printed to stderr as `Warning: <message>`.

## Testing

- All tests in `tests/` directory, mirroring `src/` structure
- Test data files in `tests/data/` with subdirectories: `basic/`, `consanguinity/`,
  `extended/`, `twins/`, `multi_partner/`, `phenotype/`, `plink/`, `format/`, `complex/`,
  `nomenclature/`
- **All fixtures are synthetic.** Never add pedigrees derived from real
  individuals, even pseudonymised — accession numbers and sample barcodes are
  pseudonyms, not anonymisation. See `tests/data/README.md`.
  - `tests/data/generate_fixtures.py` writes the generated ones from a fixed
    seed: `uv run python tests/data/generate_fixtures.py`. Among them,
    `complex/multigenerational.ped` (81 people, five generations, two loops) is
    the gallery's large example and a visual-regression golden
  - `tests/test_fixtures.py` asserts the files match the generator output and
    that no fixture carries an accession-shaped identifier or a populated DOB
- `conftest.py` sets matplotlib to Agg backend (no GUI windows)
- `tests/test_layout_large.py` — stress tests for large pedigrees (100-200 individuals)
- `conftest.py` generators: `generate_clinical_pedigree` (descent tree with outside
  partners, optional cousin marriages), `generate_realistic_pedigree` (every child
  marries into another branch: worst case for loops), deep and wide pedigrees
- `tests/visual/test_visual_regression.py` — comparison against golden images by
  changed-pixel count (a channel moving by > 0.1; fails above 0.01% of pixels or on
  any size change). Renders are pixel-reproducible across platforms and matplotlib
  3.7-3.10, so any failure is a real change; a mean over the image hid regressions.
  matplotlib 3.11 places text about a pixel lower, so outside `GOLDEN_MATPLOTLIB`
  (3.7 to <3.11) the comparison is skipped, not failed
  - Golden images in `tests/visual/golden/`; a failure writes `<name>.actual.png`
    and `<name>.diff.png` (changes in red) next to the golden
  - Regenerate: `uv run pytest tests/visual/ --regen-golden`
- `tests/test_api.py` — PedigreeBuilder programmatic API tests
- `tests/test_carrier.py` — Carrier notation: NSGC 2022 sections and `cr_` columns, legacy half-fill / dot (model, parsing, rendering, API, CLI)
- `tests/test_svg_renderer.py` — Interactive SVG renderer (XML validity, nodes, tooltips, interactivity, CLI;
  `TestSvgPhenotypes` checks each band's color against `band_cell`, legend, PGS, outline, annotations;
  `TestSvgFixtures` renders every fixture and checks ids are unique; multi-family placement)
- `tests/test_phenotype_resolution.py` — Style resolution (purity, idempotence),
  per-column AUTO mode detection, categorical color assignment, API/CLI parity,
  color validation
- `tests/test_layout_circular.py` — Circular layout: ring geometry, spacing,
  rotation, edge composition, arc waypoints, rendering across 15 fixtures
- `tests/test_legend.py` — Legend panel metrics, physical size and row composition
- `tests/test_fixtures.py` — Fixture provenance: reproducibility and identifier scan
- `tests/test_nomenclature.py` — Standard notation: parsing and codes, `symbols.py`
  decisions and geometry, layout (shape from gender, dashed adopted-in line,
  childless edges, circular direction), PNG and SVG drawing, builder
- CI also runs the suite with the newest release of every dependency
  (`uv sync --upgrade`, job `test-latest-versions`), since `uv.lock` pins older ones
- CI also runs the suite with every dependency at its declared minimum version
  (`uv pip install --resolution lowest-direct`) on Python 3.10; raising a lower
  bound in `pyproject.toml` is a user-visible change
- Run: `uv run pytest -v`

## Visual Reference

Reference images for every visual element are in `docs/imgs/`.
The `diagram_informations.png` shows the labeled structure of a
complete pedigree diagram (relationship line, descent line,
sibship line, individual lines).

## Documentation

- `docs/pypedigree_documentation.html` — Full HTML documentation
- `docs/pypedigree_documentation.pdf` — PDF version of the documentation
- `docs/pedigree_standard.md` — File format specification
- `docs/api.md` — API reference for all public functions, classes, and types
- `docs/examples_gallery.md` — Gallery of example files with CLI commands
- `docs/RELEASING.md` — Release procedure and the one-time PyPI trusted-publishing setup
- `docs/imgs/` — Reference images for visual elements
- `docs/imgs/examples/` — Generated example images for documentation

## PLINK Integration

Native support for PLINK bioinformatics formats:

- `parse_fam_file(path)` — reads `.fam` files (6-column, no header)
- `parse_plink_ped(path)` — reads PLINK `.ped` files (ignores genotype columns after col 6)
- `export_fam(pedigree, path)` — writes `.fam` format (excludes phantoms/duplicates)
- CLI: `.fam` extension auto-detected; `--plink` flag for PLINK `.ped` files

## Pandas Integration

Optional pandas support via `parse_pedigree_from_dataframe()`:

```python
import pandas as pd
from pypedigree import parse_pedigree_from_dataframe

df = pd.read_csv("family.ped", sep="\t", header=None)
ped = parse_pedigree_from_dataframe(df)
```

Install with: `pip install pypedigree[pandas]`

## Adding a New Feature

1. If it affects data representation, update `models.py` first
2. If it affects file parsing, update `parser.py` with tests
3. If it affects positioning, update `layout.py` with coordinate tests
4. If it affects drawing, put the decision in `phenotypes.py` / `carriers.py` /
   `symbols.py` and draw it in both `renderer.py` and `svg_renderer.py`, with tests
5. Update `config.py` if new style parameters are needed
6. Update CLI if new flags are needed
7. Add example `.ped` file if the feature is user-visible
