# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0] - 2026-11-09

Rendering-correctness release. Figures produced by this version differ from
0.3.0 wherever a pedigree carries phenotypes — the changes below are fixes,
but they change output, so regenerate any figure you keep under version
control.

### Breaking

- Carrier status follows the NSGC 2022 focused revision of the pedigree
  nomenclature (Bennett et al., J Genet Couns 2022;31(6):1238-48) by
  default. The symbol is divided into vertical sections, an affected one
  first, then one per carrier result, each with its own fill pattern defined
  in the legend; obligate carriers get the same pattern at half the density.
  Figures with a `carrier` column change. The previous half-filled shape and
  center dot remain available with `carrier_notation="legacy"` (style option,
  style file key, or `--carrier-notation legacy`).
- `render_pedigree()` now draws phenotype bands and a legend on a
  multi-phenotype pedigree where 0.3.0 drew plain black-and-white shapes.
  Nothing to change in calling code; the output is simply what the CLI
  already produced.
- An invalid color in a style raises `PedigreeError` instead of surfacing
  matplotlib's `ValueError`. Code with `except ValueError` around a render
  must catch `PedigreeError` (or `PedigreeError, ValueError`) instead.
- `PedigreeStyle.from_dict()` emits `PedigreeStyleWarning` for unknown keys.
  Under `-W error` this now raises; pass `warn_unknown=False` to opt out.
- The circular layout types a single-child union's edge as
  `EdgeType.DESCENT` rather than `EdgeType.INDIVIDUAL`. Rendering is
  unaffected; only callers inspecting `EdgeType` see the difference.
- Figures change even with unchanged input: legend rows for unprefixed
  continuous columns became gradient bars, categorical colors are assigned
  differently, and CLI figures with a column named exactly `PGS` are no
  longer spaced twice as wide as intended.
- `pypedigree.__version__` is read from installed package metadata. Running
  from an uninstalled source tree now reports `0.0.0.dev0`.
- A column named `gender` next to a `sex` column is read as the gender
  identity (see Added) instead of an extra phenotype column. Alone, `gender`
  is still the sex column.

### Added

- Standard notation beyond sex, affection status and death, following the
  NSGC nomenclature (Bennett et al. 2008 and its 2022 revision), in both the
  PNG and the interactive SVG, hierarchical and circular layouts. One header
  column each, case-insensitive codes:
  - `pregnancy`: `P` ongoing pregnancy ("P" inside the symbol), `SAB`
    miscarriage, `TOP` termination (struck through) and `ECT` ectopic
    pregnancy ("ECT" below), drawn as triangles, and `SB` stillbirth (struck
    through, "SB" below).
  - `adopted`: `in` or `out`. Brackets around the symbol; the line to
    adoptive parents is dashed.
  - `childless`: `choice` or `infertility`. A line below the symbol ending in
    one bar or two.
  - `art`: `D` gamete donor or `S` surrogate, written inside the symbol.
  - `gender_identity`: `man`, `woman` or `nonbinary`. The symbol follows the
    gender identity and the `sex` column becomes the sex assigned at birth,
    written below as AMAB or AFAB when the two differ.
  - `consultand`: `1` points an arrow without "P" at the symbol.

  The same values are `PedigreeBuilder.add_individual()` keywords
  (`pregnancy=`, `adopted=`, `childless=`, `art=`, `gender=`,
  `consultand=`), and `build_notation()` turns codes into the new
  `Notation` attached to every `Individual` and `NodeLayout`. New enums:
  `Pregnancy`, `Adoption`, `Childlessness`, `ReproductiveRole`, `Gender`;
  new edge type `EdgeType.CHILDLESS`. An unknown code is a
  `PedigreeParseError` naming the line. The decisions both renderers share
  live in the new `pypedigree.symbols` module. Pedigrees without these
  columns render exactly as before. Example fixture:
  `tests/data/nomenclature/nomenclature.ped`, also a visual-regression
  golden.
- The interactive SVG (`--interactive`, `render_interactive_svg()`) draws
  the same figure as the PNG: phenotype bands (binary, gradient, diverging,
  categorical), `o_` outline colors, PGS distribution curves, `a_`
  annotations and the legend, which `--no-legend` (or `show_legend=False`)
  hides. It used to draw shapes and clinical notation only, with phenotype
  values in the tooltips. Tooltips now also list categorical values, and
  "NA" rather than "nan" for missing data.
- `pypedigree.phenotypes`: the phenotype decisions both renderers share
  (effective display mode of each column, color and text of every band,
  rows and geometry of the legend). The matplotlib renderer's output is
  unchanged, pixel for pixel; the SVG can no longer diverge from it.
- `cr_` column prefix: carrier status of one gene or condition per column
  (`0` not a carrier, `1` carrier, `2` obligate carrier, `-9`/`NA` unknown).
  Each `cr_` column is its own patterned section and legend row, so a
  carrier of two conditions shows both, as NSGC 2022 recommends. With
  phenotype bands, each carrier source is an extra band below them.
  Example fixture: `tests/data/phenotype/carrier_genes.ped`.
- `PedigreeStyle.carrier_notation` (`"nsgc2022"` or `"legacy"`) and the
  matching `--carrier-notation` CLI option. Any other value raises
  `PedigreeError`.
- `resolve_style()`, `apply_pgs_spacing()` and `is_pgs_column()` are exported
  from the package for callers that need to know how a pedigree will be styled.
- `tests/test_phenotype_resolution.py`: regression tests for style resolution,
  per-column mode detection and categorical color assignment.
- `tests/test_cli.py::TestCLIWarnings`: coverage for CLI warning formatting.

- Warning categories so callers can filter pypedigree diagnostics
  independently of other libraries: `PedigreeWarning` (base),
  `PedigreeStyleWarning` and `PedigreeDataWarning`. For example
  `warnings.simplefilter("error", PedigreeDataWarning)` turns a
  sex-inconsistent parent into a hard failure.
- `PedigreeStyle.from_dict()` now reports unknown style options as a
  `PedigreeStyleWarning`, with a suggested spelling when one is close
  (`phenotype_color` -> "did you mean 'phenotype_colors'?"). Silently
  ignoring a typo was indistinguishable from the option having no effect.
  Pass `warn_unknown=False` to accept a dictionary carrying extra data.
- Invalid colors now raise `PedigreeError` naming the style option they came
  from, instead of surfacing a raw matplotlib `ValueError: Invalid RGBA
  argument`. Every color in the style is checked before any drawing starts,
  so all offending values are reported at once.
- `tests/test_layout_circular.py`: 169 tests for the circular layout, which
  previously had none — ring geometry, node spacing, rotation, edge
  composition, arc waypoints and rendering across 15 fixtures.
- `tests/test_legend.py`: legend panel geometry and row composition.

- `CITATION.cff` so the project is citable, with GitHub's "Cite this
  repository" button. Authors: Alexandre Mathieu, Freddy Cliquet.
- `.github/workflows/release.yml`: tag-driven PyPI publishing via trusted
  publishing (OIDC, no stored API token). The workflow refuses to build when
  the tag disagrees with the version in `pyproject.toml`, runs the test suite,
  checks metadata with `twine check --strict`, re-runs the tests from inside
  the extracted sdist, then publishes and creates the GitHub release.
- `docs/RELEASING.md`: release procedure and the one-time PyPI
  trusted-publishing setup (a manual step on pypi.org).
- `Alexandre Mathieu` added to the package authors, matching the commit history.

- `mypy --strict` is enforced in CI. The package ships `py.typed`, so the
  annotations are a public contract; they were never verified before. Fixing
  the 27 findings turned up one real annotation bug: `_detect_separator()` was
  annotated `-> str` but returns `None` to signal whitespace splitting.
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), issue
  forms for bugs and feature requests, and a pull-request template.
- `tests/test_svg_renderer.py::TestSvgPhenotypes`: checks that every SVG band
  has the color the matplotlib renderer gives it, and that the legend, PGS
  curves, outline colors and annotations are drawn.

### Changed

- Fewer duplicated individuals in pedigrees with loops. Siblings are now
  ordered so the two branches of a loop meet, at every generation from the
  common ancestors down (previously only the partner themselves moved to
  the end of their sibship), and when lines still overlap, one loop per
  group of overlapping lines is broken before laying out again (previously
  every loop touching an overlapping line, at once). Over 24 synthetic
  clinical pedigrees with one cousin marriage, 1 needs a duplicate instead
  of 13; with three cousin marriages each, 35 duplicates instead of 57;
  12 % fewer over pedigrees with random mating. Layout takes about 6 ms for
  200 individuals with six loops. The reordering can put a sibling out of
  birth order. `TestLoopDuplicates` guards these counts.
- Circular layout rebuilt on the hierarchical one: each individual's
  horizontal position becomes an angle (counter-clockwise, centred at the
  bottom) and each generation a ring. It inherits the hierarchical
  guarantees whatever the number of founder couples: the sector model it
  replaces put duplicates among the founders, overlapping other nodes at the
  centre, and drew `consanguinity_first_cousins.ped` as a column instead of
  a circle. A small family now occupies only the arc it needs.
- Circular layout: shapes stay upright. A square turned to face outward read
  as the unknown-sex diamond, and bands changed direction around the circle.
  ID labels are set radially outside the shape, away from incoming lines
  (`NodeLayout.label_angle`).
- Hierarchical layout: horizontal placement rewritten. Each union's
  descendants now form a block; sibling blocks are packed against each
  other generation by generation, and every couple is centred over its
  children. The sweep, crossing-reduction and refinement passes, and the
  15-individual threshold that switched them on, are gone: small and large
  pedigrees follow the same rules. Figures change for 17 of the 40 test
  fixtures, and large tree pedigrees come out about a third narrower.
- Hierarchical layout: a loop (both partners descend from people in the
  pedigree, as in consanguinity) is drawn with the partners side by side,
  the second partner's sibship line stretching to reach them. Where that
  line would overlap another sibship, the partner is drawn twice instead —
  among their siblings and, as a duplicate joined by a dotted arc, next to
  their partner — as pedigree nomenclature allows. None of the test fixtures
  needs this; densely interbred pedigrees do (305 duplicates for 1,800
  individuals in the test suite's random-mating generator).
  `compute_layout()` returns those duplicates as extra nodes; the input
  family is not modified.

- Visual-regression golden images regenerated. Besides the descent-line fix,
  four multi-phenotype goldens still dated from March and showed the legend
  and gradient bugs fixed in this release; they passed only because the
  comparison tolerated a 2% mean pixel difference.
- Visual regression tests count changed pixels instead of averaging the
  difference over the image. A pixel counts as changed when a channel moves
  by more than 0.1, and a render fails above 0.01% changed pixels or on any
  size change (sizes up to 10% off used to be cropped and compared). Every
  past regression checked had passed the old 2% mean threshold, the smallest
  moving the mean by 0.0003; the new criterion catches them all, and a
  single missing proband arrow (175 pixels). Measured noise is 0 pixels,
  across platforms and on matplotlib 3.7 to 3.10. A failure saves the
  actual render and a diff image with the changed pixels in red.
- The PGS golden test no longer widens the gaps itself: `render_pedigree`
  does it, so the golden encoded doubled spacing.
- Ten new golden images: carrier notation (NSGC 2022, `cr_` columns,
  legacy), the circular layout (rings, twins, duplicates, bands), separate
  sibships, a second partnership and the uncle-niece loop.

- `numpy` is now a declared dependency (`numpy>=1.24`). The renderer imports
  it directly but relied on matplotlib to pull it in. The HTML documentation
  already listed it as a requirement.
- CI now also runs the suite on Python 3.10 with every dependency at its
  declared minimum version (matplotlib 3.7, click 8.0, numpy 1.24, pandas 1.5),
  so the lower bounds in `pyproject.toml` are tested rather than assumed.
- The `python -m pypedigree` test runs the current interpreter instead of
  `uv run`. Going through `uv run` re-synced the project's `.venv` to the
  lockfile in the middle of the test session, which removed pandas and pytest
  for every later test when the environment held other versions, and made the
  suite require uv.

- Reference images under `tests/data/` and `docs/imgs/examples/` regenerated
  for the legend and spacing fixes above. Several of them (`annotations.png`,
  `multi_phenotype_inline.png`, `test_all.png`, `large_family.png`,
  `real_world/multi_phenotype_pgs.png`, `synesthesia.png`) were already stale
  from an earlier legend revision and are now current.

- `__version__` is read from the installed package metadata instead of being
  duplicated in `__init__.py`, so `pyproject.toml` is the single source.
- The circular layout types the edge of a single-child union as
  `EdgeType.DESCENT` instead of `EdgeType.INDIVIDUAL`, matching what the
  hierarchical layout emits for the same structure. Rendering is unaffected
  (the renderer does not branch on edge type), but callers filtering on
  `EdgeType.DESCENT` now see parent-child links in both layouts.
- Sex-inconsistency warnings now carry `PedigreeDataWarning` rather than a
  bare `UserWarning`.
- Internal: the diverging-gradient color computation was duplicated between
  the bands and the legend bar; both now call `_diverging_color_at()`, so the
  two cannot drift apart.
- Internal: `_draw_legend()` went from 472 lines to about 90, with the panel
  geometry in a `_LegendMetrics` dataclass and one function per row type.
  Verified output-identical across 159 renders.

- Source distribution shrank from 4.5 MB to 328 kB. It carried the 2 MB
  documentation PDF, all reference and example images, `CLAUDE.md` and
  `uv.lock`. It now ships source, tests, the visual-regression golden images
  and the Markdown docs — enough to build and run the full suite, verified in
  CI on every release.
- `Development Status` classifier raised from `3 - Alpha` to `4 - Beta`.
- Project URLs, README badge and documentation now point at
  `github.com/darkxander/Pypedigree`.

- **`requires-python` lowered from `>=3.12` to `>=3.10`.** The 3.12 floor was
  arbitrary: the full suite passes on 3.10 and 3.11, and the only 3.10-specific
  feature in the code is `zip(strict=True)`. This matters for bioinformatics
  clusters, many of which still ship 3.10 or 3.11. The documentation already
  claimed 3.10+, so `pyproject.toml` was the outlier.
- CI now tests Python 3.10, 3.11, 3.12, 3.13 and 3.14 on Linux and macOS (10
  jobs, up from 4). 3.14 is the version pinned in `.python-version` for
  development and had never been exercised in CI.
- `ruff` target version lowered to `py310` to match.

### Fixed

- A duplicated individual (someone with several partners, or a partner drawn
  twice to break a loop) is drawn with the original's phenotype bands. It
  was a plain symbol, so in a figure with bands it lost every phenotype.
- The interactive SVG's default size keeps the drawing's proportions. Width
  and height were clamped separately (to 600-1600 and 400-1200 px), so a
  wide pedigree was letterboxed with white bands above and below. Small
  pedigrees are still at least 600 px wide.
- Circular layout: twin connectors and duplicate-individual arcs are drawn,
  as in the hierarchical layout. A union whose only children are one twin
  group has no sibship arc, only the twins' V.
- Band text (values, "NA") on a rotated shape is never drawn upside down.
- An affected carrier lost its affected fill: the shape was drawn white
  before the carrier half-fill. With the NSGC 2022 notation both are shown;
  the legacy notation is unchanged.
- The interactive SVG was drawn upside down, founders at the bottom: the
  layout's y axis points up and SVG's points down. Coordinates are now
  mirrored before rendering, duplicate arcs included.
- Interactive SVG of a file with several families: every family was drawn at
  the same place, on top of the others. They are now side by side, each
  titled, as in the PNG.
- Interactive SVG: element ids were built from sample IDs, so two families
  sharing an ID produced duplicate ids (one clip path served both), and IDs
  with spaces made invalid references. Ids are now generated.
- Interactive SVG: the proband arrow had no head (its marker was never
  defined); the tooltip background was a large black disc (corner radius 3
  layout units); the tooltip ignored zoom and pan and could leave the
  visible area.
- Interactive SVG: text was positioned with `dominant-baseline`, which
  librsvg and older Inkscape ignore; positions are now explicit offsets.
- Hierarchical layout: sibship lines of different couples overlapped, so the
  children of two couples hung from one continuous line and could not be
  told apart. It happened in one synthetic tree pedigree in ten between 8
  and 15 individuals, and in almost every pedigree of 50 or more (1,280
  overlapping pairs over 12 pedigrees of 100-200). Sibship lines no longer
  overlap; `TestSibshipSeparation` checks it on every fixture and on
  generated pedigrees up to 200 individuals, with and without loops.
- Hierarchical layout: the duplicate drawn for a second partnership stood
  on the founders' row whatever the generation of the individual it
  repeats. In `multiple_partners_complex.ped` the father appeared in
  generations I and II and his third child one generation too high.
- Hierarchical layout: a child pulled to a lower generation by their
  partner (`consanguinity_uncle_niece.ped`) was connected by a line that
  stopped one row above them. It now reaches their row.
- Legend panel sized from the font instead of the figure. It took a fixed
  share of the figure height, its swatches were sized on the figure width
  and its text in points: a wide pedigree got huge swatches and tiny text
  (about 3 pt once the 45-individual test pedigree is scaled to a printed
  page), a circular one a panel stretched over the whole square, a legend
  with four rows shrank its text to 4 pt, and raising `font_size` made the
  column captions overlap. One legend unit is now a fixed multiple of
  `font_size`: text, swatches and rows keep their proportions and the same
  physical size in every figure, the border hugs the entries instead of
  spanning the figure, and the figure grows downward to make room rather
  than squeezing the pedigree.
- Duplicate arcs bow by at most 0.4 generation gaps. Their curvature grew
  with their length, so an arc across a wide pedigree swept over several
  generations and, through matplotlib's autoscaling, shrank the whole
  drawing inside a mostly empty figure.
- Hierarchical layout: a single child that could not be centred under its
  parents was joined to them by a diagonal line. Pedigree nomenclature uses
  horizontal and vertical lines only; the child is now reached through a
  short horizontal segment at sibship level, as siblings are. This affected
  10 of the test fixtures, including `three_gen.ped` and both consanguinity
  examples. Node positions are unchanged — only the connecting lines are.
- Hierarchical layout: when a twin group was a union's only offspring and its
  convergence point was off the parents' midpoint, the descent line stopped
  beside it, leaving a gap. The same horizontal segment now joins them.
- Hierarchical layout: a married twin's partner could be drawn between the
  two twins (the side came from the order of the partnership, e.g. the wife
  of the left twin went to his right), so the twin connector spanned across
  the spouse. The partner of the leftmost twin of a group is now drawn on
  its left, and that of the rightmost twin on its right; a middle triplet,
  a twin's second partner, and a partner from a consanguineous loop keep the
  previous placement. `TestTwinSpouses` checks that no one is drawn between
  twins, on dedicated cases and on every fixture. No fixture or golden image
  changes.
- Header columns are matched by name, as documented, even when the header
  has no extra phenotype, annotation or notation column. Such headers were
  read by position, so a structural column outside the standard order was
  misread without a warning: a `Twins` column right after `Phenotype` was
  taken as the date of birth and the twins were lost, a `Deceased` column
  with no `Twins` column before it was parsed as twins, and reordered core
  columns (`IID FID Sex ...`) swapped sex, parents and IDs. Headers in the
  standard order parse exactly as before (every fixture is unchanged).
- Circular layout: ID labels were hidden by their own shape (13 of the 45
  individuals of the large test pedigree), and the proband "P" was drawn
  sideways with an arrow pointing at an unrotated corner, because shapes were
  rotated to face outward. Shapes now stay upright and labels sit radially
  outside them (see Changed). For nodes an API caller rotates, the label
  still moves below the rotated outline, with a white backing, and the
  proband marker stays upright, pointing at the rotated outline.

- Categorical phenotype colors (`c_` columns) are now resolved once per figure
  instead of being accumulated in a module-level cache. This fixes three defects:
  - Two categories could be assigned the same color when the column's own
    `phenotype_colors` entry also appeared in the categorical palette (the two
    palettes share six colors). In `prefix_modes.ped`, `ASD` and `ADHD` were
    both rendered green.
  - Categories leaked between figures: rendering family A then family B in one
    process listed A's categories in B's legend.
  - A second render with a different `phenotype_colors` silently reused the
    colors of the first, so changing a color in a notebook had no effect.
  - Category colors now also avoid the colors of other phenotype columns, so
    two bands on the same individual stay visually distinct.
- `render_pedigree()` now produces the same figure as the CLI. Palette
  auto-assignment lived only in `cli.py` and `PedigreeBuilder`, so
  `render_pedigree(parse_pedigree("multi_pheno.ped"))` silently rendered plain
  shapes with no bands and no legend. The logic moved to
  `resolve_style()` in `config.py`, which every entry point now shares.
- PGS spacing is no longer applied twice. Both the CLI and `render_pedigree()`
  widened `couple_gap`/`sibling_gap`, so figures from the CLI with a column
  named exactly `PGS` were spaced twice as far apart as intended.
  `apply_pgs_spacing()` is now applied by `render_pedigree()` only.
- `PedigreeBuilder.render()` no longer mutates the style passed to
  `set_style()`, so rendering the same builder twice no longer compounds
  PGS spacing or leaves auto-assigned colors behind.
- The legend now matches the bands for AUTO (unprefixed) phenotype columns.
  The legend had no AUTO branch and fell through to the binary
  Affected/Unaffected/Unknown swatches, while the bands were drawn as
  gradients. An unprefixed continuous column such as `IQ` showed gradient
  bands under a binary legend.
- The effective display mode of a column is decided once from all of its
  values instead of per individual, so a column can no longer mix binary and
  gradient bands across the same diagram (e.g. values `0,1,2,5,17`).
- The legend no longer reserves rows for `PGS`-prefixed columns as if they
  were bands, which left an oversized empty legend box.
- The Affected/Unaffected/Unknown legend header is only drawn when at least
  one band is actually binary.
- CLI validation warnings are reported as `Warning: <message>` on stderr
  instead of raw Python `UserWarning` output carrying internal file paths and
  an unrelated source line.
- Sex-consistency warnings are reported once per parent instead of once per
  child of that parent.
- `docs/examples_gallery.md`: `--style` takes a path to a JSON file; the
  inline-JSON example did not work.

### Removed

- `tests/visual/golden/multi_phenotype_inline.actual.png`: a visual-diff
  failure artifact committed by accident, sitting next to the golden images it
  could be mistaken for. `*.actual.png` is now gitignored.

- `NodeLayout.phenotype_display_modes`: a field neither layout ever populated
  and nothing ever read. The live field is `Pedigree.phenotype_display_modes`,
  which the renderer uses. Removing the dead twin stops it misleading readers.

### Documented

- A larger gallery example: `tests/data/complex/multigenerational.ped`, a
  synthetic 81-person, five-generation pedigree with a first-cousin
  marriage and a marriage between two branches, both drawn side by side
  without duplicates, plus bands, carriers, twins and standard notation. It
  is written by `generate_fixtures.py` and has a visual-regression golden.
  It replaces the 100-individual gallery image, which predated the current
  layout (crossing lines, unconnected spouses).
- The examples gallery no longer describes `synesthesia.ped` as data from a
  real study: like every fixture, it is synthetic.

- The interactive SVG does not support the circular layout. Stated in the
  README, the API reference and the HTML/PDF documentation.

### Security

- **Test fixtures are now fully synthetic.** Five fixtures previously held
  data derived from real participants: pedigree topologies plus diagnoses, IQ,
  SRS scores, polygenic scores and gene-level findings. De-identification alone
  was not enough — an accession code is a pseudonym rather than anonymisation,
  and a large pedigree carrying a rare phenotype can be identifying on its own
  — so the fixtures were regenerated from scratch:
  - `tests/data/generate_fixtures.py` builds them from a fixed seed, so they are
    reproducible and demonstrably fabricated.
  - Feature coverage is unchanged: the same column sets, display-mode prefixes,
    `-9` no-data values, comma decimal separators, generation shape (2 / 20 / 23
    individuals) and union count. No test assertion needed changing.
  - The git history was rewritten with `git filter-repo` so the two cohort
    accession codes appear in no commit, message or path.
  - `tests/data/real_world/` is renamed `tests/data/complex/`; the old name
    described provenance that no longer exists and would have misled the next
    reader into handling the files as real data.
  - `tests/data/README.md` records the provenance rule, and
    `tests/test_fixtures.py` enforces it: it fails on any accession-shaped
    identifier, any populated date-of-birth column, any family ID outside the
    synthetic set, and on a generated fixture that has drifted from the
    generator output.
- `CONTRIBUTING.md` and all three issue/PR templates now require test data to
  be synthetic or fully de-identified.

### Known limitations

- The interactive SVG has no circular layout: `--circular` with
  `--interactive` warns and draws the hierarchical layout.
- Gamete donors and surrogates are marked "D" and "S", but the diagonal
  lines that join them to a pregnancy in the NSGC figures are not drawn: a
  donor is a parent in the file and is joined to the other parent like any
  couple. Not covered at all: several individuals in one symbol ("n"),
  documented evaluation ("*"), the asymptomatic-carrier vertical line, and
  childlessness of a couple (the marker hangs from one individual).
- A pregnancy-loss triangle has no room for phenotype bands; in a figure
  with bands it takes the color of the individual's first band.

## [0.3.0] - 2026-04-03

### Added

- Annotation support:
  - Columns with `a_` prefix in .ped header are annotation columns (text displayed below individual's ID label)
  - Values rendered as italic grey text; missing values (0, -9, NA, empty) produce no annotation line
  - Generation gap auto-increases when annotations are present
  - External annotation file via `--annotations PATH` CLI option
  - `parse_annotation_file(path)` and `apply_annotations_to_pedigree()` API functions
  - `Individual.annotations: dict[str, str]` field on the data model
  - `NodeLayout.annotations: list[str]` field on the layout model
  - `Pedigree.annotation_names: list[str]` field
  - `PedigreeBuilder.add_individual(..., annotations=["line1", "line2"])` parameter
- Gradient inverted display mode:
  - New `gi_` prefix for `GRADIENT_INVERTED` mode
  - Color-to-white scale (high value = white, low value = full color)
  - Useful for metrics where low = bad (e.g., IQ: low IQ = more red)
- Multiple PGS columns:
  - Any column starting with `PGS` (case-insensitive) is now treated as a PGS column (e.g. `PGS_autism`, `PGS_IQ`, `PGS_ADHD`)
  - Single gaussian curve with multiple colored markers
  - `pgs_colors` style option for per-column marker colors
  - `couple_gap` and `sibling_gap` auto-increase when PGS columns are present
- New `PedigreeStyle` fields:
  - `pgs_colors: dict[str, str]` -- PGS column name to hex color for marker line
  - `gradient_range: dict[str, list[float]]` -- Custom min/max for gradient normalization
  - `gradient_diverging: dict[str, dict]` -- Tricolor diverging gradient configuration
  - `zero_pad_ids: int` -- Zero-pad sample IDs in labels to this width (0 = no padding)
- Comma decimal separator support in parser (European locale: `0,712` parsed as `0.712`)
- Gradient band values always displayed regardless of band color intensity
- AUTO mode now uses raw values (1/2) for binary detection instead of normalized values

## [0.2.2] - 2026-03-20

### Added

- Circular/radial layout mode:
  - `compute_layout_circular(family, style)` places founders at the center and descendants on concentric rings
  - `--circular` CLI flag to activate radial layout
  - `circular=True` parameter on `render_pedigree()` and `PedigreeBuilder.render()`
  - `radial_gap` and `base_radius` style parameters in `PedigreeStyle` for tuning ring spacing
  - `waypoints` field on `EdgeLayout` for curved/polyline edge routing
  - Waypoint-based polyline drawing in the renderer for non-orthogonal connections

### Changed

- Circular layout algorithm completely rewritten: sector-based radial model with per-lineage angular sectors, auto-scaling radius, and proper pedigree structure (DESCENT + SIBSHIP arcs + INDIVIDUAL lines instead of direct star-pattern lines).
- `base_radius` default changed from 2.0 to 3.0 in `PedigreeStyle`.
- Extra phenotype value semantics: `0` = unknown (grey band), `-9` = no data (white band with "NA" text). When all bands are NA for an individual, a single centered "NA" is shown instead of repeating per band.
- Legend layout is now fully adaptive: row height, swatch sizes, font sizes, and gaps all scale based on the number of phenotype rows (no more hardcoded sizes).
- Legend text simplified: removed "Top band:" / "Middle band:" prefixes, now shows just the phenotype name.
- Legend column headers ("Affected" / "Unaffected" / "Unknown") now shown once above all phenotype rows instead of repeated per row, preventing text overlap and reducing clutter.
- Legend text positioning fully dynamic: `text_x` (swatch rows) computed from actual swatch widths, `icon_text_x` (PGS, carrier, outline rows) computed from actual icon widths — both with proportional gap that scales with row count.
- Added automated legend overlap tests (`TestLegendNoOverlap`) that verify no text bounding boxes collide.
- `synesthesia.ped` updated: uses `b_ASD`, `b_Synesthesia`, `o_Psychiatric comorbidities` prefixed columns, IDs without leading spaces, and includes Deceased and Proband columns.

## [0.2.1] - 2026-03-19

### Added

- Smart header-based column mapping:
  - Parser now recognises column names (FID, ID, Father, Mother, Sex, ASD, etc.) and maps them to their roles
  - Extra columns not matching known structural names are treated as phenotype columns
  - Works with any column order (e.g., `FID ID Father Mother Sex ASD Synesthesia`)
  - No longer requires 10 fixed-position columns before phenotypes
  - Supports column name aliases (e.g., `Pat`/`Father`/`Dad`, `Pheno`/`Phenotype`/`Affected`)
  - First extra phenotype column used as primary phenotype when no explicit `Phenotype` column present
- Synesthesia pedigree test data (45 individuals, 2 phenotypes: ASD, Synesthesia)
- Standard pedigree encoding support for extra phenotypes (1=unaffected, 2=affected, 0/-9=unknown)
- Binary phenotype rendering: color=affected, white=unaffected, grey=unknown (no text by default)
- `show_phenotype_labels` style option to show YES/NO text inside binary bands
- `--show-phenotype-labels` CLI flag
- Prefix-based phenotype display modes in column headers:
  - `b_` → binary (color/white/grey bands)
  - `g_` → gradient (white→color continuous scale with value labels)
  - `o_` → outline (colored shape border when affected)
  - `c_` → categorical (distinct color per unique string value)
  - No prefix → auto-detect from data (backward compatible)
- `PhenotypeDisplayMode` enum exported from pypedigree

## [0.2.0] - 2026-03-19

### Added

- Carrier notation support:
  - Half-filled shapes for heterozygous carriers (`is_carrier`)
  - Center dot for obligate carriers (`is_obligate_carrier`)
  - `carrier` column in header-based .ped format (values: 0/1/2)
  - `carrier` and `obligate_carrier` parameters in `PedigreeBuilder.add_individual()`
  - Auto-generated legend entries for carrier symbols
- Interactive SVG export:
  - New `svg_renderer.py` module producing standalone SVG with embedded JavaScript
  - Tooltips on hover (sample_id, sex, phenotype, DOB, extra phenotypes)
  - Click to highlight connected relatives (parents, children, siblings)
  - Mouse wheel zoom and Ctrl+drag pan for large pedigrees
  - `--interactive` CLI flag
  - `render_interactive_svg()` public API function
  - `PedigreeBuilder.render_interactive_svg()` method
  - Carrier notation rendered in SVG (clip-path half-fill and center dot)

## [0.1.0] - 2026-03-16

### Added

- Initial release
- Standard 6-column pedigree file parsing (TSV, CSV, space-separated)
- Extended format support (DOB, twins, deceased, proband)
- Multi-phenotype rendering with horizontal color bands
  - Binary phenotypes (0/1): filled or empty bands with YES/NO labels
  - Continuous phenotypes: color gradient bands with value labels
- PGS (Polygenic Score) gaussian distribution curves per individual
- Inline phenotype parsing from columns 11+ (with header row)
- Separate `.pheno` TSV file support via `--pheno` CLI option
- Auto-generated phenotype legend with proper male/female symbols
- JSON-based style customization (`--style`)
- Multiple output formats: PNG, SVG, PDF
- Interactive display mode (`--show`)
- Family filtering (`--family`)
- Twins support: monozygotic, dizygotic, unknown zygosity
- Deceased markers, proband arrows, date of birth labels
- Consanguineous relationship rendering
- Missing/inferred parent handling (dashed shapes)
- Pandas DataFrame integration (`parse_pedigree_from_dataframe`)
- CLI options: `--outline-color`, `--no-labels`, `--no-legend`, `--dpi`
- Programmatic `PedigreeBuilder` API with fluent interface and Jupyter support
- Visual regression tests with golden image comparison (`--regen-golden`)
- Large pedigree layout (100-200+ individuals):
  - O(n) generation assignment with parent-child adjacency index
  - Barycenter heuristic for edge crossing reduction
  - Bidirectional overlap resolution (left-to-right + right-to-left pullback)
  - Bottom-up parent centering over children
  - Backward-compatible: new passes only activate for >15 individuals
- Auto-scaling figure dimensions based on layout bounds
- Comprehensive HTML and PDF documentation (`docs/`)
- PLINK bioinformatics integration:
  - `parse_fam_file()` for native PLINK `.fam` reading
  - `parse_plink_ped()` for PLINK `.ped` files (ignores genotype columns)
  - `export_fam()` for exporting pedigrees to `.fam` format (PLINK/KING interop)
  - CLI `--plink` flag and auto-detection of `.fam` extension
