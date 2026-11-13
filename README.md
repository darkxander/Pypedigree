# pypedigree

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/darkxander/Pypedigree/actions/workflows/ci.yml/badge.svg)](https://github.com/darkxander/Pypedigree/actions/workflows/ci.yml)
[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)

Generate publication-quality pedigree (family tree) diagrams from standard pedigree files used in genetics research.

## Features

- **Standard pedigree symbols**: squares (male), circles (female), diamonds (unknown sex), with affected/unaffected fill
- **Multi-phenotype rendering**: horizontal color bands for multiple phenotypes per individual (binary: 0=unknown/grey, -9=NA/white; continuous: gradient). The legend always matches the bands it describes; it is sized from the font, so it reads the same under a 4-person pedigree as under a 45-person one.
- **PGS (Polygenic Score) curves**: gaussian distribution with individual score marker, drawn alongside each node. Supports multiple PGS columns (e.g. `PGS_autism`, `PGS_IQ`) with colored markers on a single curve.
- **Annotations**: columns with `a_` prefix in .ped header display italic grey text below individual ID labels. Also loadable from external file via `--annotations`.
- **Twins support**: monozygotic, dizygotic, and unknown zygosity with correct connectors
- **Extended metadata**: deceased markers, proband arrows, date of birth labels
- **Standard notation (NSGC 2008/2022)**: pregnancies ("P"), miscarriages, ectopic pregnancies and terminations (triangles), stillbirths, adoption (brackets, dashed line to adoptive parents), no children by choice and infertility, gamete donors ("D") and surrogates ("S"), gender identity (the symbol follows gender, AMAB/AFAB written below), consultand (arrow without "P"). One column each; see [Standard notation columns](#standard-notation-columns)
- **Carrier notation (NSGC 2022)**: the symbol is divided into patterned sections, one per carrier result, so an affected carrier or a carrier of several conditions shows every result; `cr_` columns give one carrier status per gene or condition. The pre-2022 half-fill / center-dot notation stays available (`carrier_notation: "legacy"`)
- **Interactive SVG**: the same figure as the PNG (phenotype bands, legend, PGS curves, annotations, carrier notation) plus tooltips, click-to-highlight relatives and zoom/pan for browser viewing
- **Multiple output formats**: PNG, SVG, PDF, interactive SVG
- **Customizable styling**: JSON-based style overrides for colors, sizes, fonts
- **Separate phenotype files**: load phenotypes from `.pheno` TSV files
- **Optional legend**: auto-generated phenotype legend with proper male/female symbols
- **Programmatic API**: fluent `PedigreeBuilder` for constructing pedigrees from code
- **Circular/radial layout**: founders at center, descendants on concentric rings (`--circular` flag), in the order of the hierarchical layout, so sibships never interleave; ring radius auto-expands so nodes never overlap; upright shapes, labels set radially outside them, twin connectors and duplicate arcs as in the default layout
- **Large pedigree support**: each couple's children hang from their own stretch of sibship line, so sibships never merge, at 10 or 500 individuals; loops (consanguinity, marriages between branches) are ordered so the partners meet side by side, and drawn with a duplicate individual only where they would still tangle lines
- **Auto-scaling figures**: figure dimensions grow with pedigree size
- **Jupyter integration**: inline display via `_repr_png_`
- **PLINK interoperability**: native `.fam` parsing, `.fam` export, PLINK `.ped` support (genotype columns ignored)

## Installation

```bash
pip install pypedigree
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv add pypedigree
```

With optional pandas support:

```bash
pip install pypedigree[pandas]
```

## Quick Start

### Command Line

```bash
# Basic PNG output
pypedigree family.ped -o family.png

# SVG for publications
pypedigree family.ped -o family.svg

# With multi-phenotype colors from a style file
pypedigree family.ped --style style.json -o family.png

# With external phenotype data
pypedigree family.ped --pheno phenotypes.tsv --style style.json -o family.png

# With external annotations
pypedigree family.ped --annotations annotations.tsv -o family.png

# Render only one family from a multi-family file
pypedigree multi_family.ped --family FAM1 -o fam1.png

# High resolution, no legend
pypedigree family.ped --dpi 300 --no-legend -o family.png

# Interactive display
pypedigree family.ped --show

# Interactive SVG with tooltips and click highlighting
pypedigree family.ped --interactive -o family.svg

# Circular/radial layout
pypedigree family.ped --circular -o family_radial.png
```


> `--interactive` draws the same figure as the PNG/PDF output — phenotype
> bands, outline colors, PGS curves, `a_` annotations, carrier notation and
> the legend (`--no-legend` hides it) — and adds a tooltip per individual
> listing every phenotype value. Several families are drawn side by side.
> Only the circular layout (`--circular`) is not available in interactive
> mode.

### Python API

`render_pedigree()` produces the same figure as the CLI, including
auto-assigned phenotype colors and the legend — no extra setup needed.

```python
from pypedigree import parse_pedigree, render_pedigree, PedigreeStyle

# Parse and render
pedigree = parse_pedigree("family.ped")
fig = render_pedigree(pedigree, output="family.png")

# Custom style
style = PedigreeStyle(dpi=300, fill_color="red")
fig = render_pedigree(pedigree, style=style, output="custom.png")

# Circular layout
fig = render_pedigree(pedigree, circular=True, output="radial.png")

# With phenotype colors
style = PedigreeStyle()
style.phenotype_colors = {"epilepsy": "#e74c3c", "IQ": "#3498db"}
fig = render_pedigree(pedigree, style=style, output="pheno.png")
```

### Programmatic Builder

```python
from pypedigree import PedigreeBuilder

builder = PedigreeBuilder("FAM1")
builder.add_individual("father", sex="male", carrier=True)
builder.add_individual("mother", sex="female")
builder.add_individual("child", sex="male", father="father", mother="mother",
                       phenotype="affected", phenotypes={"IQ": 85, "epilepsy": 1})
builder.add_individual("sibling", father="father", mother="mother",
                       pregnancy="SAB", annotations=["8 wk"])
fig = builder.render(output="family.png")

# Interactive SVG output
svg = builder.render_interactive_svg(output="family.svg")
```

### PLINK Integration

```python
from pypedigree import parse_fam_file, parse_plink_ped, export_fam, render_pedigree

# Read PLINK .fam file
pedigree = parse_fam_file("family.fam")
render_pedigree(pedigree, output="family.png")

# Read PLINK .ped (genotype columns are ignored)
pedigree = parse_plink_ped("family.ped")

# Export to .fam for use with PLINK/KING/etc.
export_fam(pedigree, "output.fam")
```

```bash
# CLI: .fam files are auto-detected
pypedigree family.fam -o family.png

# Force PLINK mode for .ped files with genotype data
pypedigree family.ped --plink -o family.png
```

### Pandas Integration

```python
import pandas as pd
from pypedigree import parse_pedigree_from_dataframe, render_pedigree

df = pd.read_csv("family.ped", sep="\t")
pedigree = parse_pedigree_from_dataframe(df)
render_pedigree(pedigree, output="family.png")
```

## File Formats

### Pedigree File (.ped)

Tab, comma, or space-separated with 6+ columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | family_id | Family identifier |
| 2 | sample_id | Individual identifier |
| 3 | paternal_id | Father's sample_id (0 = unknown) |
| 4 | maternal_id | Mother's sample_id (0 = unknown) |
| 5 | sex | 1=male, 2=female, 0=unknown |
| 6 | phenotype | 1=unaffected, 2=affected, 0=unknown |
| 7 | date_of_birth | YYYYMMDD / YYYYMM / YYYY (optional) |
| 8 | twins | TYPE:sample_id entries (optional) |
| 9 | deceased | 1=deceased, 0=alive (optional) |
| 10 | proband | 1=proband, 0=not (optional) |
| 11+ | phenotypes | Extra phenotype columns (requires header) |
| 11+ | carrier | Special column: 0=not, 1=carrier, 2=obligate carrier (detected by header name) |
| 11+ | cr_GENE | One carrier column per gene or condition, same 0/1/2 coding (e.g. `cr_CFTR`) |

**Column name prefixes** control how phenotype columns are rendered: `b_` (binary), `g_` (gradient), `gi_` (gradient inverted: color-to-white, for metrics where low = bad), `o_` (outline), `c_` (categorical), `cr_` (carrier status: 0=not, 1=carrier, 2=obligate carrier). Columns with `a_` prefix are annotation columns (text displayed below the individual's ID label, not phenotypes).

A column with no prefix is auto-detected once per column, from all of its values: any string value makes it categorical, any numeric value outside `0`/`1`/`2` makes it a gradient, otherwise it is binary.

**Smart header mapping**: When a header row is present, columns are matched by name (not position). Any unrecognised column name after the core 5 (FID, ID, Father, Mother, Sex) is treated as a phenotype column. This means you can write compact files without the DOB/twins/deceased/proband columns:

```
FID	ID	Father	Mother	Sex	ASD	Synesthesia
FAM1	001	002	003	1	2	2
FAM1	002	0	0	1	1	2
```

Full 10-column format with extras:

```
FID	IID	PAT	MAT	SEX	PHENO	DOB	TWINS	DECEASED	PROBAND	epilepsy	IQ	SRS	PGS
mockup	1	0	0	1	2	0	0	0	0	1	105	62	0.72
mockup	2	0	0	2	1	0	0	0	0	0	118	45	0.35
```

### Standard notation columns

Each of these columns is recognised by name in a header row (case-insensitive;
aliases in brackets). Codes are case-insensitive, and `0`, `-9`, `NA`, `.` or an
empty cell mean nothing to note. An unknown code is an error naming the line.
Symbols follow the NSGC nomenclature (Bennett et al., J Genet Couns 2008 and
its 2022 revision).

| Column | Codes | Drawn as |
|--------|-------|----------|
| `pregnancy` (`pregnancy_outcome`, `outcome`) | `P` ongoing pregnancy, `SAB` miscarriage, `TOP` termination, `ECT` ectopic, `SB` stillbirth | `P`: "P" inside the symbol (a diamond when the sex is unknown). `SAB`, `TOP`, `ECT`: a triangle, filled when affected, struck through for `TOP`, "ECT" below. `SB`: the symbol struck through, "SB" below |
| `adopted` (`adoption`) | `in` (or `1`): the parents given are the adoptive parents; `out` (or `2`): they are the biological parents | Brackets around the symbol; the line to adoptive parents is dashed |
| `childless` (`childlessness`, `infertility`) | `choice` (or `1`): no children by choice or for an unknown reason; `infertility` (or `2`) | A line below the symbol ending in one bar, or two for infertility. Write the reason in an `a_` column |
| `art` (`assisted_reproduction`) | `D` gamete donor, `S` surrogate (also `GC`) | "D" or "S" inside the symbol |
| `gender_identity` (`identity`) | `man`/`M`/`1`, `woman`/`F`/`W`/`2`, `nonbinary`/`NB`/`X`/`3` | The symbol follows the gender identity (square, circle, diamond) and the `sex` column is taken as the sex assigned at birth, written below as AMAB or AFAB when the two differ |
| `consultand` | `1` | An arrow without "P" (the proband's arrow has one) |

A column named `gender` is the sex column, as before, unless the header also
has a `sex` column: then it is the gender identity.

Gestational age, karyotype and similar details go in `a_` annotation columns,
written below the symbol. A pregnancy-loss triangle has no room for phenotype
bands: in a figure with bands it takes the color of the individual's first band.
On a symbol with bands, the "P", "D" or "S" moves beside its upper-right corner.

```
FID	IID	Father	Mother	Sex	Phenotype	Pregnancy	Adopted	Gender_identity	a_Note
FAM1	f	0	0	1	1	0	0	0	0
FAM1	m	0	0	2	1	0	0	0	0
FAM1	c1	f	m	2	1	0	in	0	0
FAM1	c2	f	m	1	1	0	0	woman	0
FAM1	c3	f	m	0	2	SAB	0	0	8 wk
```

Example fixture with every symbol: `tests/data/nomenclature/nomenclature.ped`.
The diagonal lines that connect gamete donors and surrogates in the NSGC
figures are not drawn: a donor is a parent in the file and is joined to the
other parent like any couple.

### PLINK FAM File (.fam)

PLINK `.fam` files use the same 6-column format (FID, IID, PAT, MAT, SEX, PHENO), whitespace-delimited, no header. Files with `.fam` extension are auto-detected by the CLI.

### Style File (.json)

```json
{
    "phenotype_colors": {
        "epilepsy": "#e74c3c",
        "IQ": "#3498db",
        "SRS": "#e67e22"
    },
    "pgs_colors": {
        "PGS_autism": "#e74c3c",
        "PGS_IQ": "#3498db"
    },
    "gradient_range": {
        "SRS": [40, 90]
    },
    "gradient_diverging": {
        "IQ": {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}
    },
    "zero_pad_ids": 3
}
```

### Phenotype File (.pheno)

Tab-separated, first column is sample_id, remaining columns are phenotype values:

```
sample_id	epilepsy	IQ	SRS
1	1	105	62
2	0	118	45
```

## CLI Reference

```
pypedigree INPUT.ped [OPTIONS]

Options:
  -o, --output PATH              Output file path
  --format [png|svg|pdf]         Output format (auto-detected from extension)
  --family ID                    Render only this family ID
  --dpi INT                      Output resolution (default: 150)
  --show                         Display interactively instead of saving
  --style PATH                   JSON file with style overrides
  --separator [tab|comma|space]  Force field separator
  --pheno PATH                   Separate .pheno TSV file with extra phenotypes
  --annotations PATH             Separate TSV file with annotation text (displayed below IDs)
  --outline-color COLOR          Shape outline color (e.g. 'grey', '#888888')
  --carrier-notation [nsgc2022|legacy]  Carrier notation (default: nsgc2022)
  --no-labels                    Hide sample ID labels
  --no-legend                    Hide the phenotype legend
  --show-phenotype-labels        Show YES/NO text labels inside binary phenotype bands
  --interactive                  Generate interactive SVG with tooltips and click highlighting
  --circular                     Use circular/radial layout (founders at center)
  --plink                        Parse as PLINK format (ignore genotype columns)
  --version                      Show version
```

Data-quality warnings (a father recorded as female, a parent missing from the
family, ...) are printed to stderr as `Warning: <message>` and do not stop
rendering. Unknown options in a style file are reported the same way, with a
suggested spelling. An invalid color is a hard error naming the option it came
from.

In Python these are the `PedigreeDataWarning` / `PedigreeStyleWarning`
categories, so they can be filtered or escalated:

```python
import warnings
from pypedigree import PedigreeDataWarning

warnings.simplefilter("error", PedigreeDataWarning)  # fail on questionable data
```

## Diagram Symbols

- **Square** = male, **Circle** = female, **Diamond** = unknown sex
- **Filled shape** = affected, **Empty shape** = unaffected
- **Patterned section** = carrier, one pattern per carrier column, sparser for obligate carriers (NSGC 2022: the symbol is divided into as many sections as there are results, an affected section first). With `carrier_notation: "legacy"`: **half-filled shape** = carrier, **center dot** = obligate carrier
- **Color bands** = multiple phenotypes (binary: color/white/grey/NA; continuous: gradient)
- **Gaussian curve** = PGS distribution with individual score marker
- **Diagonal line** = deceased
- **"P" with arrow** = proband; **arrow without "P"** = consultand
- **"P" inside the symbol** = pregnancy; **triangle** = miscarriage, ectopic pregnancy ("ECT"), or termination (struck through); **struck-through symbol with "SB"** = stillbirth
- **Brackets** = adopted; a **dashed line** joins an adopted child to the adoptive parents
- **Line ending in a bar** below the symbol = no children by choice; **two bars** = infertility
- **"D" / "S" inside the symbol** = gamete donor / surrogate
- **AMAB / AFAB** below the symbol = sex assigned at birth, when the symbol (gender identity) differs
- **Double horizontal line** = consanguineous relationship
- **Dashed shape** = missing/inferred parent
- **Connecting lines** are horizontal or vertical only, as the standard nomenclature requires

## Development

```bash
git clone https://github.com/darkxander/Pypedigree.git
cd pypedigree
uv sync --extra dev

# Run tests
uv run pytest

# Run tests with coverage
uv run pytest --cov=pypedigree --cov-report=term-missing

# Lint & format
uv run ruff check src/ tests/
uv run ruff format src/ tests/
```

## Citing pypedigree

If you use pypedigree in your research, please cite it. Citation metadata is
in [CITATION.cff](CITATION.cff); GitHub renders it as a **Cite this
repository** button.

## Releasing

See [docs/RELEASING.md](docs/RELEASING.md). Releases are published to PyPI
automatically from a version tag, using trusted publishing (no API token).

## License

MIT
