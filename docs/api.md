# API Reference

## Parsing

### `parse_pedigree(file_path)`

Parse a pedigree from a file. Auto-detects separator (tab, comma, space) and header rows. When a header is present, columns are matched by name (not position), so phenotype columns can appear directly after Sex without needing the full 10-column format.

```python
from pypedigree import parse_pedigree

ped = parse_pedigree("family.ped")
```

**Args**: `file_path` — Path to the pedigree file (str or Path).

**Returns**: `Pedigree` object containing all families.

**Raises**: `PedigreeParseError`, `PedigreeValidationError`, `FileNotFoundError`.

---

### `parse_pedigree_from_string(text)`

Parse a pedigree from a string (same format as file contents).

```python
from pypedigree import parse_pedigree_from_string

text = "FAM1\tfather\t0\t0\t1\t1\nFAM1\tmother\t0\t0\t2\t1\nFAM1\tchild\tfather\tmother\t2\t2\n"
ped = parse_pedigree_from_string(text)
```

**Args**: `text` — Pedigree data as a string.

**Returns**: `Pedigree` object.

---

### `parse_pedigree_from_dataframe(df)`

Parse a pedigree from a pandas DataFrame. Requires `pip install pypedigree[pandas]`.

Column position determines meaning (names are ignored):
1. family_id, 2. sample_id, 3. paternal_id, 4. maternal_id, 5. sex, 6. phenotype,
7. date_of_birth (opt), 8. twins (opt), 9. deceased (opt), 10. proband (opt).

```python
import pandas as pd
from pypedigree import parse_pedigree_from_dataframe

df = pd.DataFrame({
    "fid": ["FAM1", "FAM1", "FAM1"],
    "iid": ["father", "mother", "child"],
    "pat": ["0", "0", "father"],
    "mat": ["0", "0", "mother"],
    "sex": [1, 2, 2],
    "pheno": [1, 1, 2],
})
ped = parse_pedigree_from_dataframe(df)
```

**Args**: `df` — pandas DataFrame with at least 6 columns.

**Returns**: `Pedigree` object.

---

### `parse_fam_file(file_path)`

Parse a PLINK `.fam` file. The `.fam` format is whitespace-delimited with exactly 6 columns (FID, IID, PAT, MAT, SEX, PHENO) and no header row.

```python
from pypedigree import parse_fam_file

ped = parse_fam_file("family.fam")
```

**Args**: `file_path` — Path to the `.fam` file (str or Path).

**Returns**: `Pedigree` object.

---

### `parse_plink_ped(file_path)`

Parse a PLINK `.ped` file, ignoring genotype columns after column 6. Only the standard 6 pedigree columns are read; no header row is expected.

```python
from pypedigree import parse_plink_ped

ped = parse_plink_ped("family.ped")
```

**Args**: `file_path` — Path to the PLINK `.ped` file (str or Path).

**Returns**: `Pedigree` object.

---

### `export_fam(pedigree, file_path)`

Export a `Pedigree` to PLINK `.fam` format. Writes tab-delimited, 6 columns, no header. Phantom parents and duplicate nodes are excluded; internal references (`_dup_`, `_phantom_`) are resolved back to original IDs.

```python
from pypedigree import parse_pedigree, export_fam

ped = parse_pedigree("family.ped")
export_fam(ped, "output.fam")
```

**Args**:
- `pedigree` — `Pedigree` to export.
- `file_path` — Output file path (str or Path).

**Returns**: `Path` of the written file.

---

### `parse_phenotype_file(file_path)`

Parse a separate phenotype file (.pheno TSV). First row is a header (`sample_id`, `pheno1`, `pheno2`, ...), subsequent rows are data.

```python
from pypedigree import parse_phenotype_file

pheno_data, pheno_names = parse_phenotype_file("data.pheno")
# pheno_data: {"1": {"epilepsy": 1.0, "IQ": 105.0}, ...}
# pheno_names: ["epilepsy", "IQ", "SRS"]
```

**Args**: `file_path` — Path to the .pheno file (str or Path).

**Returns**: Tuple of `(pheno_data, pheno_names)` where `pheno_data` maps `sample_id -> {pheno_name: raw_value}`.

---

### `apply_phenotypes_to_pedigree(pedigree, pheno_data, pheno_names)`

Apply external phenotype data to individuals in a pedigree and normalize values.

```python
from pypedigree import parse_pedigree, parse_phenotype_file, apply_phenotypes_to_pedigree

ped = parse_pedigree("family.ped")
pheno_data, pheno_names = parse_phenotype_file("data.pheno")
apply_phenotypes_to_pedigree(ped, pheno_data, pheno_names)
```

**Args**:
- `pedigree` — `Pedigree` to update (modified in place).
- `pheno_data` — Map of `sample_id -> {pheno_name: raw_value}`.
- `pheno_names` — Ordered list of phenotype names.

---

### `parse_annotation_file(file_path)`

Parse an annotation file (TSV). First row is a header (`sample_id`, followed by annotation column names), subsequent rows are data. Annotation columns provide text displayed below the individual's ID label (not phenotypes).

```python
from pypedigree import parse_annotation_file

anno_data, anno_names = parse_annotation_file("annotations.tsv")
# anno_data: {"1": {"Role": "Proband", "Clinic": "Paris"}, ...}
# anno_names: ["Role", "Clinic"]
```

**Args**: `file_path` -- Path to the annotation file (str or Path).

**Returns**: Tuple of `(anno_data, anno_names)` where `anno_data` maps `sample_id -> {anno_name: value}` and `anno_names` is the ordered list of annotation column names.

---

### `apply_annotations_to_pedigree(pedigree, anno_data, anno_names)`

Apply external annotation data to individuals in a pedigree.

```python
from pypedigree import parse_pedigree, parse_annotation_file, apply_annotations_to_pedigree

ped = parse_pedigree("family.ped")
anno_data, anno_names = parse_annotation_file("annotations.tsv")
apply_annotations_to_pedigree(ped, anno_data, anno_names)
```

**Args**:
- `pedigree` -- `Pedigree` to update (modified in place).
- `anno_data` -- Map of `sample_id -> {anno_name: value}`.
- `anno_names` -- Ordered list of annotation column names.

### `build_notation(*, pregnancy=None, adopted=None, childless=None, art=None, gender=None, consultand=False)`

Build a `Notation` from enum values or the codes of the `.ped` notation columns
(case-insensitive): `pregnancy` `P`/`SAB`/`TOP`/`ECT`/`SB`, `adopted` `in`/`out`,
`childless` `choice`/`infertility`, `art` `D`/`S`, `gender` `man`/`woman`/`nonbinary`.
Raises `PedigreeParseError` on an unknown code. `PedigreeBuilder.add_individual()`
uses it.

```python
from pypedigree import build_notation

build_notation(pregnancy="SAB", adopted="in")
# Notation(pregnancy=<Pregnancy.MISCARRIAGE: 'SAB'>, adoption=<Adoption.IN: 'in'>, ...)
```

---

## Layout

### `compute_layout(family, style=None)`

Compute (x, y) coordinates for all individuals and connection geometry. No matplotlib dependency.

```python
from pypedigree import compute_layout, PedigreeStyle

layout = compute_layout(family, style=PedigreeStyle())
```

**Args**: `family` — `Family` object. `style` — `PedigreeStyle` (optional, defaults used if None).

**Returns**: `FamilyLayout` with `.nodes`, `.edges`, `.twin_connectors`, `.duplicate_arcs`.

#### Layout pipeline

1. **Assign generations** — BFS from founders using a parent→child index, then partners are aligned and children propagated until every child sits below both parents. A duplicate node stands on the generation of the individual it repeats.
2. **Compute horizontal positions** — Each union's descendants form a block. Sibling blocks are packed left to right against each other's contours, generation by generation, and every couple is centred over its children. Nodes on a row stay at least `sibling_gap` apart (`couple_gap` between partners), and the children of one union occupy a stretch of their row that no other sibship enters.
3. **Break tangled loops** — A loop partner (both partners descend from people in the pedigree) is drawn next to their partner, with their parents' sibship line stretched to reach them. Step 2 orders siblings so the two branches of a loop meet: each ancestor on one side is placed last among their siblings, each on the other side first. Where the stretched line would still overlap another sibship, the partner is drawn twice instead: among their siblings, and as a duplicate next to their partner, joined by a dotted arc. One loop is broken per group of overlapping lines, then steps 1-2 are repeated, so no loop is broken that another break made unnecessary. The input family is never modified; the duplicates appear as extra nodes with `is_duplicate=True`.
4. **Build geometry** — Nodes, edges, twin connectors, duplicate arcs. Every edge is horizontal or vertical: a single child that is not centred under its parents (possible only in a loop) hangs from a short sibship-level segment rather than a diagonal. Sibship lines of different unions never overlap. Duplicate arcs bow by at most 0.4 generation gaps.

---

### `compute_layout_circular(family, style=None)`

Compute a circular/radial layout. Founders sit at the centre and each later generation on its own concentric ring, in the left-to-right order of the hierarchical layout wrapped around the circle (counter-clockwise, centred at the bottom). It therefore inherits that layout's guarantees: siblings stay together, sibship arcs of different couples never overlap, partners sit side by side, and tangled loops are drawn with duplicated individuals. A small family occupies only the arc it needs; a large one the whole circle.

```python
from pypedigree import compute_layout_circular, PedigreeStyle

layout = compute_layout_circular(family, style=PedigreeStyle())
```

**Args**: `family` — `Family` object. `style` — `PedigreeStyle` (optional; `base_radius` and `radial_gap` control ring sizing).

`base_radius` is expanded automatically when a ring would be too crowded for
the shapes to fit, so nodes never overlap. Shapes stay upright (a square turned
to face outward would read as the unknown-sex diamond), so phenotype bands are
horizontal as in the hierarchical layout. ID labels are set radially outside
each shape, away from the lines arriving from the centre, on a white backing.
Twins hang from one point of the sibship arc, as a V (with a bar for MZ
twins), and duplicates are joined to their original by a dotted arc.

**Returns**: `FamilyLayout` with `.nodes`, `.edges` (edges may contain `waypoints`), `.twin_connectors`, `.duplicate_arcs`.

---

## PedigreeBuilder

### `PedigreeBuilder(family_id="FAM1")`

Fluent builder for constructing pedigrees programmatically. All mutation methods return `self` for chaining.

```python
from pypedigree import PedigreeBuilder

builder = PedigreeBuilder("FAM1")
builder.add_individual("father", sex="male")
builder.add_individual("mother", sex="female")
builder.add_individual("child", sex="male", father="father", mother="mother",
                       phenotype="affected", phenotypes={"IQ": 85})
fig = builder.render(output="family.png")
```

### `add_individual(sample_id, ...)`

Add an individual to the pedigree.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sample_id` | `str` | required | Unique identifier |
| `sex` | `str \| Sex` | `"unknown"` | `"male"`, `"female"`, `"unknown"` (or `"m"`, `"f"`, `"1"`, `"2"`) |
| `phenotype` | `str \| Phenotype` | `"unknown"` | `"affected"`, `"unaffected"`, `"unknown"` |
| `father` | `str \| None` | `None` | Father's sample_id |
| `mother` | `str \| None` | `None` | Mother's sample_id |
| `family_id` | `str \| None` | `None` | Override default family ID |
| `deceased` | `bool` | `False` | Deceased marker |
| `proband` | `bool` | `False` | Proband marker |
| `carrier` | `bool` | `False` | Carrier: patterned section (NSGC 2022 notation) or half-filled shape (legacy) |
| `obligate_carrier` | `bool` | `False` | Obligate carrier: sparser pattern (NSGC 2022) or center dot (legacy) |
| `date_of_birth` | `str \| None` | `None` | Date of birth |
| `phenotypes` | `dict \| None` | `None` | Extra phenotypes as `{name: value}` |
| `annotations` | `list[str] \| None` | `None` | Annotation lines displayed below the ID label (italic grey text) |
| `pregnancy` | `str \| Pregnancy \| None` | `None` | `"P"` ongoing pregnancy, `"SAB"` miscarriage, `"TOP"` termination, `"ECT"` ectopic, `"SB"` stillbirth |
| `adopted` | `str \| Adoption \| None` | `None` | `"in"` (`father`/`mother` are the adoptive parents; dashed line) or `"out"` (they are the biological parents); brackets around the symbol |
| `childless` | `str \| Childlessness \| None` | `None` | `"choice"` (one bar) or `"infertility"` (two bars) below the symbol |
| `art` | `str \| ReproductiveRole \| None` | `None` | `"D"` gamete donor or `"S"` surrogate, written inside the symbol |
| `gender` | `str \| Gender \| None` | `None` | Gender identity `"man"`, `"woman"`, `"nonbinary"`: the symbol follows it and `sex` is the sex assigned at birth (AMAB/AFAB below when they differ) |
| `consultand` | `bool` | `False` | Arrow without "P" |

An unknown notation code raises `PedigreeError`.

### `add_union(partner1, partner2, ...)`

Explicitly add a union. Only needed for childless couples or consanguineous/broken flags.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `partner1` | `str` | required | First partner's sample_id |
| `partner2` | `str` | required | Second partner's sample_id |
| `consanguineous` | `bool` | `False` | Double-line (shared ancestors) |
| `broken` | `bool` | `False` | Broken-line (dissolved) |

### `set_style(style)` / `build()` / `render(output=None, ...)` / `render_interactive_svg(output=None, ...)`

- `set_style(style)` — Set custom `PedigreeStyle`.
- `build()` — Build a `Pedigree` object. Auto-infers unions and phantom parents.
- `render(output, *, style, show_labels, show_legend, circular=False)` — Build and render as matplotlib Figure. Auto-assigns phenotype colors. Pass `circular=True` for radial layout.
- `render_interactive_svg(output, *, style, show_labels, show_legend)` — Build and render as interactive SVG string with tooltips and click highlighting.

### Jupyter integration

The builder supports inline display in Jupyter notebooks:

```python
builder  # displays the pedigree inline
```

---

## Rendering

### `render_pedigree(pedigree, style=None, output=None, ax=None, show_labels=True, show_legend=True, circular=False)`

Render all families in a pedigree to a matplotlib figure.

```python
from pypedigree import parse_pedigree, render_pedigree

ped = parse_pedigree("family.ped")
fig = render_pedigree(ped, output="family.png")
```

**Args**:
- `pedigree` — `Pedigree` object.
- `style` — `PedigreeStyle` (optional).
- `output` — File path to save (optional). Format inferred from extension.
- `ax` — Existing matplotlib `Axes` to render onto (optional).
- `show_labels` — Whether to draw sample ID labels (default: True).
- `show_legend` — Whether to draw phenotype legend (default: True).
- `circular` — Use circular/radial layout instead of the default top-down layout (default: False).

**Returns**: matplotlib `Figure`.

**Note**: For large pedigrees, figure dimensions are auto-scaled based on layout bounds.

---

### `render_family(family, style=None, ax=None, show_labels=True)`

Render a single family, returning both the figure and the computed layout.

```python
from pypedigree import parse_pedigree, render_family

ped = parse_pedigree("family.ped")
fam = ped.get_family("FAM1")
fig, layout = render_family(fam)
```

**Returns**: Tuple of (`Figure`, `FamilyLayout`).

---

### `render_interactive_svg(pedigree, style=None, output=None, show_labels=True, show_legend=True)`

Render all families as an interactive SVG with embedded JavaScript for tooltips, click highlighting, zoom, and pan. The SVG is written directly (matplotlib is only used to read color names), and draws the same figure as `render_pedigree()`: both renderers take their decisions from `pypedigree.phenotypes` and `pypedigree.carriers`.

```python
from pypedigree import parse_pedigree, render_interactive_svg

ped = parse_pedigree("family.ped")
svg_str = render_interactive_svg(ped, output="family.svg")
```

**Args**:
- `pedigree` — `Pedigree` object.
- `style` — `PedigreeStyle` (optional).
- `output` — File path to save (optional).
- `show_labels` — Whether to draw sample ID labels (default: True).
- `show_legend` — Whether to draw the phenotype and carrier legend (default: True).

**Returns**: SVG content as a string.

**Features**:
- Hover tooltips (sample_id, sex, phenotype, DOB, extra phenotypes, carrier status)
- Click node to highlight parents, children, siblings (click again to reset)
- Mouse wheel zoom, Ctrl+drag to pan
- Carrier notation: patterned sections defined once as SVG `<pattern>`s (NSGC 2022), or half-fill and center dot with `carrier_notation="legacy"`
- Phenotype bands (binary, gradient, diverging, categorical), `o_` outline colors, PGS distribution curves and `a_` annotations, from the same `style.phenotype_colors`, `gradient_range`, `gradient_diverging` and `pgs_colors` as the PNG
- The legend below the pedigree, with the same rows and proportions as the PNG legend
- Several families side by side, each titled

The circular layout is not available in the interactive SVG. Text is placed
with explicit offsets rather than `dominant-baseline`, so the file renders the
same in browsers, Inkscape and librsvg.
---

## Data Model

### Enums

| Enum | Values | Description |
|------|--------|-------------|
| `Sex` | `MALE (1)`, `FEMALE (2)`, `UNKNOWN (0)` | Biological sex |
| `Phenotype` | `UNAFFECTED (1)`, `AFFECTED (2)`, `UNKNOWN (0)` | Phenotypic status |
| `TwinType` | `MONOZYGOTIC`, `DIZYGOTIC`, `UNKNOWN` | Twin zygosity |
| `EdgeType` | `RELATIONSHIP`, `DESCENT`, `SIBSHIP`, `INDIVIDUAL`, `CHILDLESS` | Connection line type (`CHILDLESS`: the stem and bars of a childlessness marker) |
| `Pregnancy` | `ONGOING ("P")`, `MISCARRIAGE ("SAB")`, `TERMINATION ("TOP")`, `ECTOPIC ("ECT")`, `STILLBIRTH ("SB")` | A pregnancy or how it ended |
| `Adoption` | `IN`, `OUT` | Adopted in (parents given are adoptive) or out (they are biological) |
| `Childlessness` | `BY_CHOICE`, `INFERTILITY` | No children by choice or reason unknown, or infertility |
| `ReproductiveRole` | `DONOR ("D")`, `SURROGATE ("S")` | Role in assisted reproduction |
| `Gender` | `MAN`, `WOMAN`, `NONBINARY` | Gender identity; decides the symbol |
| `LineStyle` | `SOLID`, `DOUBLE`, `BROKEN`, `DASHED` | Line rendering style |

### `Individual`

Represents one person in the pedigree.

| Field | Type | Description |
|-------|------|-------------|
| `sample_id` | `str` | Unique identifier within family |
| `family_id` | `str` | Family this individual belongs to |
| `sex` | `Sex` | Biological sex |
| `phenotype` | `Phenotype` | Affected/unaffected status |
| `father_id` | `str \| None` | Father's sample_id |
| `mother_id` | `str \| None` | Mother's sample_id |
| `date_of_birth` | `str \| None` | DOB string (YYYYMMDD/YYYYMM/YYYY) |
| `twin_info` | `list[TwinInfo]` | Twin relationship declarations |
| `is_deceased` | `bool` | Whether the individual is deceased |
| `is_proband` | `bool` | Whether this is the proband |
| `is_phantom` | `bool` | Whether this is a generated phantom parent |
| `is_duplicate` | `bool` | Whether this is a duplicate node |
| `is_carrier` | `bool` | Whether this is a carrier (`carrier` column = 1) |
| `is_obligate_carrier` | `bool` | Whether this is an obligate carrier (`carrier` column = 2) |
| `duplicate_of` | `str \| None` | Original sample_id if duplicate |
| `extra_phenotypes` | `dict[str, float \| None]` | Normalized phenotype values [0, 1] |
| `extra_phenotypes_raw` | `dict[str, float \| None]` | Original raw phenotype values |
| `annotations` | `dict[str, str]` | Annotation key-value pairs (displayed as italic grey text below ID label) |
| `notation` | `Notation` | Pregnancy, adoption, childlessness, assisted reproduction, gender identity and consultand (see below) |

### `Notation`

Frozen dataclass of the standard notation beyond sex, affection status and death
(NSGC nomenclature, Bennett et al. 2008 and 2022). Every field defaults to "nothing
to note". How each is drawn is decided in `pypedigree.symbols`.

| Field | Type | Drawn as |
|-------|------|----------|
| `pregnancy` | `Pregnancy \| None` | "P" inside the symbol; a triangle for SAB, ECT ("ECT" below) and TOP (struck through); the symbol struck through with "SB" below for a stillbirth |
| `adoption` | `Adoption \| None` | Brackets; the line to adoptive parents is dashed |
| `childlessness` | `Childlessness \| None` | A line below the symbol ending in one bar, or two |
| `reproductive_role` | `ReproductiveRole \| None` | "D" or "S" inside the symbol |
| `gender` | `Gender \| None` | Symbol shape; AMAB/AFAB below when it differs from `Individual.sex` |
| `is_consultand` | `bool` | Arrow without "P" |

### `Union`

Represents a mating pair and their children.

| Field | Type | Description |
|-------|------|-------------|
| `partner_ids` | `tuple[str, str]` | (male_id, female_id) |
| `children_ids` | `list[str]` | Children of this union |
| `is_consanguinous` | `bool` | Auto-detected shared ancestry |
| `is_broken` | `bool` | Separated/divorced (API-only) |

### `Family`

Collection of individuals and unions sharing a family_id.

| Method | Returns | Description |
|--------|---------|-------------|
| `get_founders()` | `list[Individual]` | Individuals with no parents |
| `get_children(sample_id)` | `list[Individual]` | Children of a given parent |
| `get_union_for_parents(id1, id2)` | `Union \| None` | Find union by partner IDs |

### `Pedigree`

Top-level container for multiple families.

| Property / Method | Returns | Description |
|-------------------|---------|-------------|
| `family_ids` | `list[str]` | List of family identifiers |
| `get_family(family_id)` | `Family \| None` | Look up a family by ID |
| `phenotype_names` | `list[str]` | Names of extra phenotype columns |
| `annotation_names` | `list[str]` | Names of annotation columns (from `a_` prefixed headers) |

### Layout Types

| Type | Key Fields | Description |
|------|------------|-------------|
| `NodeLayout` | `individual_id, x, y, sex, phenotype, is_deceased, is_proband, is_phantom, is_duplicate, is_carrier, is_obligate_carrier, label, extra_phenotypes, extra_phenotypes_raw, annotations, pgs_side, notation, notes, marker_depth` | Position and attributes of one node. `sex` is the shape drawn (from the gender identity when given); `is_deceased` includes stillbirths. `annotations: list[str]` contains annotation lines to display below the ID label, after `notes` ("SB", "ECT", "AMAB", "AFAB"). `marker_depth` is the room a childlessness marker takes below the symbol; the label moves down by it. `pgs_side: str` indicates PGS curve placement ("left" or "right"). |
| `EdgeLayout` | `start, end, edge_type, line_style, waypoints` | A connection line segment (waypoints for polyline/curved edges) |
| `TwinConnector` | `convergence_point, child_positions, twin_type` | V-lines for twin groups |
| `DuplicateArc` | `from_pos, to_pos` | Dotted arc between duplicate nodes |
| `FamilyLayout` | `nodes, edges, twin_connectors, duplicate_arcs` | Complete layout for one family |

---

## Configuration

### `PedigreeStyle`

All visual parameters in one dataclass. Pass to `compute_layout()` or `render_pedigree()`.

```python
from pypedigree import PedigreeStyle

style = PedigreeStyle(shape_size=1.0, dpi=300)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `shape_size` | `0.8` | Size of individual shapes |
| `shape_line_width` | `1.5` | Shape border width |
| `fill_color` | `"black"` | Fill color for affected individuals |
| `outline_color` | `"black"` | Shape outline and deceased marker color |
| `phantom_line_style` | `"--"` | Line style for phantom parents |
| `line_color` | `"black"` | Connection line color |
| `line_width` | `1.0` | Connection line width |
| `double_line_gap` | `0.06` | Gap between double lines (consanguinity) |
| `broken_gap_size` | `0.15` | Gap size in broken lines |
| `carrier_notation` | `"nsgc2022"` | `"nsgc2022"`: the symbol is divided into patterned sections, one per carrier result, defined in the legend (Bennett et al., J Genet Couns 2022). `"legacy"`: half-filled shape for carriers, center dot for obligate carriers. Any other value raises `PedigreeError` |
| `carrier_dot_radius` | `0.08` | Radius of the obligate carrier center dot (legacy notation) |
| `font_size` | `8.0` | Label font size, in points. The legend panel is sized from it (text, swatches and row height), so raise it to keep a figure legible once scaled down for print |
| `font_family` | `"sans-serif"` | Label font family |
| `label_offset` | `0.15` | Vertical offset for labels below shapes |
| `generation_gap` | `2.0` | Vertical distance between generations |
| `sibling_gap` | `1.5` | Horizontal distance between siblings |
| `couple_gap` | `1.5` | Horizontal distance between partners |
| `min_node_distance` | `1.2` | Minimum distance between any two nodes |
| `radial_gap` | `2.5` | Distance between concentric rings in circular layout |
| `base_radius` | `3.0` | Radius of the innermost ring in circular layout |
| `phenotype_colors` | `{}` | Map of phenotype name to hex color (e.g. `{"IQ": "#3498db"}`) |
| `show_phenotype_labels` | `False` | Show YES/NO text labels inside binary phenotype bands |
| `pgs_colors` | `{}` | Map of PGS column name to hex color for the marker line (default marker color: `#e74c3c`) |
| `gradient_range` | `{}` | Custom min/max for gradient normalization, e.g. `{"SRS": [40, 90]}` |
| `gradient_diverging` | `{}` | Tricolor diverging gradient config, e.g. `{"IQ": {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}}` |
| `zero_pad_ids` | `0` | Zero-pad sample IDs to this width in labels (0 = no padding) |
| `dpi` | `150` | Output resolution |
| `fig_width` | `10.0` | Figure width in inches |
| `fig_height` | `8.0` | Figure height in inches |
| `margin` | `0.5` | Margin around the diagram |

**Methods**: `PedigreeStyle.from_dict(d)` — create from dictionary. `style.to_dict()` — convert to dictionary.

### `DEFAULT_PHENOTYPE_PALETTE`

List of 8 default hex colors used to auto-assign phenotype band colors when no `phenotype_colors` are specified: red, blue, orange, green, purple, teal, yellow, pink.

### `PhenotypeDisplayMode`

Enum controlling how extra phenotype columns are rendered. Set via column name prefix in the .ped header:

| Prefix | Mode | Rendering |
|--------|------|-----------|
| `b_` | `BINARY` | Color (affected) / white (unaffected) / grey (unknown) / NA (no data) bands. When all bands are NA, a single centered "NA" is shown. |
| `g_` | `GRADIENT` | White to color continuous gradient with value labels |
| `gi_` | `GRADIENT_INVERTED` | Color to white gradient (high value = white, low value = full color). Useful for metrics where low = bad (e.g., IQ: low IQ = more red). |
| `o_` | `OUTLINE` | Colored shape border (not a band) |
| `c_` | `CATEGORICAL` | Distinct color per unique string value |
| `cr_` | `CARRIER` | Carrier status of one gene or condition (0 = not, 1 = carrier, 2 = obligate carrier), drawn as a patterned section, not a phenotype band |
| *(none)* | `AUTO` | Auto-detect from data (default, backward compatible) |

**Example header:**
```
FID  ID  Father  Mother  Sex  b_ASD  g_IQ  o_Proband  c_Diagnosis
```

The prefix is stripped from the display name (e.g. `b_ASD` is shown as `ASD`).

**Multiple PGS columns**: Any column whose name starts with `PGS` (case-insensitive) is treated as a PGS column -- e.g. `PGS_autism`, `PGS_IQ`, `PGS_ADHD`. All PGS values are rendered on a single gaussian curve with multiple colored markers. Marker colors are set via the `pgs_colors` style option. When PGS columns are present, `couple_gap` and `sibling_gap` are auto-increased to accommodate the curve.

**Gradient value display**: Gradient band values (e.g. "IQ: 90", "SRS: 76") are always displayed regardless of band color intensity.

**AUTO mode resolution**: A column with no prefix is resolved **once per column**, from all of its values across the pedigree:

| Values found in the column | Effective mode |
|---------------------------|----------------|
| Any string value | categorical |
| Any numeric value outside `0`, `1`, `2` | gradient |
| Otherwise (only `0`/`1`/`2`, or empty) | binary |

The decision is never taken per individual, so a column cannot mix binary and
gradient bands within the same diagram, and its legend row always matches its
bands. An explicit `gradient_range` or `gradient_diverging` entry for the
column also forces gradient.

---

### `resolve_style(style, phenotype_names=None)`

Fill in the pedigree-dependent parts of a style: phenotype columns with no
explicit color get one from `DEFAULT_PHENOTYPE_PALETTE`. Bands are only drawn
when `phenotype_colors` is populated, so this is what makes a multi-phenotype
figure render. `render_pedigree()`, the CLI and `PedigreeBuilder` all call it,
which is why they produce identical figures.

```python
from pypedigree import parse_pedigree, resolve_style

ped = parse_pedigree("multi_pheno.ped")
style = resolve_style(None, ped.phenotype_names)
print(style.phenotype_colors)   # {'epilepsy': '#e74c3c', 'IQ': '#3498db', ...}
```

**Args**: `style` — `PedigreeStyle` or None for defaults. `phenotype_names` — ordered phenotype column names (e.g. `pedigree.phenotype_names`).

**Returns**: a new `PedigreeStyle`. The input is never mutated, and the
function is idempotent.

---

### `apply_pgs_spacing(style, phenotype_names=None)`

Widen `couple_gap` and `sibling_gap` to leave room for PGS distribution
curves drawn beside each node. Returns `style` unchanged when no PGS column
is present.

Unlike `resolve_style()` this is **not** idempotent — each call widens
further. `render_pedigree()` is its only caller; do not apply it yourself
before calling `render_pedigree()`.

**Args**: `style` — `PedigreeStyle`. `phenotype_names` — ordered phenotype column names.

**Returns**: a new `PedigreeStyle`, or the input unchanged.

---

### `is_pgs_column(name)`

Whether a phenotype column name denotes a polygenic score, i.e. starts with
`PGS` (case-insensitive).

**Returns**: `bool`.

---

## Exceptions

| Exception | Parent | When raised |
|-----------|--------|-------------|
| `PedigreeError` | `Exception` | Base class for all pypedigree errors; also raised for an invalid color in a style |
| `PedigreeParseError` | `PedigreeError` | File format or parsing issues |
| `PedigreeValidationError` | `PedigreeError` | Data validation failures (cycles, invalid references) |
| `PedigreeLayoutError` | `PedigreeError` | Layout computation failures |

Every color in a style is validated before drawing starts, and all offending
values are reported in one error:

```
PedigreeError: Invalid color(s) in style: outline_color = '#12345';
phenotype_colors['Epilepsy'] = '#zzzzzz'. Expected hex values like
'#e74c3c' or named colors like 'grey'.
```

## Warnings

Non-fatal diagnostics use dedicated categories, so they can be filtered or
escalated independently of other libraries' warnings.

| Warning | Parent | When emitted |
|---------|--------|--------------|
| `PedigreeWarning` | `UserWarning` | Base class for all pypedigree warnings |
| `PedigreeStyleWarning` | `PedigreeWarning` | A style definition contains something unusable, e.g. an unknown option |
| `PedigreeDataWarning` | `PedigreeWarning` | Pedigree data is questionable but still renderable, e.g. a father recorded as female |

```python
import warnings
from pypedigree import PedigreeDataWarning, parse_pedigree

# Treat questionable data as a hard failure
warnings.simplefilter("error", PedigreeDataWarning)
parse_pedigree("family.ped")

# Or silence only pypedigree's data warnings
warnings.simplefilter("ignore", PedigreeDataWarning)
```

An unknown style option is ignored but reported, with a suggestion when one is
close:

```python
PedigreeStyle.from_dict({"phenotype_color": {...}})
# PedigreeStyleWarning: Unknown style option 'phenotype_color' ignored;
# did you mean 'phenotype_colors'?
```

Pass `from_dict(d, warn_unknown=False)` to accept a dictionary that
legitimately carries extra data.

The CLI captures these and prints them as `Warning: <message>` on stderr.
