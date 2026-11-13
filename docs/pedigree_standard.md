# Pedigree Standard

## File format

Pedigree are usually stored as text file, tab-separated value id the most common, but space, multiple-spaces or commas are also possible separators. There is one row per individual.

They are made of 6 columns:

- the **family_id** or **fid**: the family unique identifier
- the **sample_id** or **iid**: an individual sample identifier
- the **paternal_id** or **father**: the sample identifier of the father
- the **maternal_id** or **mother**: the sample identifier of the mother
- the **sex**: the biological sex of the individual designed by the *sample_id* (1 = male, 2 = female, 0 or -9 = unknown)
- the **phenotype**: a phenotypic value that can be binary (1 = unaffected, 2 = affected, 0 or -9 = unknown)

For robustness we will treat `0` or `nan` in the paternal_id and/or maternal_id column as no information available.

### Extended format

Additional columns can be appended to the standard 6-column format. All extended columns are optional — the parser auto-detects how many columns are present.

- 7th column: **date_of_birth** (dob) using a format that could be YYYYMMDD or YYYYMM or YYYY but should be uniform across the file. Use `0` or `-9` for unknown.
- 8th column: **twins** — a comma-separated list of values describing twin relationships:
  - `MZ:sample_id` to indicate that `sample_id` is a monozygotic twin to the current individual
  - `DZ:sample_id` to indicate that `sample_id` is a dizygotic twin to the current individual
  - `TW:sample_id` to indicate that `sample_id` is a twin of unknown zygosity to the current individual
  - Use `0` or `-9` for no twin information.
- 9th column: **deceased** — whether the individual is deceased (1 = deceased, 0 = alive). Use `0` or `-9` for unknown/alive.
- 10th column: **proband** — whether the individual is the proband, i.e. the family member at the origin of the medical consultation (1 = proband, 0 = not proband). Use `0` or `-9` for not proband.

### Inline phenotypes

When a **header row is present**, the parser uses **smart column mapping** to detect phenotype columns:

- Column names are matched against known structural names (FID, ID, Father, Mother, Sex, Phenotype, DOB, Twins, Deceased, Proband, Carrier) with multiple aliases supported.
- Any unrecognised column name after the core 5 (FID, ID, Father, Mother, Sex) is treated as an **extra phenotype column**.
- If no explicit `Phenotype` column is present, the first extra phenotype column is used as the primary affected/unaffected status.
- This means **you don't need the full 10-column format** — phenotype columns can appear directly after Sex.

Phenotype types:

- **Binary phenotypes**: values use standard pedigree encoding — `2`=affected (color band), `1`=unaffected (white band), `0`=unknown (grey band), `-9`=no data (white band with "NA" text). When all bands for an individual are NA (`-9`), a single centered "NA" is shown instead of repeating per band.
- **Continuous phenotypes**: numeric values > 2 — rendered as color gradients (min→white, max→full color).
- **Missing values**: empty, `.`, `NA`, `NaN` — treated as `-9` (no data).
- **PGS columns**: any column whose name starts with `PGS` (case-insensitive) is treated as a PGS column -- e.g. `PGS`, `PGS_autism`, `PGS_IQ`. Multiple PGS values are displayed on a single gaussian curve with distinct colored markers. Marker colors are controlled by the `pgs_colors` style option. When PGS columns are present, `couple_gap` and `sibling_gap` auto-increase to accommodate the curve.
- **Annotation columns**: columns with the `a_` prefix (e.g. `a_Role`, `a_Clinic`) are annotation columns. Their values are displayed as italic grey text below the individual's ID label. Missing values (0, -9, NA, empty) produce no annotation line. Annotations can also be loaded from a separate file via the `--annotations` CLI option.
- **Comma decimal separator**: the parser accepts comma as decimal separator for European locales (e.g. `0,712` is parsed as `0.712`).

### Column name prefixes

Column name prefixes control how phenotype columns are rendered:

| Prefix | Mode | Rendering |
|--------|------|-----------|
| `b_` | Binary | Color (affected) / white (unaffected) / grey (unknown) / NA (no data) bands |
| `g_` | Gradient | White to color continuous gradient with value labels |
| `gi_` | Gradient Inverted | Color to white gradient (high value = white, low value = full color). Useful for metrics where low = bad (e.g., IQ). |
| `o_` | Outline | Colored shape border (not a band) |
| `c_` | Categorical | Distinct color per unique string value |
| `cr_` | Carrier | Carrier status of one gene or condition: `0` not a carrier, `1` carrier, `2` obligate carrier, `-9`/`NA` unknown. Drawn as a patterned section of the symbol |
| `a_` | Annotation | Text displayed below the ID label (italic grey), not a phenotype band |
| *(none)* | Auto | Auto-detect from data (default) |

The prefix is stripped from the display name (e.g. `b_ASD` is shown as `ASD`).

### Style JSON options for gradients

The style JSON supports additional gradient configuration:

- **`gradient_range`**: custom min/max for gradient normalization. Example: `{"SRS": [40, 90]}` normalizes SRS values between 40 and 90.
- **`gradient_diverging`**: tricolor diverging gradient. Values below the midpoint blend from `color_low` to white; values above blend from white to `color_high`. Example: `{"IQ": {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}}`.

Gradient band values (e.g. "IQ: 90") are always displayed regardless of band color intensity.

Example — compact format (phenotypes directly after Sex):

```
FID	ID	Father	Mother	Sex	ASD	Synesthesia
FAM1	001	002	003	1	2	2
FAM1	002	0	0	1	1	2
```

Example — full 10-column format with extras:

```
FID	IID	PAT	MAT	SEX	PHENO	DOB	TWINS	DECEASED	PROBAND	epilepsy	IQ	SRS	PGS
mockup	1	0	0	1	2	0	0	0	0	1	105	62	0.72
mockup	2	0	0	2	1	0	0	0	0	0	118	45	0.35
```

### Separate phenotype file (.pheno)

Phenotype data can also be provided as a separate TSV file via the `--pheno` CLI option. Format: first column is `sample_id`, subsequent columns are phenotype values. A header row is required.

```
sample_id	epilepsy	IQ	SRS
1	1	105	62
2	0	118	45
```

### Standard notation columns

With a header row, these columns are recognised by name (case-insensitive).
Their codes are case-insensitive; `0`, `-9`, `NA`, `.` or an empty cell mean
nothing to note, and any other value is an error naming the line.

- **pregnancy** (aliases `pregnancy_outcome`, `outcome`): `P` ongoing pregnancy,
  `SAB` spontaneous abortion (miscarriage), `TOP` termination of pregnancy, `ECT`
  ectopic pregnancy, `SB` stillbirth. Long forms are accepted (`miscarriage`,
  `termination`, `ectopic`, `stillbirth`).
- **adopted** (alias `adoption`): `in` (or `1`) when the parents given are the
  adoptive parents, `out` (or `2`) when they are the biological parents of a
  child adopted out.
- **childless** (aliases `childlessness`, `infertility`): `choice` (or `1`) for
  no children by choice or for an unknown reason, `infertility` (or `2`).
- **art** (alias `assisted_reproduction`): `D` gamete (egg or sperm) donor, `S`
  surrogate (`GC`, gestational carrier, is accepted too).
- **gender_identity** (alias `identity`): `man` (`M`, `1`), `woman` (`F`, `W`,
  `2`), `nonbinary` (`NB`, `X`, `3`). When it is given, the `sex` column is read
  as the sex assigned at birth. A column named `gender` is the gender identity
  only when the header also has a `sex` column; alone it is the sex column.
- **consultand**: `1` for the consultand (the person seeking counselling, when
  that is not the proband).

Gestational age, karyotype, the reason for infertility and similar details go in
`a_` annotation columns.

## Diagram representation

### 1. Individual

**sex**: A male individual is coded in the diagram by a square, a female by a circle, and an individual with unspecified sex is a diamond.

:::image type="content" source="imgs/individual_unaffected_male.png" alt-text="male individual":::

:::image type="content" source="imgs/individual_unaffected_female.png" alt-text="female individual":::

:::image type="content" source="imgs/individual_unaffected_unspecified_sex.png" alt-text="unspecified sex individual":::

**phenotype**: unaffected individuals are symbolized by empty shapes while affected one are symbolized by filled shapes.

:::image type="content" source="imgs/individual_affected_male.png" alt-text="affected male":::

:::image type="content" source="imgs/individual_unaffected_male.png" alt-text="unaffected male":::

**deceased**: deceased people are marked with a diagonal line starting at the bottom left of the shape and going to the top right of it. The line should stick out a little bit on both sides of the shape.

:::image type="content" source="imgs/deceased_female.png" alt-text="deceased female":::

:::image type="content" source="imgs/deceased_male.png" alt-text="deceased male":::

:::image type="content" source="imgs/deceased_unspecified_sex.png" alt-text="deceased person of unspecified sex":::

**proband**: the proband (the family member at the origin of the medical consultations) is marker by a small diagonal arrow pointing up-right at the bottom left corner of the shape.

:::image type="content" source="imgs/male_proband.png" alt-text="male proband":::

:::image type="content" source="imgs/female_proband.png" alt-text="female proband":::

### 2. Relationship line

A horizontal line between a male individual and a female individual symbolize a relationship. That is the case if both individuals appears together as parents of a third individual.
If the relationship is consanguinous, then it should be a double horizontal line.
Male is expected to be on the left side of a relationship line.
this line is located verticaly in the middle of the "individual shapes".
A past relationship is mark by a "break" on the relationship line

:::image type="content" source="imgs/diagram_informations.png" alt-text="diagram informations":::

:::image type="content" source="imgs/consanguinous_relationship.png" alt-text="consanguinous relationship":::

:::image type="content" source="imgs/broken_relationship.png" alt-text="broken relationship":::

#### Consanguinity detection

Consanguinity (related partners) is automatically detected during parsing.
A union is marked as consanguinous if both partners share at least one
common ancestor. The detection algorithm performs ancestor traversal
(BFS upward through parent links) for each partner in every union,
then checks if the two ancestor sets intersect.

This covers all degrees of consanguinity:
- Sibling marriages (shared parents)
- First-cousin marriages (shared grandparents)
- Uncle-niece/aunt-nephew marriages (one partner is ancestor of the other's parent)
- Second-cousin and more distant relationships (shared great-grandparents, etc.)

No additional file column is needed — consanguinity is inferred from the
existing parent-child structure.

#### Broken relationships

Broken (past/dissolved) relationships are currently API-only. The
`Union.is_broken` field can be set programmatically for rendering
a relationship line with a gap and "//" marks. No file format column
is currently defined for broken relationships.

### 3. Line of Descent

a line of descent start at the midpoint of a relationship line and go down to a point below the parents.

:::image type="content" source="imgs/diagram_informations.png" alt-text="diagram informations":::

### 4. Sibship line

when there is more than one birth coming from the union of 2 individuals, there is a sibshipline. this line is horizontal, parallel to the relationship line and located at the bottom of the line of descent. this line is of variable length depending on the number of children and shoule be center on the bottom of the descent line.

:::image type="content" source="imgs/diagram_informations.png" alt-text="diagram informations":::

siblings should be listed from left to rigth from the oldest to the younger

### 5. Individual line

Each children symbolized as an individual, is linked by a small vertical line at the top to a sibship line. if there is a single kid in the relationship then the small vertical line connect to the bottom of the descent line.

When a single child cannot be placed directly below its parents (because neighbouring individuals take the space), a short horizontal segment at sibship level joins the bottom of the line of descent to the child's vertical line. Lines between individuals are always horizontal or vertical, never diagonal. The same segment joins a twin group that is the only offspring of a union to the line of descent.

:::image type="content" source="imgs/diagram_informations.png" alt-text="diagram informations":::

**multiple gestation**: twins are symbolized side by side by 2 shapes connected from their top middle point to the

:::image type="content" source="imgs/monozygotic_twins.png" alt-text="monozygotic twins":::

:::image type="content" source="imgs/dizygotic_twins.png" alt-text="dizygotic twins":::

:::image type="content" source="imgs/unknown_zygocity_twins.png" alt-text="unknown zygocity twins":::

:::image type="content" source="imgs/trizygotic_twins.png" alt-text="trizygotic":::

### 6. Complementary informations

If only a parent is indicated in the pedigree file (so nothing, 0 or -9 for the other parent), then a "fake" parent made of a dashed border can be created in its place.

If one individual is involved in multiple relationship, we should plot complete tree side by side for each relationship and connect the multiple isntance of the same individual through dotted arcs.

If incestuous consiguinous relationship crossing generations are present, the older individual (higher up in the tree) can be duplicated to be put in a relationship as a "new" individual but with a dotted arc to link it to its duplicate.

When both partners of a relationship descend from people in the pedigree (a loop: consanguinity, or a marriage between two branches), they are drawn side by side and the sibship line of the second partner's parents stretches to reach them. Siblings are ordered so that the two branches of a loop meet (the ancestors on one side placed last among their siblings, those on the other side first), which keeps that line clear of other sibships in most pedigrees. If it would still overlap the sibship line of another couple, the second partner is drawn twice instead: once among their siblings, and once as a duplicate next to their partner, linked by a dotted arc. Such reordering takes precedence over birth order. Sibship lines of two different couples must never overlap, since a child hanging from merged lines cannot be traced to its parents.

The sample_id should be written under each individual.

#### Multiple relationships (duplicate nodes)

When an individual participates in more than one union (has children with
multiple partners), the pedigree is drawn with separate subtrees for each
relationship. The individual appears once as the "original" in the first
union and as a "duplicate" node in each subsequent union. Duplicates are
connected to the original by a dotted arc.

This is handled automatically during parsing:

1. The parser detects individuals who appear in more than one inferred union.
2. For each additional union, a duplicate node is created with:
   - The same shape, fill (sex, phenotype), and deceased marker as the original
   - A label showing the original sample_id (so both instances show the same name)
   - The `is_duplicate` flag set and `duplicate_of` pointing to the original
3. The union's partner_ids and children's parent references are updated to
   point to the duplicate instead of the original.
4. The layout engine positions each subtree independently and computes
   `DuplicateArc` geometry between the original and duplicate positions.
5. The renderer draws a dotted curved arc connecting the two instances.

No additional file column is needed — multiple relationships are inferred
from the existing parent-child structure (each unique parent pair forms a
separate union).

#### Carrier notation

Carrier status is drawn following the NSGC focused revision of the nomenclature
(Bennett et al., J Genet Couns 2022;31(6):1238-48), which no longer uses a dot:
the symbol is divided into as many vertical sections as needed, each with its
own fill pattern explained in the legend. An affected section, filled, comes
first, then one patterned section per carrier result. A carrier of a single
result therefore has a fully patterned symbol; an affected carrier, or a carrier
of two conditions, a symbol split in two.

- The `carrier` column (0/1/2) and each `cr_` column are separate carrier
  sources, each with its own pattern. Obligate carriers (value `2`) get the same
  pattern at half the density.
- When phenotype bands are drawn, each carrier source is an extra band below them.
- `carrier_notation: "legacy"` restores the pre-2022 notation for the `carrier`
  column: a half-filled symbol for carriers and a central dot for obligate
  carriers.

### 7. Standard notation

Following the NSGC nomenclature (Bennett et al., J Genet Couns 2008;17:424-33,
and its 2022 revision centred on sex and gender, 31(6):1238-48). Text below a
symbol is written under the ID label: notation notes first ("SB", "ECT",
"AMAB"...), then annotations.

**Gender identity**: the symbol represents gender, not the sex assigned at
birth: a square for a man, a circle for a woman, a diamond for a non-binary or
gender-diverse person. When the symbol differs from the sex assigned at birth,
AMAB (assigned male at birth) or AFAB (assigned female at birth) is written
below it. A cisgender person gets no annotation.

**Pregnancy**: a "P" inside the symbol (a diamond when the sex is unknown);
white on a filled (affected) symbol.

**Pregnancy loss**: a triangle whose apex meets the individual's line, filled
when affected. A termination (TOP) is struck through; an ectopic pregnancy has
"ECT" below. A triangle is too small for phenotype bands: in a figure with
bands it takes the color the individual's first band would have.

**Stillbirth**: the symbol of the fetus' sex (a diamond if unknown), struck
through like a deceased individual, with "SB" below.

**Adoption**: brackets around the symbol for every adoption. The line to the
adoptive parents (adopted in) is dashed; the line to the biological parents
(adopted out) is solid.

**No children by choice or for an unknown reason**: a short vertical line below
the symbol ending in a horizontal bar. **Infertility**: the same with two bars.
The label is moved below the marker. In the circular layout the marker points
outward, on the side of the label.

**Assisted reproduction**: "D" (gamete donor) or "S" (surrogate) inside the
symbol. On a symbol with phenotype bands, which carry text of their own, the
letter is written beside the upper-right corner instead. The diagonal lines that
join a donor or a surrogate to the pregnancy in the NSGC figures are not drawn:
in a pedigree file a donor is a parent, joined to the other parent like any
couple.

**Consultand**: the proband's arrow without the "P".

Not covered: several individuals in one symbol ("n"), documented evaluation
("*"), the asymptomatic-carrier vertical line, and childlessness of a couple
(drawn from the relationship line in the NSGC figures; here the marker hangs
from one individual).
