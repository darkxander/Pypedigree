# Test fixtures

**Every pedigree in this directory is synthetic.** None of it derives from a
real family, cohort or participant. Topologies are constructed by hand or by
`generate_fixtures.py`, and all phenotype values come from a seeded random
number generator.

Gene symbols appearing in `a_` annotation columns (`DMXL2`, `CAMK2B`,
`SHANK3`, `SCN2A`) are real, publicly documented symbols used as illustrative
labels on invented individuals. They carry no information about anyone.

## If you contribute a fixture

Contribute synthetic data only. Do not add pedigrees derived from real
individuals, even pseudonymised: cohort accession numbers, study identifiers
and sample barcodes are pseudonyms, not anonymisation, and a large pedigree
carrying a rare phenotype can be identifying on its own. Reduce a
reproduction case to a fabricated file that still shows the bug.

## Layout

| Directory | What it covers |
|-----------|----------------|
| `basic/` | Trios, sibships, three generations, multiple families |
| `extended/` | Deceased, proband, date-of-birth ordering, partial column sets |
| `twins/` | Monozygotic, dizygotic, unknown zygosity, triplets, one-sided |
| `consanguinity/` | First cousins, uncle–niece, sibling union |
| `nomenclature/` | Pregnancies and pregnancy losses, adoption, childlessness, gamete donor, gender identity, consultand |
| `multi_partner/` | Individuals with several partners, duplicate nodes |
| `phenotype/` | Multi-phenotype bands, prefixes, PGS, annotations, carriers |
| `plink/` | `.fam` files and PLINK `.ped` with genotype columns |
| `format/` | Separator detection, headers, unicode IDs, malformed input |
| `complex/` | Large and structurally awkward pedigrees, including the gallery's five-generation example with two loops (`multigenerational.ped`) |

Each `.ped` is paired with a `.png` rendered at the default 150 dpi with no
style file, as a visual reference. These are documentation, not assertions —
the pixel-comparison suite uses `tests/visual/golden/` instead.

## Regenerating

Fixtures written by the generator are reproducible:

```bash
uv run python tests/data/generate_fixtures.py
```

The seed is fixed, so this is a no-op unless the generator is edited.
`tests/test_fixtures.py` asserts the files on disk match what the generator
produces, so an accidental hand-edit is caught by the test suite.

Hand-written fixtures (everything outside `generate_fixtures.py`'s output) are
edited directly. After changing any `.ped`, re-render its `.png`:

```bash
uv run pypedigree tests/data/<dir>/<name>.ped -o tests/data/<dir>/<name>.png
```
