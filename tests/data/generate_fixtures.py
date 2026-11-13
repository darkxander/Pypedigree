#!/usr/bin/env python3
"""Generate the synthetic test fixtures under `complex/` and `phenotype/`.

Every pedigree this script writes is **fabricated**. Topologies are built from
the loop structures below and all phenotype values come from a seeded RNG, so
nothing here derives from a real family, cohort or participant. Gene symbols
used as annotation values are real, publicly documented symbols attached to
invented individuals; they say nothing about anyone.

The script is deterministic: the same seed always produces the same files, so
regenerating is a no-op unless the structure below is edited.

    uv run python tests/data/generate_fixtures.py

Structures are chosen to exercise specific renderer features, not to resemble
any particular study. See `tests/data/README.md`.
"""

from __future__ import annotations

import random
from pathlib import Path

DATA_DIR = Path(__file__).parent
SEED = 20260903

# Real gene symbols used as illustrative annotation labels on invented people.
GENE_LABELS = ["DMXL2", "CAMK2B", "SHANK3", "SCN2A"]


# Column sets, one per fixture. Kept byte-identical to what the fixtures had
# before, so parsing behaviour and feature coverage are unchanged.
HEADER_EXTENDED = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PHENO",
    "DOB",
    "TWINS",
    "DECEASED",
    "PROBAND",
    "epilepsy",
    "IQ",
    "SRS",
    "PGS",
]
HEADER_TWO_BINARY_PGS = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PROBAND",
    "b_ASD",
    "b_ADHD",
    "IQ",
    "SRS",
    "PGS",
]
HEADER_ALL_MODES = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PROBAND",
    "b_Autism",
    "IQ",
    "g_SRS",
    "PGS_autism",
    "PGS_int",
    "PGS_ADHD",
    "a_Variant1",
    "a_Variant2",
]
HEADER_ANNOTATIONS = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PROBAND",
    "b_TSA",
    "IQ",
    "g_SRS",
    "PGS_autism",
    "PGS_IQ",
    "PGS_ADHD",
    "a_Variant1",
    "a_Variant2",
]


def _write(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Write a tab-separated pedigree file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines += ["\t".join(str(f) for f in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    print(f"  wrote {path.relative_to(DATA_DIR.parent.parent)} ({len(rows)} individuals)")


def build_large_multi_union() -> tuple[list[str], list[list[str]]]:
    """A founder couple whose ten children each marry out and have children.

    Exercises: `b_` and `o_` prefixes, structural Deceased/Proband columns,
    `-9` no-data values, a crowded middle generation (20 individuals on one
    ring, which forces the circular layout to expand its base radius), and
    more than ten unions.

    IDs are laid out in blocks of four — child of the founders, their spouse,
    then two children — with the last three individuals added as extra
    children, so ID 1 is a child of the founders and ID 3 a grandchild.
    """
    rng = random.Random(SEED)
    header = [
        "FID",
        "ID",
        "Father",
        "Mother",
        "Sex",
        "b_ASD",
        "b_Synesthesia",
        "Deceased",
        "Proband",
        "o_Psychiatric comorbidities",
    ]
    fid = "FAM_SYN01"
    father, mother = "45", "44"

    rows: list[list[str]] = []

    def add(sid: str, pat: str, mat: str, sex: int, *, deceased: int = 0, proband: int = 0) -> None:
        """Append one individual with fabricated phenotype values."""
        # A tenth of the pedigree carries no phenotype data at all.
        if rng.random() < 0.16:
            asd, syn = "-9", "-9"
        else:
            asd, syn = str(rng.choice([1, 1, 2])), str(rng.choice([1, 1, 2]))
        rows.append([fid, sid, pat, mat, sex, asd, syn, deceased, proband, rng.choice([1, 2])])

    # Ten sibships, each a child of the founder couple plus a married-in spouse.
    next_id = 1
    couples: list[tuple[str, str]] = []
    for block in range(10):
        child = str(next_id)
        spouse = str(next_id + 1)
        grandchildren = [str(next_id + 2), str(next_id + 3)]
        next_id += 4

        child_sex = 1 if block % 2 == 0 else 2
        spouse_sex = 2 if child_sex == 1 else 1
        # ID 1 is the proband, which the integration tests rely on.
        add(child, father, mother, child_sex, proband=1 if child == "1" else 0)
        add(spouse, "0", "0", spouse_sex)

        pat, mat = (child, spouse) if child_sex == 1 else (spouse, child)
        couples.append((pat, mat))
        for i, gc in enumerate(grandchildren):
            add(gc, pat, mat, 1 if i % 2 == 0 else 2)

    # Three extra grandchildren, so the youngest generation is not uniform.
    for i, (pat, mat) in enumerate(couples[:3]):
        add(str(next_id), pat, mat, 2 if i % 2 == 0 else 1)
        next_id += 1

    # The founder couple last, so they take the highest IDs (44 and 45).
    add(mother, "0", "0", 2)
    add(father, "0", "0", 1, deceased=1)

    # ID 1 must read as affected and ID 3 as unaffected: the first phenotype
    # column drives Individual.phenotype, which the integration tests assert.
    by_id = {r[1]: r for r in rows}
    by_id["1"][5] = "2"
    by_id["3"][5] = "1"

    return header, rows


def build_multigenerational() -> tuple[list[str], list[list[str]]]:
    """Five generations descending from one founder couple, with two loops.

    The gallery's large example: about seventy individuals, a first-cousin
    marriage and a marriage between two branches (second cousins), dizygotic
    and monozygotic twins, three phenotype bands (`b_`, `g_`, `cr_`),
    deceased and proband markers, a miscarriage, an adoption and annotations.
    Identifiers are cohort-style, generation then rank in the file (203 is
    the third individual of generation 2). They are not the clinical II-3
    numbering, which runs left to right across the drawing: the layout
    reorders sibships so the branches of a loop meet, so file order and
    drawing order differ.
    """
    rng = random.Random(SEED + 101)
    people: list[dict] = []

    def new(gen: int, sex: int, father: dict | None = None, mother: dict | None = None) -> dict:
        person = {"gen": gen, "sex": sex, "father": father, "mother": mother}
        person.update(deceased=0, proband=0, twins=[], pregnancy="0", adopted="0", note="0")
        people.append(person)
        return person

    def children(father: dict, mother: dict, sexes: list[int]) -> list[dict]:
        return [new(father["gen"] + 1, s, father, mother) for s in sexes]

    def marry_out(person: dict) -> tuple[dict, dict]:
        # The spouse comes right after the person in the file, so numbering
        # in file order runs roughly left to right across the drawing.
        spouse = new(person["gen"], 3 - person["sex"])
        people.pop()  # dicts compare by value: move it by position, not with remove()
        at = next(i for i, q in enumerate(people) if q is person)
        people.insert(at + 1, spouse)
        return (person, spouse) if person["sex"] == 1 else (spouse, person)

    def twins(a: dict, b: dict, kind: str) -> None:
        a["twins"].append((kind, b))
        b["twins"].append((kind, a))

    # I: the founders.
    f1, f2 = new(1, 1), new(1, 2)
    # II: four children, each married to someone from outside the family.
    gen2 = children(f1, f2, [1, 2, 1, 2])
    couples2 = [marry_out(p) for p in gen2]
    # III: 4, 3, 4 and 3 children; the third sibship ends with DZ twins.
    gen3 = [
        children(pat, mat, sexes)
        for (pat, mat), sexes in zip(
            couples2, [[1, 2, 2, 1], [1, 2, 1], [2, 1, 1, 2], [2, 1, 2]], strict=True
        )
    ]
    twins(gen3[2][2], gen3[2][3], "DZ")
    # The first two of each sibship marry out, the others more often than
    # not, and each couple has two to four children (IV). The twins stay
    # single: a spouse drawn beside one twin would sit inside their connector.
    couples3: list[tuple[dict, dict]] = []
    for sibship in gen3:
        for i, person in enumerate(sibship):
            if person["twins"]:
                continue
            if i < 2 or rng.random() < 0.6:
                couples3.append(marry_out(person))
    gen4 = {
        id(pat): children(pat, mat, [rng.choice([1, 2]) for _ in range(rng.randint(2, 4))])
        for pat, mat in couples3
    }
    sibships4 = list(gen4.values())
    twins(sibships4[2][0], sibships4[2][1], "MZ")
    sibships4[2][1]["sex"] = sibships4[2][0]["sex"]

    def descends_from(person: dict, ancestor: dict) -> bool:
        stack = [person]
        while stack:
            p = stack.pop()
            if p is ancestor:
                return True
            stack += [q for q in (p["father"], p["mother"]) if q is not None]
        return False

    def child_of(branch: dict, sex: int, skip: set[int]) -> dict:
        """A generation-IV descendant of `branch` (a generation-II person)."""
        for sibship in sibships4:
            for p in sibship:
                if id(p) not in skip and descends_from(p, branch) and not p["twins"]:
                    p["sex"] = sex
                    skip.add(id(p))
                    return p
        raise RuntimeError("no candidate")

    used: set[int] = set()
    # Loop 1, first cousins: two grandchildren of the first generation-II couple
    # through different children.
    a = child_of(gen3[0][0], 1, used)
    b = child_of(gen3[0][1], 2, used)
    # Loop 2, second cousins: descendants of the second and fourth branches.
    c = child_of(gen2[1], 1, used)
    d = child_of(gen2[3], 2, used)
    loop1, loop2 = (a, b), (c, d)
    # Two more couples in generation IV marry out.
    out4 = [marry_out(child_of(gen2[2], 1, used)), marry_out(child_of(gen2[0], 2, used))]
    # V: their children.
    gen5_loop1 = children(*loop1, [2, 1, 0, 1])
    gen5_loop2 = children(*loop2, [1, 2])
    for pat, mat in out4:
        children(pat, mat, [rng.choice([1, 2]) for _ in range(2)])

    # Markers and notation.
    for p in (f1, f2, gen2[0], couples2[2][1], gen3[3][2]):
        p["deceased"] = 1
        p["note"] = f"d. {rng.randint(55, 92)}y"
    gen5_loop1[0]["proband"] = 1
    gen5_loop1[2]["pregnancy"] = "SAB"
    gen5_loop1[2]["note"] = "10 wk"
    adoptee = sibships4[-1][-1]
    adoptee["adopted"] = "in"

    # Phenotypes: a recessive-looking condition, affected mostly below loops.
    loop_children = {id(p) for p in gen5_loop1 + gen5_loop2}
    for p in people:
        p["asd"] = (
            "2"
            if (id(p) in loop_children and rng.random() < 0.6)
            else ("2" if rng.random() < 0.08 else "1")
        )
        if p["gen"] <= 2 and rng.random() < 0.4:
            p["asd"] = "-9"
        affected = p["asd"] == "2"
        p["iq"] = (
            "-9"
            if p["gen"] <= 2
            else str(rng.randint(55, 85) if affected else rng.randint(85, 130))
        )
        p["carrier"] = "0"
    gen5_loop1[0]["asd"] = "2"
    gen5_loop1[0]["iq"] = "68"
    gen5_loop1[2]["asd"] = gen5_loop1[2]["iq"] = "-9"  # the miscarriage
    # Recessive story: the partners of each loop and one of the proband's
    # brothers are unaffected carriers.
    for p in (*loop1, *loop2, gen5_loop1[1]):
        if p["asd"] == "2":
            p["asd"], p["iq"] = "1", str(rng.randint(90, 125))
    for p in (*loop1, *loop2):
        p["carrier"] = "2"  # obligate: parents of an affected child
    for p in gen5_loop1 + gen5_loop2:
        if p["asd"] != "2" and p["pregnancy"] == "0":
            p["carrier"] = rng.choice(["0", "1"])
    adoptee["carrier"] = "0"
    gen5_loop1[1]["carrier"] = "1"

    # Number generation by generation, in file order.
    counters = [0] * 6
    for p in people:
        counters[p["gen"]] += 1
        p["id"] = f"{p['gen']}{counters[p['gen']]:02d}"

    header = [
        "FID",
        "IID",
        "Father",
        "Mother",
        "Sex",
        "Deceased",
        "Proband",
        "Twins",
        "b_ASD",
        "g_IQ",
        "cr_HEXA",
        "Pregnancy",
        "Adopted",
        "a_Note",
    ]
    rows = []
    for p in sorted(people, key=lambda q: q["gen"]):
        rows.append(
            [
                "FAM1",
                p["id"],
                p["father"]["id"] if p["father"] else "0",
                p["mother"]["id"] if p["mother"] else "0",
                p["sex"],
                p["deceased"],
                p["proband"],
                ",".join(f"{kind}:{q['id']}" for kind, q in p["twins"]) or "0",
                p["asd"],
                p["iq"],
                p["carrier"],
                p["pregnancy"],
                p["adopted"],
                p["note"],
            ]
        )
    return header, rows


def build_trio_with_siblings(
    fid: str,
    header: list[str],
    n_individuals: int,
    value_fn,
) -> tuple[list[str], list[list[str]]]:
    """Two parents and their children, `n_individuals` in total.

    `value_fn(rng, sid, is_child, proband=0)` returns the trailing phenotype
    fields, so each fixture decides its own column set.
    """
    rng = random.Random(SEED + len(fid))
    rows: list[list[str]] = []
    pat, mat = "002", "003"

    # Proband first, so ID 001 carries the arrow.
    rows.append([fid, "001", pat, mat, 1, *value_fn(rng, "001", True, proband=1)])
    rows.append([fid, pat, "0", "0", 1, *value_fn(rng, pat, False)])
    rows.append([fid, mat, "0", "0", 2, *value_fn(rng, mat, False)])
    for i in range(n_individuals - 3):
        sid = f"{i + 4:03d}"
        rows.append([fid, sid, pat, mat, 1 if i % 2 else 2, *value_fn(rng, sid, True)])
    return header, rows


def _maybe_missing(rng: random.Random, value: object, chance: float = 0.18) -> str:
    """Return the value, or the -9 no-data sentinel."""
    return "-9" if rng.random() < chance else str(value)


def _comma_decimal(value: float) -> str:
    """Format a score with a comma decimal separator (European locale)."""
    return f"{value:.10f}".replace(".", ",")


def main() -> None:
    print(f"Generating synthetic fixtures (seed {SEED})")

    header, rows = build_large_multi_union()
    _write(DATA_DIR / "complex/synesthesia.ped", header, rows)

    # Extended 10-column format plus four phenotype columns, one a PGS.
    def extended_values(rng, sid, is_child, proband=0):
        return [
            rng.choice([1, 2]),  # PHENO
            0,  # DOB
            0,  # TWINS
            0,  # DECEASED
            proband,
            rng.choice([0, 1]),  # epilepsy
            rng.randint(75, 125),  # IQ
            rng.randint(40, 90),  # SRS
            f"{rng.uniform(-1.5, 1.5):.2f}",  # PGS
        ]

    header, rows = build_trio_with_siblings(
        "FAM1",
        HEADER_EXTENDED,
        7,
        extended_values,
    )
    _write(DATA_DIR / "complex/large_family.ped", header, rows)

    # Two binary columns and two continuous ones, with no-data values.
    def pgs_values(rng, sid, is_child, proband=0):
        return [
            proband,
            _maybe_missing(rng, rng.choice([1, 2])),  # b_ASD
            _maybe_missing(rng, rng.choice([1, 2])),  # b_ADHD
            _maybe_missing(rng, rng.randint(70, 120)),  # IQ
            _maybe_missing(rng, rng.randint(45, 155)),  # SRS
            f"{rng.uniform(-1.2, 1.2):.3f}",  # PGS
        ]

    header, rows = build_trio_with_siblings("FAM1", HEADER_TWO_BINARY_PGS, 7, pgs_values)
    _write(DATA_DIR / "complex/multi_phenotype_pgs.ped", header, rows)

    # Every display mode at once: binary, AUTO-gradient, explicit gradient,
    # three PGS columns and two annotation columns. Scores use comma decimals.
    def all_modes_values(gene_col_name):
        def _fn(rng, sid, is_child, proband=0):
            genes = rng.sample(GENE_LABELS, k=rng.choice([0, 1, 1, 2]))
            return [
                proband,
                rng.choice([1, 2]),  # binary
                _maybe_missing(rng, rng.randint(70, 120), 0.12),  # IQ (AUTO -> gradient)
                _maybe_missing(rng, rng.randint(40, 95), 0.12),  # g_SRS
                _comma_decimal(rng.uniform(-1, 1)),
                _comma_decimal(rng.uniform(-1, 1)),
                _comma_decimal(rng.uniform(-1, 1)),
                genes[0] if genes else "",
                genes[1] if len(genes) > 1 else "",
            ]

        del gene_col_name
        return _fn

    header, rows = build_trio_with_siblings(
        "FAM1",
        HEADER_ALL_MODES,
        6,
        all_modes_values("PGS_int"),
    )
    _write(DATA_DIR / "phenotype/test_all.ped", header, rows)

    header, rows = build_trio_with_siblings(
        "FAM1",
        HEADER_ANNOTATIONS,
        7,
        all_modes_values("PGS_IQ"),
    )
    _write(DATA_DIR / "phenotype/annotations.ped", header, rows)

    header, rows = build_multigenerational()
    _write(DATA_DIR / "complex/multigenerational.ped", header, rows)


if __name__ == "__main__":
    main()
