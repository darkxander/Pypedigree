"""Shared test fixtures and configuration."""

from __future__ import annotations

import random
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib
import pytest

if TYPE_CHECKING:
    from pypedigree.models import Family

matplotlib.use("Agg")

TEST_DATA_DIR = Path(__file__).parent / "data"
BASIC_DIR = TEST_DATA_DIR / "basic"
FORMAT_DIR = TEST_DATA_DIR / "format"
EXTENDED_DIR = TEST_DATA_DIR / "extended"
TWINS_DIR = TEST_DATA_DIR / "twins"
CONSANG_DIR = TEST_DATA_DIR / "consanguinity"
PARTNERS_DIR = TEST_DATA_DIR / "multi_partner"
PHENOTYPE_DIR = TEST_DATA_DIR / "phenotype"
PLINK_DIR = TEST_DATA_DIR / "plink"
COMPLEX_DIR = TEST_DATA_DIR / "complex"


@pytest.fixture
def style():
    from pypedigree.config import PedigreeStyle

    return PedigreeStyle()


# ---------------------------------------------------------------------------
# Synthetic pedigree generators for stress testing
# ---------------------------------------------------------------------------


def generate_deep_pedigree(generations: int = 10, children_per_union: int = 2) -> Family:
    """Generate a deep linear pedigree.

    Creates a chain: founders → children → their children → etc.
    Each generation has one couple producing `children_per_union` children,
    and the first child continues the chain.
    """
    from pypedigree.models import Family, Individual, Sex, Union

    fam = Family(family_id="DEEP")
    unions: list[Union] = []
    prev_child_id: str | None = None

    for g in range(generations):
        if g == 0:
            father_id = "f_0"
            mother_id = "m_0"
        else:
            father_id = prev_child_id
            mother_id = f"spouse_{g}"
            fam.individuals[mother_id] = Individual(
                sample_id=mother_id,
                family_id="DEEP",
                sex=Sex.FEMALE,
            )

        if father_id not in fam.individuals:
            fam.individuals[father_id] = Individual(
                sample_id=father_id,
                family_id="DEEP",
                sex=Sex.MALE,
            )
        if mother_id not in fam.individuals:
            fam.individuals[mother_id] = Individual(
                sample_id=mother_id,
                family_id="DEEP",
                sex=Sex.FEMALE,
            )

        children = []
        for c in range(children_per_union):
            child_id = f"c_{g}_{c}"
            sex = Sex.MALE if c == 0 else Sex.FEMALE
            fam.individuals[child_id] = Individual(
                sample_id=child_id,
                family_id="DEEP",
                sex=sex,
                father_id=father_id,
                mother_id=mother_id,
            )
            children.append(child_id)

        unions.append(Union(partner_ids=(father_id, mother_id), children_ids=children))
        prev_child_id = children[0]

    fam.unions = unions
    return fam


def generate_wide_pedigree(num_couples: int = 20, children_per_couple: int = 3) -> Family:
    """Generate a wide single-generation pedigree.

    Creates `num_couples` founder couples, each with `children_per_couple` children.
    """
    from pypedigree.models import Family, Individual, Sex, Union

    fam = Family(family_id="WIDE")
    unions: list[Union] = []

    for i in range(num_couples):
        father_id = f"f_{i}"
        mother_id = f"m_{i}"
        fam.individuals[father_id] = Individual(sample_id=father_id, family_id="WIDE", sex=Sex.MALE)
        fam.individuals[mother_id] = Individual(
            sample_id=mother_id, family_id="WIDE", sex=Sex.FEMALE
        )

        children = []
        for c in range(children_per_couple):
            child_id = f"c_{i}_{c}"
            sex = Sex.MALE if c % 2 == 0 else Sex.FEMALE
            fam.individuals[child_id] = Individual(
                sample_id=child_id,
                family_id="WIDE",
                sex=sex,
                father_id=father_id,
                mother_id=mother_id,
            )
            children.append(child_id)

        unions.append(Union(partner_ids=(father_id, mother_id), children_ids=children))

    fam.unions = unions
    return fam


def generate_clinical_pedigree(
    target_size: int = 100, seed: int = 0, cousin_marriages: int = 0
) -> Family:
    """Generate a descent tree with partners from outside the family.

    Starts from one founder couple. Each child has 1-4 siblings and, with
    probability 0.6, a partner with no parents in the pedigree, as in most
    clinical pedigrees. ``cousin_marriages`` then joins pairs of unpartnered
    relatives of one generation from different sibships (consanguinous
    loops), each with 1-3 children.
    """
    from pypedigree.models import Family, Individual, Sex, Union

    def new(fam: Family, sex: Sex, father: str | None = None, mother: str | None = None) -> str:
        sid = f"i{len(fam.individuals)}"
        fam.individuals[sid] = Individual(
            sample_id=sid, family_id="CLIN", sex=sex, father_id=father, mother_id=mother
        )
        return sid

    for attempt in range(100):
        rng = random.Random(seed * 1000 + attempt)
        fam = Family(family_id="CLIN")
        couples = [(new(fam, Sex.MALE), new(fam, Sex.FEMALE))]
        generation_of: dict[str, int] = {}
        while couples and len(fam.individuals) < target_size:
            next_couples = []
            for dad, mum in couples:
                kids = []
                for _ in range(rng.randint(1, 4)):
                    if len(fam.individuals) >= target_size:
                        break
                    kids.append(new(fam, rng.choice([Sex.MALE, Sex.FEMALE]), dad, mum))
                if kids:
                    fam.unions.append(Union(partner_ids=(dad, mum), children_ids=kids))
                for kid in kids:
                    generation_of[kid] = generation_of.get(dad, 0) + 1
                    if len(fam.individuals) < target_size and rng.random() < 0.6:
                        male = fam.individuals[kid].sex == Sex.MALE
                        partner = new(fam, Sex.FEMALE if male else Sex.MALE)
                        generation_of[partner] = generation_of[kid]
                        next_couples.append((kid, partner) if male else (partner, kid))
            couples = next_couples
        if len(fam.individuals) >= target_size:
            break

    partnered = {p for u in fam.unions for p in u.partner_ids}
    sibship = {c: i for i, u in enumerate(fam.unions) for c in u.children_ids}
    pool = [s for s in sibship if s not in partnered]
    made = 0
    for _ in range(2000):
        if made >= cousin_marriages or len(pool) < 2:
            break
        a, b = rng.sample(pool, 2)
        ind_a, ind_b = fam.individuals[a], fam.individuals[b]
        if (
            generation_of[a] != generation_of[b]
            or sibship[a] == sibship[b]
            or ind_a.sex == ind_b.sex
        ):
            continue
        dad, mum = (a, b) if ind_a.sex == Sex.MALE else (b, a)
        kids = [
            new(fam, rng.choice([Sex.MALE, Sex.FEMALE]), dad, mum) for _ in range(rng.randint(1, 3))
        ]
        fam.unions.append(Union(partner_ids=(dad, mum), children_ids=kids, is_consanguinous=True))
        pool.remove(a)
        pool.remove(b)
        made += 1
    return fam


def generate_realistic_pedigree(target_size: int = 100, seed: int = 42) -> Family:
    """Generate a multi-generation pedigree with random mating.

    Creates a realistic pedigree by starting with founder couples and
    randomly mating children across generations until target_size is reached.
    """
    from pypedigree.models import Family, Individual, Sex, Union

    rng = random.Random(seed)
    fam = Family(family_id="REAL")
    unions: list[Union] = []

    # Start with founder couples
    num_founders = max(4, target_size // 15)
    if num_founders % 2 != 0:
        num_founders += 1

    males: list[str] = []
    females: list[str] = []

    for i in range(num_founders):
        sid = f"founder_{i}"
        sex = Sex.MALE if i % 2 == 0 else Sex.FEMALE
        fam.individuals[sid] = Individual(sample_id=sid, family_id="REAL", sex=sex)
        if sex == Sex.MALE:
            males.append(sid)
        else:
            females.append(sid)

    # Create initial couples
    available_males = list(males)
    available_females = list(females)
    rng.shuffle(available_males)
    rng.shuffle(available_females)

    gen = 0

    spouse_counter = 0
    while len(fam.individuals) < target_size and (available_males or available_females):
        # Inject outside spouses when one pool is exhausted
        while len(available_males) < len(available_females) and len(fam.individuals) < target_size:
            sid = f"spouse_m_{spouse_counter}"
            spouse_counter += 1
            fam.individuals[sid] = Individual(sample_id=sid, family_id="REAL", sex=Sex.MALE)
            available_males.append(sid)

        while len(available_females) < len(available_males) and len(fam.individuals) < target_size:
            sid = f"spouse_f_{spouse_counter}"
            spouse_counter += 1
            fam.individuals[sid] = Individual(sample_id=sid, family_id="REAL", sex=Sex.FEMALE)
            available_females.append(sid)

        if not available_males or not available_females:
            break

        pairs = min(len(available_males), len(available_females))
        new_children_m: list[str] = []
        new_children_f: list[str] = []

        for p in range(pairs):
            if len(fam.individuals) >= target_size:
                break
            father_id = available_males[p]
            mother_id = available_females[p]
            num_children = rng.randint(2, 4)
            children = []
            for c in range(num_children):
                if len(fam.individuals) >= target_size:
                    break
                child_id = f"g{gen}_c{p}_{c}"
                sex = rng.choice([Sex.MALE, Sex.FEMALE])
                fam.individuals[child_id] = Individual(
                    sample_id=child_id,
                    family_id="REAL",
                    sex=sex,
                    father_id=father_id,
                    mother_id=mother_id,
                )
                children.append(child_id)
                if sex == Sex.MALE:
                    new_children_m.append(child_id)
                else:
                    new_children_f.append(child_id)

            if children:
                unions.append(
                    Union(
                        partner_ids=(father_id, mother_id),
                        children_ids=children,
                    )
                )

        gen += 1
        available_males = new_children_m
        available_females = new_children_f
        rng.shuffle(available_males)
        rng.shuffle(available_females)

    fam.unions = unions
    return fam
