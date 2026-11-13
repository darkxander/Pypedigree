"""Tests for pypedigree.models."""

from pypedigree.models import (
    EdgeLayout,
    EdgeType,
    Family,
    FamilyLayout,
    Individual,
    LineStyle,
    NodeLayout,
    Pedigree,
    Phenotype,
    Sex,
    TwinInfo,
    TwinType,
    Union,
)


class TestSexEnum:
    def test_from_raw_male(self):
        assert Sex.from_raw(1) == Sex.MALE
        assert Sex.from_raw("1") == Sex.MALE

    def test_from_raw_female(self):
        assert Sex.from_raw(2) == Sex.FEMALE
        assert Sex.from_raw("2") == Sex.FEMALE

    def test_from_raw_unknown_zero(self):
        assert Sex.from_raw(0) == Sex.UNKNOWN

    def test_from_raw_unknown_minus9(self):
        assert Sex.from_raw(-9) == Sex.UNKNOWN
        assert Sex.from_raw("-9") == Sex.UNKNOWN

    def test_from_raw_invalid(self):
        assert Sex.from_raw("abc") == Sex.UNKNOWN
        assert Sex.from_raw(None) == Sex.UNKNOWN


class TestPhenotypeEnum:
    def test_from_raw_unaffected(self):
        assert Phenotype.from_raw(1) == Phenotype.UNAFFECTED
        assert Phenotype.from_raw("1") == Phenotype.UNAFFECTED

    def test_from_raw_affected(self):
        assert Phenotype.from_raw(2) == Phenotype.AFFECTED
        assert Phenotype.from_raw("2") == Phenotype.AFFECTED

    def test_from_raw_unknown(self):
        assert Phenotype.from_raw(0) == Phenotype.UNKNOWN
        assert Phenotype.from_raw(-9) == Phenotype.UNKNOWN

    def test_from_raw_invalid(self):
        assert Phenotype.from_raw("abc") == Phenotype.UNKNOWN


class TestIndividual:
    def test_founder_no_parents(self):
        ind = Individual(sample_id="A", family_id="F1")
        assert ind.is_founder is True

    def test_not_founder_with_parents(self):
        ind = Individual(sample_id="C", family_id="F1", father_id="A", mother_id="B")
        assert ind.is_founder is False

    def test_not_founder_with_one_parent(self):
        ind = Individual(sample_id="C", family_id="F1", father_id="A")
        assert ind.is_founder is False

    def test_defaults(self):
        ind = Individual(sample_id="A", family_id="F1")
        assert ind.sex == Sex.UNKNOWN
        assert ind.phenotype == Phenotype.UNKNOWN
        assert ind.father_id is None
        assert ind.mother_id is None
        assert ind.date_of_birth is None
        assert ind.twin_info == []
        assert ind.is_deceased is False
        assert ind.is_proband is False
        assert ind.is_phantom is False
        assert ind.is_duplicate is False
        assert ind.duplicate_of is None

    def test_twin_info(self):
        twin = TwinInfo(twin_type=TwinType.MONOZYGOTIC, sibling_id="B")
        ind = Individual(sample_id="A", family_id="F1", twin_info=[twin])
        assert len(ind.twin_info) == 1
        assert ind.twin_info[0].twin_type == TwinType.MONOZYGOTIC
        assert ind.twin_info[0].sibling_id == "B"


class TestUnion:
    def test_basic_union(self):
        u = Union(partner_ids=("A", "B"), children_ids=["C", "D"])
        assert u.partner_ids == ("A", "B")
        assert u.children_ids == ["C", "D"]
        assert u.is_consanguinous is False
        assert u.is_broken is False

    def test_consanguinous_union(self):
        u = Union(partner_ids=("A", "B"), is_consanguinous=True)
        assert u.is_consanguinous is True


class TestFamily:
    def _make_family(self) -> Family:
        fam = Family(family_id="F1")
        fam.individuals["father"] = Individual(sample_id="father", family_id="F1", sex=Sex.MALE)
        fam.individuals["mother"] = Individual(sample_id="mother", family_id="F1", sex=Sex.FEMALE)
        fam.individuals["child"] = Individual(
            sample_id="child",
            family_id="F1",
            sex=Sex.FEMALE,
            father_id="father",
            mother_id="mother",
            phenotype=Phenotype.AFFECTED,
        )
        fam.unions.append(Union(partner_ids=("father", "mother"), children_ids=["child"]))
        return fam

    def test_get_founders(self):
        fam = self._make_family()
        founders = fam.get_founders()
        founder_ids = {f.sample_id for f in founders}
        assert founder_ids == {"father", "mother"}

    def test_get_children(self):
        fam = self._make_family()
        children = fam.get_children("father")
        assert len(children) == 1
        assert children[0].sample_id == "child"

    def test_get_children_by_mother(self):
        fam = self._make_family()
        children = fam.get_children("mother")
        assert len(children) == 1
        assert children[0].sample_id == "child"

    def test_get_union_for_parents(self):
        fam = self._make_family()
        union = fam.get_union_for_parents("father", "mother")
        assert union is not None
        assert "child" in union.children_ids

    def test_get_union_for_parents_reversed(self):
        fam = self._make_family()
        union = fam.get_union_for_parents("mother", "father")
        assert union is not None

    def test_get_union_for_parents_not_found(self):
        fam = self._make_family()
        assert fam.get_union_for_parents("father", "nonexistent") is None

    def test_get_unions_for_individual(self):
        fam = self._make_family()
        unions = fam.get_unions_for_individual("father")
        assert len(unions) == 1


class TestPedigree:
    def test_get_family(self):
        ped = Pedigree()
        fam = Family(family_id="F1")
        ped.families["F1"] = fam
        assert ped.get_family("F1") is fam

    def test_get_family_not_found(self):
        ped = Pedigree()
        assert ped.get_family("MISSING") is None

    def test_family_ids(self):
        ped = Pedigree()
        ped.families["F1"] = Family(family_id="F1")
        ped.families["F2"] = Family(family_id="F2")
        assert set(ped.family_ids) == {"F1", "F2"}


class TestLayoutTypes:
    def test_node_layout(self):
        node = NodeLayout(
            individual_id="A",
            x=1.0,
            y=2.0,
            sex=Sex.MALE,
            phenotype=Phenotype.UNAFFECTED,
            label="A",
        )
        assert node.x == 1.0
        assert node.y == 2.0
        assert node.is_deceased is False

    def test_edge_layout(self):
        edge = EdgeLayout(
            start=(0.0, 0.0),
            end=(1.0, 0.0),
            edge_type=EdgeType.RELATIONSHIP,
            line_style=LineStyle.SOLID,
        )
        assert edge.edge_type == EdgeType.RELATIONSHIP

    def test_family_layout_defaults(self):
        layout = FamilyLayout(family_id="F1")
        assert layout.nodes == []
        assert layout.edges == []
        assert layout.twin_connectors == []
        assert layout.duplicate_arcs == []
