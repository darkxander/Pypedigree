"""Tests for pypedigree.parser."""

import pytest
from conftest import TEST_DATA_DIR

from pypedigree.errors import PedigreeParseError, PedigreeValidationError
from pypedigree.models import Phenotype, PhenotypeDisplayMode, Sex
from pypedigree.parser import (
    _detect_header,
    _detect_separator,
    _strip_display_prefix,
    parse_pedigree,
    parse_pedigree_from_string,
)


class TestSeparatorDetection:
    def test_tab_separated(self):
        assert _detect_separator("FAM1\tfather\t0\t0\t1\t1") == "\t"

    def test_comma_separated(self):
        assert _detect_separator("FAM1,father,0,0,1,1") == ","

    def test_space_separated(self):
        assert _detect_separator("FAM1 father 0 0 1 1") is None

    def test_too_few_fields(self):
        with pytest.raises(PedigreeParseError, match="fewer than 6"):
            _detect_separator("FAM1 father")


class TestHeaderDetection:
    def test_detects_fid_header(self):
        assert _detect_header("fid\tiid\tfather\tmother\tsex\tphenotype", "\t") is True

    def test_detects_family_id_header(self):
        assert _detect_header("family_id\tsample_id\tpat\tmat\tsex\tphenotype", "\t") is True

    def test_no_header(self):
        assert _detect_header("FAM1\tfather\t0\t0\t1\t1", "\t") is False

    def test_numeric_first_field_not_header(self):
        assert _detect_header("1\tfather\t0\t0\t1\t1", "\t") is False


class TestParsePedigree:
    def test_simple_trio(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        assert "FAM1" in ped.family_ids
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3
        assert fam.individuals["father"].sex == Sex.MALE
        assert fam.individuals["mother"].sex == Sex.FEMALE
        assert fam.individuals["child"].sex == Sex.FEMALE
        assert fam.individuals["child"].phenotype == Phenotype.AFFECTED
        assert fam.individuals["child"].father_id == "father"
        assert fam.individuals["child"].mother_id == "mother"

    def test_two_gen_siblings(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/two_gen_siblings.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 5
        children = fam.get_children("father")
        assert len(children) == 3

    def test_three_gen(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 5
        assert fam.individuals["father"].father_id == "gfather"
        assert fam.individuals["father"].mother_id == "gmother"

    def test_comma_separated(self):
        ped = parse_pedigree(TEST_DATA_DIR / "format/comma_separated.csv")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3

    def test_space_separated(self):
        ped = parse_pedigree(TEST_DATA_DIR / "format/space_separated.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3

    def test_with_header(self):
        ped = parse_pedigree(TEST_DATA_DIR / "format/with_header.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3
        # Header should not be parsed as an individual
        assert "fid" not in fam.individuals

    def test_multi_family(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/multi_family.ped")
        assert len(ped.family_ids) == 2
        assert "FAM1" in ped.family_ids
        assert "FAM2" in ped.family_ids
        assert len(ped.get_family("FAM1").individuals) == 3
        assert len(ped.get_family("FAM2").individuals) == 3

    def test_affected_mixed(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/affected_mixed.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].phenotype == Phenotype.UNAFFECTED
        assert fam.individuals["mother"].phenotype == Phenotype.AFFECTED
        assert fam.individuals["son"].phenotype == Phenotype.AFFECTED
        assert fam.individuals["daughter"].phenotype == Phenotype.UNAFFECTED
        assert fam.individuals["unknown_child"].phenotype == Phenotype.UNKNOWN

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_pedigree("/nonexistent/file.ped")

    def test_malformed_missing_cols(self):
        with pytest.raises(PedigreeParseError, match="fewer than 6"):
            parse_pedigree(TEST_DATA_DIR / "format/malformed_missing_cols.ped")


class TestUnionInference:
    def test_trio_has_one_union(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.unions) == 1
        union = fam.unions[0]
        assert set(union.partner_ids) == {"father", "mother"}
        assert union.children_ids == ["child"]

    def test_siblings_share_union(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/two_gen_siblings.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.unions) == 1
        union = fam.unions[0]
        assert len(union.children_ids) == 3

    def test_three_gen_has_two_unions(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.unions) == 2

    def test_union_male_first(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        union = fam.unions[0]
        assert union.partner_ids[0] == "father"
        assert union.partner_ids[1] == "mother"


class TestPhantomParents:
    def test_single_parent_creates_phantom(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/single_parent.ped")
        fam = ped.get_family("FAM1")
        # Should have: mother, child, phantom_father
        assert len(fam.individuals) == 3
        child = fam.individuals["child"]
        assert child.mother_id == "mother"
        assert child.father_id is not None
        phantom = fam.individuals[child.father_id]
        assert phantom.is_phantom is True
        assert phantom.sex == Sex.MALE

    def test_phantom_has_union(self):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/single_parent.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.unions) == 1


class TestExtendedFormat:
    def test_parse_deceased_column(self):
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_deceased.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].is_deceased is True
        assert fam.individuals["mother"].is_deceased is False
        assert fam.individuals["child"].is_deceased is False

    def test_parse_proband_column(self):
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_proband.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].is_proband is False
        assert fam.individuals["mother"].is_proband is False
        assert fam.individuals["child"].is_proband is True

    def test_extended_full_format(self):
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_full.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["gfather"].is_deceased is True
        assert fam.individuals["gfather"].date_of_birth == "19400515"
        assert fam.individuals["child1"].is_proband is True
        assert fam.individuals["child1"].phenotype == Phenotype.AFFECTED
        assert fam.individuals["child2"].is_proband is False
        assert fam.individuals["child2"].is_deceased is False

    def test_extended_with_missing_values(self):
        """0 and -9 in extended columns should be treated as missing/default."""
        text = "FAM1\tf\t0\t0\t1\t1\t-9\t-9\t-9\t-9\n"
        ped = parse_pedigree_from_string(text)
        fam = ped.get_family("FAM1")
        assert fam.individuals["f"].is_deceased is False
        assert fam.individuals["f"].is_proband is False
        assert fam.individuals["f"].date_of_birth is None

    def test_six_column_still_works(self):
        """Existing 6-column files must continue to parse correctly."""
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].is_deceased is False
        assert fam.individuals["father"].is_proband is False
        assert fam.individuals["father"].date_of_birth is None

    def test_partial_extended_9_columns(self):
        """File with 9 columns (no proband column) still works."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_partial_9col.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].is_deceased is True
        assert fam.individuals["mother"].is_deceased is False
        # No proband column — should default to False
        assert fam.individuals["father"].is_proband is False

    def test_date_of_birth_parsed(self):
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_full.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["gfather"].date_of_birth == "19400515"
        assert fam.individuals["gmother"].date_of_birth == "19420301"
        assert fam.individuals["mother"].date_of_birth == "19680101"


class TestParseFromString:
    def test_basic_string(self):
        text = "FAM1\tf\t0\t0\t1\t1\nFAM1\tm\t0\t0\t2\t1\nFAM1\tc\tf\tm\t2\t2\n"
        ped = parse_pedigree_from_string(text)
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3

    def test_empty_string(self):
        with pytest.raises(PedigreeParseError, match="empty"):
            parse_pedigree_from_string("")

    def test_whitespace_only(self):
        with pytest.raises(PedigreeParseError, match="empty"):
            parse_pedigree_from_string("   \n  \n  ")


class TestGetAncestors:
    def test_founder_has_no_ancestors(self):
        """Founders (no parents) should have empty ancestor set."""
        from pypedigree.parser import _get_ancestors

        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        ancestors = _get_ancestors(fam, "father")
        assert ancestors == set()

    def test_child_has_parents_as_ancestors(self):
        """Child should have both parents as ancestors."""
        from pypedigree.parser import _get_ancestors

        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        ancestors = _get_ancestors(fam, "child")
        assert "father" in ancestors
        assert "mother" in ancestors

    def test_grandchild_has_all_ancestors(self):
        """Grandchild should have parents and grandparents."""
        from pypedigree.parser import _get_ancestors

        ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
        fam = ped.get_family("FAM1")
        ancestors = _get_ancestors(fam, "child")
        assert ancestors == {"father", "mother", "gfather", "gmother"}

    def test_nonexistent_individual(self):
        """Nonexistent individual should return empty set."""
        from pypedigree.parser import _get_ancestors

        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        ancestors = _get_ancestors(fam, "nonexistent")
        assert ancestors == set()


class TestConsanguinityDetection:
    def test_first_cousin_marriage_detected(self):
        """First cousins sharing grandparents -> consanguinous union."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        fam = ped.get_family("FAM1")
        cousin_union = fam.get_union_for_parents("cousin1", "cousin2")
        assert cousin_union is not None
        assert cousin_union.is_consanguinous is True

    def test_first_cousin_other_unions_not_consanguinous(self):
        """Non-consanguinous unions in the same family should remain unmarked."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        fam = ped.get_family("FAM1")
        for union in fam.unions:
            if set(union.partner_ids) != {"cousin1", "cousin2"}:
                assert union.is_consanguinous is False, (
                    f"Union {union.partner_ids} should not be consanguinous"
                )

    def test_uncle_niece_detected(self):
        """Uncle-niece marriage -> consanguinous (cross-generational)."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_uncle_niece.ped")
        fam = ped.get_family("FAM1")
        uncle_niece_union = fam.get_union_for_parents("uncle", "daughter")
        assert uncle_niece_union is not None
        assert uncle_niece_union.is_consanguinous is True

    def test_sibling_marriage_detected(self):
        """Sibling marriage -> consanguinous (direct shared parents)."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_sibling_marriage.ped")
        fam = ped.get_family("FAM1")
        sibling_union = fam.get_union_for_parents("son", "daughter")
        assert sibling_union is not None
        assert sibling_union.is_consanguinous is True

    def test_no_consanguinity_in_simple_pedigree(self):
        """Three-gen pedigree with no shared ancestors -> no consanguinity."""
        ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
        fam = ped.get_family("FAM1")
        for union in fam.unions:
            assert union.is_consanguinous is False

    def test_consanguinity_from_string(self):
        """Consanguinity detection works with parse_pedigree_from_string."""
        text = (
            "FAM1\tgf\t0\t0\t1\t1\n"
            "FAM1\tgm\t0\t0\t2\t1\n"
            "FAM1\tson\tgf\tgm\t1\t1\n"
            "FAM1\tdaughter\tgf\tgm\t2\t1\n"
            "FAM1\tchild\tson\tdaughter\t1\t2\n"
        )
        ped = parse_pedigree_from_string(text)
        fam = ped.get_family("FAM1")
        sibling_union = fam.get_union_for_parents("son", "daughter")
        assert sibling_union is not None
        assert sibling_union.is_consanguinous is True


class TestTwinParsing:
    def test_mz_twins_parsed(self):
        """MZ twin declarations are parsed into TwinInfo objects."""
        from pypedigree.models import TwinType

        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")
        twin1 = fam.individuals["twin1"]
        assert len(twin1.twin_info) == 1
        assert twin1.twin_info[0].twin_type == TwinType.MONOZYGOTIC
        assert twin1.twin_info[0].sibling_id == "twin2"

    def test_dz_twins_parsed(self):
        """DZ twin declarations are parsed into TwinInfo objects."""
        from pypedigree.models import TwinType

        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_dz.ped")
        fam = ped.get_family("FAM1")
        twin1 = fam.individuals["twin1"]
        assert len(twin1.twin_info) == 1
        assert twin1.twin_info[0].twin_type == TwinType.DIZYGOTIC

    def test_unknown_twins_parsed(self):
        """TW (unknown zygosity) twin declarations are parsed."""
        from pypedigree.models import TwinType

        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_unknown.ped")
        fam = ped.get_family("FAM1")
        twin1 = fam.individuals["twin1"]
        assert len(twin1.twin_info) == 1
        assert twin1.twin_info[0].twin_type == TwinType.UNKNOWN

    def test_triplets_parsed(self):
        """Triplet declarations with multiple comma-separated entries."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_triplets.ped")
        fam = ped.get_family("FAM1")
        trip1 = fam.individuals["trip1"]
        assert len(trip1.twin_info) == 2
        sibling_ids = {ti.sibling_id for ti in trip1.twin_info}
        assert sibling_ids == {"trip2", "trip3"}

    def test_one_sided_twin_resolved_bidirectionally(self):
        """One-sided twin declaration is resolved bidirectionally."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_one_sided.ped")
        fam = ped.get_family("FAM1")
        # twin1 declares MZ:twin2
        twin1 = fam.individuals["twin1"]
        assert len(twin1.twin_info) == 1
        # twin2 had 0 in column 8, but should get reciprocal
        twin2 = fam.individuals["twin2"]
        assert len(twin2.twin_info) == 1
        assert twin2.twin_info[0].sibling_id == "twin1"

    def test_no_twin_info_when_missing(self):
        """Individuals with 0 in twin column have empty twin_info."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].twin_info == []
        assert fam.individuals["mother"].twin_info == []

    def test_invalid_twin_format_raises(self):
        """Invalid twin entry (missing colon) raises PedigreeParseError."""
        text = "FAM1\tf\t0\t0\t1\t1\t0\tINVALID\t0\t0\n"
        with pytest.raises(PedigreeParseError, match="invalid twin entry"):
            parse_pedigree_from_string(text)

    def test_unknown_twin_type_raises(self):
        """Unknown twin type (e.g., XX:id) raises PedigreeParseError."""
        text = "FAM1\tf\t0\t0\t1\t1\t0\tXX:sib\t0\t0\n"
        with pytest.raises(PedigreeParseError, match="unknown twin type"):
            parse_pedigree_from_string(text)

    def test_six_column_file_no_twins(self):
        """Existing 6-column files still produce empty twin_info."""
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        for ind in fam.individuals.values():
            assert ind.twin_info == []

    def test_dob_parsed_with_twins(self):
        """DOB in column 7 is parsed alongside twin info in column 8."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_with_siblings.ped")
        fam = ped.get_family("FAM1")
        assert fam.individuals["older"].date_of_birth == "19900101"
        assert fam.individuals["twin1"].date_of_birth == "19930601"
        assert fam.individuals["younger"].date_of_birth == "19960301"


def _tsv(rows: str) -> str:
    """Tab-separated text from whitespace-aligned rows; '-' is an empty field."""
    lines = [r.split() for r in rows.strip().splitlines()]
    return "\n".join("\t".join("" if f == "-" else f for f in fields) for fields in lines) + "\n"


class TestHeaderColumnsOutOfOrder:
    """A header naming structural columns outside the standard order is read by name.

    Such headers used to be read by position whenever they had no extra
    phenotype, annotation or notation column: a ``Twins`` column given as
    column 7 was read as the date of birth and the twins were silently lost.
    """

    def test_twins_as_seventh_column(self):
        from pypedigree.models import TwinType

        text = _tsv("""
            FID IID Father Mother Sex Phenotype Twins
            F   f   0      0      1   1         -
            F   m   0      0      2   1         -
            F   t1  f      m      1   1         DZ:t2
            F   t2  f      m      2   1         DZ:t1
        """)
        fam = parse_pedigree_from_string(text).get_family("F")
        t1 = fam.individuals["t1"]
        assert [(ti.twin_type, ti.sibling_id) for ti in t1.twin_info] == [
            (TwinType.DIZYGOTIC, "t2")
        ]
        assert t1.date_of_birth is None
        assert fam.individuals["t2"].twin_info[0].sibling_id == "t1"

    def test_deceased_without_twins_column(self):
        text = _tsv("""
            FID IID Father Mother Sex Phenotype DOB  Deceased
            F   f   0      0      1   1         1950 1
            F   m   0      0      2   1         1952 0
        """)
        fam = parse_pedigree_from_string(text).get_family("F")
        assert fam.individuals["f"].is_deceased
        assert not fam.individuals["m"].is_deceased
        assert fam.individuals["f"].date_of_birth == "1950"

    def test_reordered_core_columns(self):
        text = _tsv("""
            IID FID Sex Father Mother Phenotype
            f   F   1   0      0      1
            m   F   2   0      0      1
            c   F   2   f      m      2
        """)
        c = parse_pedigree_from_string(text).get_family("F").individuals["c"]
        assert (c.sex, c.father_id, c.mother_id) == (Sex.FEMALE, "f", "m")
        assert c.phenotype == Phenotype.AFFECTED

    def test_no_phenotype_column(self):
        text = _tsv("""
            FID IID Father Mother Sex Proband
            F   f   0      0      1   0
            F   m   0      0      2   0
            F   c   f      m      1   1
        """)
        c = parse_pedigree_from_string(text).get_family("F").individuals["c"]
        assert c.is_proband
        assert c.phenotype == Phenotype.UNKNOWN

    def test_standard_order_header_unchanged(self):
        """A header in the standard order reads the same as the headerless file."""
        rows = """
            F f  0 0 1 1 19600101 0     1 0
            F m  0 0 2 1 19620101 0     0 0
            F t1 f m 1 2 19900101 MZ:t2 0 1
            F t2 f m 1 1 19900101 MZ:t1 0 0
        """
        header = "FID IID Father Mother Sex Phenotype DOB Twins Deceased Proband\n"
        bare = parse_pedigree_from_string(_tsv(rows)).get_family("F")
        named = parse_pedigree_from_string(_tsv(header + rows)).get_family("F")
        assert named.individuals == bare.individuals


class TestDuplicateNodes:
    def test_female_multi_partner_creates_duplicate(self):
        """Female with two partners should create one duplicate node."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        duplicates = [ind for ind in fam.individuals.values() if ind.is_duplicate]
        assert len(duplicates) == 1
        dup = duplicates[0]
        assert dup.duplicate_of == "mother"
        assert dup.sex == Sex.FEMALE

    def test_male_multi_partner_creates_duplicate(self):
        """Male with two partners should create one duplicate node."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_male.ped")
        fam = ped.get_family("FAM1")
        duplicates = [ind for ind in fam.individuals.values() if ind.is_duplicate]
        assert len(duplicates) == 1
        dup = duplicates[0]
        assert dup.duplicate_of == "father"
        assert dup.sex == Sex.MALE

    def test_duplicate_union_updated(self):
        """Union partner_ids should reference duplicate, not original."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        # First union: original mother with father1
        union1 = fam.unions[0]
        assert "mother" in union1.partner_ids
        # Second union: duplicate mother with father2
        union2 = fam.unions[1]
        assert "mother" not in union2.partner_ids
        dup_ids = [sid for sid in union2.partner_ids if sid.startswith("_dup_")]
        assert len(dup_ids) == 1

    def test_children_parent_refs_updated(self):
        """Children in duplicate union should reference the duplicate, not original."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        child2 = fam.individuals["child2"]
        # child2's mother should be the duplicate, not original
        assert child2.mother_id != "mother"
        assert child2.mother_id.startswith("_dup_")

    def test_duplicate_preserves_attributes(self):
        """Duplicate should copy sex and phenotype from original."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        original = fam.individuals["mother"]
        dup = [ind for ind in fam.individuals.values() if ind.is_duplicate][0]
        assert dup.sex == original.sex
        assert dup.phenotype == original.phenotype

    def test_duplicate_keeps_phenotype_values(self):
        """A duplicate is drawn like the original, bands included."""
        text = (
            "FID\tIID\tFather\tMother\tSex\tb_ASD\tg_IQ\tc_Site\n"
            "F\tf\t0\t0\t1\t2\t95\tParis\n"
            "F\tm1\t0\t0\t2\t1\t110\tLyon\n"
            "F\tm2\t0\t0\t2\t1\t120\tLyon\n"
            "F\tc1\tf\tm1\t1\t1\t100\tLyon\n"
            "F\tc2\tf\tm2\t1\t1\t105\tLyon\n"
        )
        fam = parse_pedigree_from_string(text).get_family("F")
        original = fam.individuals["f"]
        dup = next(ind for ind in fam.individuals.values() if ind.is_duplicate)
        assert dup.extra_phenotypes == original.extra_phenotypes
        assert (
            dup.extra_phenotypes_raw
            == original.extra_phenotypes_raw
            == {"ASD": 2.0, "IQ": 95.0, "Site": None}
        )
        assert dup.extra_phenotypes_categorical == {"Site": "Paris"}

    def test_no_duplicates_for_single_union(self):
        """Individuals with single union should not be duplicated."""
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        duplicates = [ind for ind in fam.individuals.values() if ind.is_duplicate]
        assert len(duplicates) == 0

    def test_complex_multi_partner(self):
        """Multi-gen pedigree: father has two unions, one via parents."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_complex.ped")
        fam = ped.get_family("FAM1")
        duplicates = [ind for ind in fam.individuals.values() if ind.is_duplicate]
        assert len(duplicates) == 1
        assert duplicates[0].duplicate_of == "father"

    def test_duplicate_proband_not_copied(self):
        """Proband marker should only be on original, not duplicate."""
        text = (
            "FAM1\tf\t0\t0\t1\t1\t0\t0\t0\t1\n"  # proband
            "FAM1\tm1\t0\t0\t2\t1\t0\t0\t0\t0\n"
            "FAM1\tm2\t0\t0\t2\t1\t0\t0\t0\t0\n"
            "FAM1\tc1\tf\tm1\t1\t1\t0\t0\t0\t0\n"
            "FAM1\tc2\tf\tm2\t1\t1\t0\t0\t0\t0\n"
        )
        ped = parse_pedigree_from_string(text)
        fam = ped.get_family("FAM1")
        original = fam.individuals["f"]
        assert original.is_proband is True
        dup = [ind for ind in fam.individuals.values() if ind.is_duplicate][0]
        assert dup.is_proband is False


class TestEdgeCases:
    def test_cycle_detection(self):
        """A cycle in parent-child relationships should raise PedigreeValidationError."""
        # a's parent is c, c's parent is a → mutual cycle
        text_cycle = (
            "FAM1\tm1\t0\t0\t2\t1\n"
            "FAM1\tm2\t0\t0\t2\t1\n"
            "FAM1\ta\tc\tm1\t1\t1\n"
            "FAM1\tc\ta\tm2\t1\t1\n"
        )
        with pytest.raises(PedigreeValidationError, match="Cycle detected"):
            parse_pedigree_from_string(text_cycle)

    def test_unicode_sample_ids(self):
        """Pedigree with unicode IDs should parse correctly."""
        ped = parse_pedigree(TEST_DATA_DIR / "format/unicode_ids.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3
        assert "père" in fam.individuals
        assert "mère" in fam.individuals
        assert "enfant" in fam.individuals
        assert fam.individuals["enfant"].father_id == "père"
        assert fam.individuals["enfant"].mother_id == "mère"
        assert fam.individuals["enfant"].phenotype == Phenotype.AFFECTED

    def test_single_individual_no_parents(self):
        """A single founder with no parents or children should parse successfully."""
        text = "FAM1\tfounder\t0\t0\t1\t1\n"
        ped = parse_pedigree_from_string(text)
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 1
        assert fam.individuals["founder"].father_id is None
        assert fam.individuals["founder"].mother_id is None
        assert len(fam.unions) == 0

    def test_header_only_file(self, tmp_path):
        """File with only a header line and no data should raise PedigreeParseError."""
        p = tmp_path / "header_only.ped"
        p.write_text("fid\tiid\tfather\tmother\tsex\tphenotype\n")
        with pytest.raises(PedigreeParseError, match="only a header"):
            parse_pedigree(p)

    def test_duplicate_sample_id(self):
        """Two rows with same family_id + sample_id should raise PedigreeParseError."""
        text = "FAM1\ta\t0\t0\t1\t1\nFAM1\ta\t0\t0\t2\t1\n"
        with pytest.raises(PedigreeParseError, match="duplicate sample_id"):
            parse_pedigree_from_string(text)

    def test_empty_string_raises(self):
        """Empty string should raise PedigreeParseError."""
        with pytest.raises(PedigreeParseError, match="empty"):
            parse_pedigree_from_string("")

    def test_whitespace_only_string_raises(self):
        """Whitespace-only string should raise PedigreeParseError."""
        with pytest.raises(PedigreeParseError, match="empty"):
            parse_pedigree_from_string("  \n  \n  ")


class TestPandasIntegration:
    @pytest.fixture(autouse=True)
    def _require_pandas(self):
        pytest.importorskip("pandas")

    def test_dataframe_simple_trio(self):
        """Parse a simple trio from a DataFrame."""
        import pandas as pd

        from pypedigree.parser import parse_pedigree_from_dataframe

        df = pd.DataFrame(
            {
                "fid": ["FAM1", "FAM1", "FAM1"],
                "iid": ["father", "mother", "child"],
                "pat": ["0", "0", "father"],
                "mat": ["0", "0", "mother"],
                "sex": [1, 2, 2],
                "pheno": [1, 1, 2],
            }
        )
        ped = parse_pedigree_from_dataframe(df)
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3
        assert fam.individuals["child"].father_id == "father"
        assert fam.individuals["child"].mother_id == "mother"

    def test_dataframe_with_nan(self):
        """DataFrame with NaN/None for missing parents should parse as founders."""
        import numpy as np
        import pandas as pd

        from pypedigree.parser import parse_pedigree_from_dataframe

        df = pd.DataFrame(
            {
                "fid": ["FAM1", "FAM1", "FAM1"],
                "iid": ["father", "mother", "child"],
                "pat": [np.nan, np.nan, "father"],
                "mat": [np.nan, np.nan, "mother"],
                "sex": [1, 2, 2],
                "pheno": [1, 1, 2],
            }
        )
        ped = parse_pedigree_from_dataframe(df)
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].father_id is None
        assert fam.individuals["mother"].mother_id is None
        assert fam.individuals["child"].father_id == "father"

    def test_dataframe_extended_columns(self):
        """10-column DataFrame with deceased/proband flags should parse correctly."""
        import pandas as pd

        from pypedigree.parser import parse_pedigree_from_dataframe

        df = pd.DataFrame(
            {
                "fid": ["FAM1", "FAM1", "FAM1"],
                "iid": ["father", "mother", "child"],
                "pat": ["0", "0", "father"],
                "mat": ["0", "0", "mother"],
                "sex": [1, 2, 2],
                "pheno": [1, 1, 2],
                "dob": ["19600101", "0", "19900101"],
                "twins": ["0", "0", "0"],
                "deceased": [1, 0, 0],
                "proband": [0, 0, 1],
            }
        )
        ped = parse_pedigree_from_dataframe(df)
        fam = ped.get_family("FAM1")
        assert fam.individuals["father"].is_deceased is True
        assert fam.individuals["child"].is_proband is True
        assert fam.individuals["father"].date_of_birth == "19600101"

    def test_dataframe_too_few_columns(self):
        """DataFrame with fewer than 6 columns should raise PedigreeParseError."""
        import pandas as pd

        from pypedigree.parser import parse_pedigree_from_dataframe

        df = pd.DataFrame(
            {
                "fid": ["FAM1"],
                "iid": ["father"],
                "pat": ["0"],
                "mat": ["0"],
            }
        )
        with pytest.raises(PedigreeParseError, match="at least 6 columns"):
            parse_pedigree_from_dataframe(df)

    def test_dataframe_empty(self):
        """Empty DataFrame (0 rows) should raise PedigreeParseError."""
        import pandas as pd

        from pypedigree.parser import parse_pedigree_from_dataframe

        df = pd.DataFrame(columns=["fid", "iid", "pat", "mat", "sex", "pheno"])
        with pytest.raises(PedigreeParseError, match="empty"):
            parse_pedigree_from_dataframe(df)


class TestPhenotypeDisplayPrefixes:
    """Tests for prefix-based phenotype display mode detection."""

    def test_strip_display_prefix_binary(self):
        name, mode = _strip_display_prefix("b_ASD")
        assert name == "ASD"
        assert mode == PhenotypeDisplayMode.BINARY

    def test_strip_display_prefix_gradient(self):
        name, mode = _strip_display_prefix("g_IQ")
        assert name == "IQ"
        assert mode == PhenotypeDisplayMode.GRADIENT

    def test_strip_display_prefix_outline(self):
        name, mode = _strip_display_prefix("o_Proband")
        assert name == "Proband"
        assert mode == PhenotypeDisplayMode.OUTLINE

    def test_strip_display_prefix_categorical(self):
        name, mode = _strip_display_prefix("c_Diagnosis")
        assert name == "Diagnosis"
        assert mode == PhenotypeDisplayMode.CATEGORICAL

    def test_strip_display_prefix_no_prefix(self):
        name, mode = _strip_display_prefix("epilepsy")
        assert name == "epilepsy"
        assert mode == PhenotypeDisplayMode.AUTO

    def test_strip_display_prefix_case_insensitive(self):
        name, mode = _strip_display_prefix("B_ASD")
        assert name == "ASD"
        assert mode == PhenotypeDisplayMode.BINARY

    def test_strip_display_prefix_short_name(self):
        """Prefix alone (e.g. 'b_') without a name should not match."""
        name, mode = _strip_display_prefix("b_")
        assert name == "b_"
        assert mode == PhenotypeDisplayMode.AUTO

    def test_parse_prefix_modes_file(self):
        """Parse the prefix_modes.ped test file."""
        ped = parse_pedigree(str(TEST_DATA_DIR / "phenotype" / "prefix_modes.ped"))
        assert ped.phenotype_names == ["Epilepsy", "IQ", "Proband", "Diagnosis"]
        assert ped.phenotype_display_modes["Epilepsy"] == PhenotypeDisplayMode.BINARY
        assert ped.phenotype_display_modes["IQ"] == PhenotypeDisplayMode.GRADIENT
        assert ped.phenotype_display_modes["Proband"] == PhenotypeDisplayMode.OUTLINE
        assert ped.phenotype_display_modes["Diagnosis"] == PhenotypeDisplayMode.CATEGORICAL

    def test_prefix_binary_normalization(self):
        """b_ prefix forces binary normalization (1→0, 2→1)."""
        ped = parse_pedigree(str(TEST_DATA_DIR / "phenotype" / "prefix_modes.ped"))
        fam = list(ped.families.values())[0]
        # father: b_Epilepsy=1 → norm=0 (unaffected)
        assert fam.individuals["father"].extra_phenotypes["Epilepsy"] == 0.0
        # mother: b_Epilepsy=2 → norm=1 (affected)
        assert fam.individuals["mother"].extra_phenotypes["Epilepsy"] == 1.0

    def test_prefix_gradient_normalization(self):
        """g_ prefix forces min-max normalization."""
        ped = parse_pedigree(str(TEST_DATA_DIR / "phenotype" / "prefix_modes.ped"))
        fam = list(ped.families.values())[0]
        # IQ values: 92, 105, 110, 118 → min=92, max=118
        assert fam.individuals["son1"].extra_phenotypes["IQ"] == 0.0  # min
        assert fam.individuals["mother"].extra_phenotypes["IQ"] == 1.0  # max

    def test_prefix_categorical_values(self):
        """c_ prefix stores string category values."""
        ped = parse_pedigree(str(TEST_DATA_DIR / "phenotype" / "prefix_modes.ped"))
        fam = list(ped.families.values())[0]
        assert fam.individuals["father"].extra_phenotypes_categorical["Diagnosis"] == "ASD"
        assert fam.individuals["mother"].extra_phenotypes_categorical["Diagnosis"] == "Control"
        assert fam.individuals["son1"].extra_phenotypes_categorical["Diagnosis"] == "ADHD"

    def test_no_prefix_backward_compatible(self):
        """Files without prefixes parse identically to before."""
        ped = parse_pedigree(str(TEST_DATA_DIR / "phenotype" / "multi_phenotype_pgs.ped"))
        assert "epilepsy" in ped.phenotype_names
        assert len(ped.phenotype_display_modes) == 0  # No prefixes

    def test_synesthesia_prefixes(self):
        """Synesthesia file uses b_ and o_ prefixes."""
        ped = parse_pedigree(str(TEST_DATA_DIR / "complex" / "synesthesia.ped"))
        assert "ASD" in ped.phenotype_names  # b_ stripped
        assert "Synesthesia" in ped.phenotype_names
        assert "Psychiatric comorbidities" in ped.phenotype_names  # o_ stripped
        from pypedigree.models import PhenotypeDisplayMode

        assert ped.phenotype_display_modes["ASD"] == PhenotypeDisplayMode.BINARY
        assert (
            ped.phenotype_display_modes["Psychiatric comorbidities"] == PhenotypeDisplayMode.OUTLINE
        )
