"""Tests for PLINK .fam/.ped integration."""

from __future__ import annotations

from pathlib import Path

import pytest

from pypedigree import (
    Phenotype,
    Sex,
    export_fam,
    parse_fam_file,
    parse_pedigree,
    parse_plink_ped,
)

DATA_DIR = Path(__file__).parent / "data"


class TestParseFamFile:
    """Tests for parse_fam_file()."""

    def test_simple_trio(self):
        ped = parse_fam_file(DATA_DIR / "plink/simple_trio.fam")
        fam = ped.get_family("FAM1")
        assert fam is not None
        assert len(fam.individuals) == 3

        father = fam.individuals["father"]
        assert father.sex == Sex.MALE
        assert father.phenotype == Phenotype.UNAFFECTED

        child = fam.individuals["child"]
        assert child.father_id == "father"
        assert child.mother_id == "mother"
        assert child.phenotype == Phenotype.AFFECTED

    def test_multi_family(self):
        ped = parse_fam_file(DATA_DIR / "plink/multi_family.fam")
        assert len(ped.family_ids) == 2
        assert ped.get_family("FAM1") is not None
        assert ped.get_family("FAM2") is not None

        child2 = ped.get_family("FAM2").individuals["child2"]
        assert child2.sex == Sex.FEMALE
        assert child2.phenotype == Phenotype.UNAFFECTED

    def test_fam_ignores_extra_columns(self):
        """A .fam file with extra whitespace or columns should only parse 6."""
        ped = parse_fam_file(DATA_DIR / "plink/simple_trio.fam")
        fam = ped.get_family("FAM1")
        father = fam.individuals["father"]
        # Should not have extended fields parsed
        assert father.date_of_birth is None
        assert father.is_deceased is False
        assert father.is_proband is False

    def test_fam_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_fam_file("/nonexistent/path.fam")

    def test_fam_unions_inferred(self):
        ped = parse_fam_file(DATA_DIR / "plink/simple_trio.fam")
        fam = ped.get_family("FAM1")
        assert len(fam.unions) == 1
        union = fam.unions[0]
        assert "child" in union.children_ids

    def test_existing_ped_as_fam(self):
        """Standard .ped files (6 columns) should parse fine as .fam."""
        ped = parse_fam_file(DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3


class TestParsePlinkPed:
    """Tests for parse_plink_ped() — ignores genotype columns."""

    def test_genotype_columns_ignored(self):
        ped = parse_plink_ped(DATA_DIR / "plink/plink_with_genotypes.ped")
        fam = ped.get_family("FAM1")
        assert fam is not None
        assert len(fam.individuals) == 3

        father = fam.individuals["father"]
        assert father.sex == Sex.MALE
        assert father.phenotype == Phenotype.UNAFFECTED
        # Genotype data should NOT be parsed as DOB/twins/etc
        assert father.date_of_birth is None
        assert father.is_deceased is False
        assert len(father.twin_info) == 0

    def test_plink_ped_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_plink_ped("/nonexistent/path.ped")

    def test_standard_ped_with_plink_flag(self):
        """Using parse_plink_ped on a standard 6-col .ped should work."""
        ped = parse_plink_ped(DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        assert len(fam.individuals) == 3

    def test_extended_ped_truncated(self):
        """Extended .ped (10 cols) parsed as PLINK should ignore cols 7+."""
        ped = parse_plink_ped(DATA_DIR / "extended/extended_full.ped")
        fam = list(ped.families.values())[0]
        for ind in fam.individuals.values():
            if ind.is_phantom:
                continue
            # Extended fields should NOT be parsed
            assert ind.date_of_birth is None
            assert ind.is_deceased is False
            assert ind.is_proband is False


class TestExportFam:
    """Tests for export_fam()."""

    def test_roundtrip(self, tmp_path):
        """Parse a .fam, export, re-parse — should be identical."""
        ped = parse_fam_file(DATA_DIR / "plink/simple_trio.fam")
        out_path = tmp_path / "exported.fam"
        export_fam(ped, out_path)

        ped2 = parse_fam_file(out_path)
        fam = ped2.get_family("FAM1")
        assert len(fam.individuals) == 3
        assert fam.individuals["child"].father_id == "father"
        assert fam.individuals["child"].phenotype == Phenotype.AFFECTED

    def test_export_excludes_phantoms(self, tmp_path):
        """Phantom parents should not appear in exported .fam."""
        ped = parse_pedigree(DATA_DIR / "basic/single_parent.ped")
        out_path = tmp_path / "no_phantoms.fam"
        export_fam(ped, out_path)

        content = out_path.read_text()
        assert "_phantom_" not in content

    def test_export_excludes_duplicates(self, tmp_path):
        """Duplicate nodes should not appear in exported .fam."""
        ped = parse_pedigree(DATA_DIR / "multi_partner/multiple_partners.ped")
        out_path = tmp_path / "no_dups.fam"
        export_fam(ped, out_path)

        content = out_path.read_text()
        assert "_dup_" not in content

    def test_multi_family_export(self, tmp_path):
        ped = parse_fam_file(DATA_DIR / "plink/multi_family.fam")
        out_path = tmp_path / "multi.fam"
        export_fam(ped, out_path)

        ped2 = parse_fam_file(out_path)
        assert len(ped2.family_ids) == 2

    def test_export_format(self, tmp_path):
        """Verify output is tab-separated with 6 columns."""
        ped = parse_fam_file(DATA_DIR / "plink/simple_trio.fam")
        out_path = tmp_path / "format_check.fam"
        export_fam(ped, out_path)

        lines = out_path.read_text().strip().splitlines()
        for line in lines:
            fields = line.split("\t")
            assert len(fields) == 6, f"Expected 6 columns, got {len(fields)}: {line}"
