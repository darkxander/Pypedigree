"""Tests for pypedigree.api — programmatic builder API."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from pypedigree.api import PedigreeBuilder
from pypedigree.errors import PedigreeError
from pypedigree.models import Phenotype, Sex


class TestPedigreeBuilderBasic:
    """Test basic builder functionality."""

    def test_add_individual(self):
        b = PedigreeBuilder()
        result = b.add_individual("father", sex="male")
        assert result is b  # fluent chaining

    def test_build_simple_trio(self):
        b = PedigreeBuilder("FAM1")
        b.add_individual("father", sex="male", phenotype="unaffected")
        b.add_individual("mother", sex="female", phenotype="unaffected")
        b.add_individual(
            "child", sex="male", father="father", mother="mother", phenotype="affected"
        )

        ped = b.build()
        assert len(ped.families) == 1
        fam = ped.get_family("FAM1")
        assert fam is not None
        assert len(fam.individuals) == 3
        assert len(fam.unions) == 1
        assert fam.unions[0].partner_ids == ("father", "mother")
        assert fam.unions[0].children_ids == ["child"]

    def test_sex_parsing(self):
        b = PedigreeBuilder()
        b.add_individual("a", sex="male")
        b.add_individual("b", sex="female")
        b.add_individual("c", sex="unknown")
        b.add_individual("d", sex="m")
        b.add_individual("e", sex="f")
        b.add_individual("f", sex=Sex.MALE)

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert fam.individuals["a"].sex == Sex.MALE
        assert fam.individuals["b"].sex == Sex.FEMALE
        assert fam.individuals["c"].sex == Sex.UNKNOWN
        assert fam.individuals["d"].sex == Sex.MALE
        assert fam.individuals["e"].sex == Sex.FEMALE
        assert fam.individuals["f"].sex == Sex.MALE

    def test_phenotype_parsing(self):
        b = PedigreeBuilder()
        b.add_individual("a", phenotype="affected")
        b.add_individual("b", phenotype="unaffected")
        b.add_individual("c", phenotype=Phenotype.AFFECTED)

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert fam.individuals["a"].phenotype == Phenotype.AFFECTED
        assert fam.individuals["b"].phenotype == Phenotype.UNAFFECTED
        assert fam.individuals["c"].phenotype == Phenotype.AFFECTED

    def test_duplicate_individual_raises(self):
        b = PedigreeBuilder()
        b.add_individual("a", sex="male")
        with pytest.raises(PedigreeError, match="already added"):
            b.add_individual("a", sex="female")

    def test_deceased_and_proband(self):
        b = PedigreeBuilder()
        b.add_individual("a", deceased=True)
        b.add_individual("b", proband=True)

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert fam.individuals["a"].is_deceased is True
        assert fam.individuals["b"].is_proband is True

    def test_date_of_birth(self):
        b = PedigreeBuilder()
        b.add_individual("a", date_of_birth="19900101")

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert fam.individuals["a"].date_of_birth == "19900101"


class TestPedigreeBuilderUnions:
    """Test union inference and explicit unions."""

    def test_auto_infer_union(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("m", sex="female")
        b.add_individual("c1", sex="male", father="f", mother="m")
        b.add_individual("c2", sex="female", father="f", mother="m")

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert len(fam.unions) == 1
        assert set(fam.unions[0].children_ids) == {"c1", "c2"}

    def test_explicit_consanguineous_union(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("m", sex="female")
        b.add_individual("c", sex="male", father="f", mother="m")
        b.add_union("f", "m", consanguineous=True)

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert fam.unions[0].is_consanguinous is True

    def test_explicit_broken_union(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("m", sex="female")
        b.add_individual("c", sex="male", father="f", mother="m")
        b.add_union("f", "m", broken=True)

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert fam.unions[0].is_broken is True

    def test_childless_union(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("m", sex="female")
        b.add_union("f", "m")

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert len(fam.unions) == 1
        assert fam.unions[0].children_ids == []


class TestPedigreeBuilderPhenotypes:
    """Test extra phenotype handling."""

    def test_extra_phenotypes_normalized(self):
        b = PedigreeBuilder()
        b.add_individual("a", sex="male", phenotypes={"IQ": 80})
        b.add_individual("b", sex="female", phenotypes={"IQ": 120})

        ped = b.build()
        fam = ped.get_family("FAM1")

        # Raw values preserved
        assert fam.individuals["a"].extra_phenotypes_raw["IQ"] == 80.0
        assert fam.individuals["b"].extra_phenotypes_raw["IQ"] == 120.0

        # Normalized to [0, 1]
        assert fam.individuals["a"].extra_phenotypes["IQ"] == pytest.approx(0.0)
        assert fam.individuals["b"].extra_phenotypes["IQ"] == pytest.approx(1.0)

    def test_binary_phenotypes(self):
        b = PedigreeBuilder()
        b.add_individual("a", sex="male", phenotypes={"epilepsy": 1})
        b.add_individual("b", sex="female", phenotypes={"epilepsy": 0})

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert fam.individuals["a"].extra_phenotypes_raw["epilepsy"] == 1.0
        assert fam.individuals["b"].extra_phenotypes_raw["epilepsy"] == 0.0

    def test_phenotype_names_tracked(self):
        b = PedigreeBuilder()
        b.add_individual("a", phenotypes={"epilepsy": 1, "IQ": 100})
        b.add_individual("b", phenotypes={"epilepsy": 0, "IQ": 110})

        ped = b.build()
        assert ped.phenotype_names == ["epilepsy", "IQ"]


class TestPedigreeBuilderMultiFamily:
    """Test multi-family support."""

    def test_two_families(self):
        b = PedigreeBuilder("FAM1")
        b.add_individual("f1", sex="male", family_id="FAM1")
        b.add_individual("m1", sex="female", family_id="FAM1")
        b.add_individual("f2", sex="male", family_id="FAM2")
        b.add_individual("m2", sex="female", family_id="FAM2")

        ped = b.build()
        assert len(ped.families) == 2
        assert "FAM1" in ped.families
        assert "FAM2" in ped.families


class TestPedigreeBuilderPhantoms:
    """Test phantom parent inference."""

    def test_phantom_father_created(self):
        b = PedigreeBuilder()
        b.add_individual("m", sex="female")
        b.add_individual("c", sex="male", father="phantom_f", mother="m")

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert "phantom_f" in fam.individuals
        assert fam.individuals["phantom_f"].is_phantom is True
        assert fam.individuals["phantom_f"].sex == Sex.MALE

    def test_phantom_mother_created(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("c", sex="male", father="f", mother="phantom_m")

        ped = b.build()
        fam = ped.get_family("FAM1")
        assert "phantom_m" in fam.individuals
        assert fam.individuals["phantom_m"].is_phantom is True
        assert fam.individuals["phantom_m"].sex == Sex.FEMALE


class TestPedigreeBuilderRender:
    """Test rendering via builder."""

    def test_render_simple_trio(self):
        b = PedigreeBuilder()
        b.add_individual("father", sex="male", phenotype="unaffected")
        b.add_individual("mother", sex="female", phenotype="unaffected")
        b.add_individual(
            "child", sex="male", father="father", mother="mother", phenotype="affected"
        )

        fig = b.render()
        assert fig is not None
        plt.close(fig)

    def test_render_to_file(self, tmp_path):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("m", sex="female")
        b.add_individual("c", sex="male", father="f", mother="m")

        out = tmp_path / "test.png"
        fig = b.render(output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)

    def test_render_with_phenotypes(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male", phenotypes={"epilepsy": 1, "IQ": 100})
        b.add_individual("m", sex="female", phenotypes={"epilepsy": 0, "IQ": 120})
        b.add_individual(
            "c", sex="male", father="f", mother="m", phenotypes={"epilepsy": 1, "IQ": 90}
        )

        fig = b.render()
        assert fig is not None
        plt.close(fig)

    def test_render_no_labels(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("m", sex="female")
        b.add_individual("c", sex="male", father="f", mother="m")

        fig = b.render(show_labels=False)
        ax = fig.axes[0]
        # No ID labels should be present
        label_texts = [t.get_text() for t in ax.texts]
        assert "f" not in label_texts
        assert "m" not in label_texts
        assert "c" not in label_texts
        plt.close(fig)

    def test_render_no_legend(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male", phenotypes={"epilepsy": 1})
        b.add_individual("m", sex="female", phenotypes={"epilepsy": 0})
        b.add_individual("c", sex="male", father="f", mother="m", phenotypes={"epilepsy": 1})

        fig = b.render(show_legend=False)
        # Only 1 axes (no legend axes)
        assert len(fig.axes) == 1
        plt.close(fig)

    def test_repr_png(self):
        b = PedigreeBuilder()
        b.add_individual("f", sex="male")
        b.add_individual("m", sex="female")
        b.add_individual("c", sex="male", father="f", mother="m")

        png_bytes = b._repr_png_()
        assert isinstance(png_bytes, bytes)
        assert len(png_bytes) > 100
        # PNG magic bytes
        assert png_bytes[:4] == b"\x89PNG"

    def test_method_chaining(self):
        fig = (
            PedigreeBuilder("FAM1")
            .add_individual("f", sex="male")
            .add_individual("m", sex="female")
            .add_individual("c", sex="male", father="f", mother="m")
            .render()
        )
        assert fig is not None
        plt.close(fig)
