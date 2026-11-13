"""Tests for carrier notation: NSGC 2022 patterned sections and the legacy half-fill / dot."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from pypedigree.api import PedigreeBuilder
from pypedigree.models import Individual, NodeLayout, Phenotype, Sex
from pypedigree.parser import parse_pedigree


class TestCarrierModel:
    """Test carrier fields on data models."""

    def test_individual_defaults(self) -> None:
        ind = Individual(sample_id="A", family_id="F1")
        assert ind.is_carrier is False
        assert ind.is_obligate_carrier is False

    def test_individual_carrier(self) -> None:
        ind = Individual(sample_id="A", family_id="F1", is_carrier=True)
        assert ind.is_carrier is True
        assert ind.is_obligate_carrier is False

    def test_individual_obligate_carrier(self) -> None:
        ind = Individual(sample_id="A", family_id="F1", is_obligate_carrier=True)
        assert ind.is_carrier is False
        assert ind.is_obligate_carrier is True

    def test_node_layout_defaults(self) -> None:
        node = NodeLayout(individual_id="A", x=0, y=0, sex=Sex.MALE, phenotype=Phenotype.UNKNOWN)
        assert node.is_carrier is False
        assert node.is_obligate_carrier is False

    def test_node_layout_carrier(self) -> None:
        node = NodeLayout(
            individual_id="A",
            x=0,
            y=0,
            sex=Sex.FEMALE,
            phenotype=Phenotype.UNAFFECTED,
            is_carrier=True,
        )
        assert node.is_carrier is True


class TestCarrierParsing:
    """Test parsing carrier column from .ped files."""

    def test_parse_carrier_file(self) -> None:
        carrier_file = Path(__file__).parent / "data" / "phenotype" / "carrier.ped"
        pedigree = parse_pedigree(str(carrier_file))
        fam = pedigree.families["FAM1"]

        # grandmother: carrier=1
        assert fam.individuals["grandmother"].is_carrier is True
        assert fam.individuals["grandmother"].is_obligate_carrier is False

        # father: carrier=1
        assert fam.individuals["father"].is_carrier is True

        # aunt: carrier=2 (obligate carrier)
        assert fam.individuals["aunt"].is_obligate_carrier is True
        assert fam.individuals["aunt"].is_carrier is False

        # grandfather: carrier=0
        assert fam.individuals["grandfather"].is_carrier is False
        assert fam.individuals["grandfather"].is_obligate_carrier is False

        # child1: carrier=0
        assert fam.individuals["child1"].is_carrier is False

        # child2: carrier=1
        assert fam.individuals["child2"].is_carrier is True

        # cousin1: carrier=2
        assert fam.individuals["cousin1"].is_obligate_carrier is True

    def test_carrier_column_not_in_phenotype_names(self) -> None:
        """The 'carrier' column should not appear in phenotype_names."""
        carrier_file = Path(__file__).parent / "data" / "phenotype" / "carrier.ped"
        pedigree = parse_pedigree(str(carrier_file))
        assert "carrier" not in pedigree.phenotype_names


class TestCarrierLayout:
    """Test carrier flags propagate through layout."""

    def test_carrier_in_node_layout(self) -> None:
        from pypedigree.config import PedigreeStyle
        from pypedigree.layout import compute_layout
        from pypedigree.models import Family

        fam = Family(family_id="F1")
        fam.individuals["A"] = Individual(
            sample_id="A",
            family_id="F1",
            sex=Sex.MALE,
            is_carrier=True,
        )
        fam.individuals["B"] = Individual(
            sample_id="B",
            family_id="F1",
            sex=Sex.FEMALE,
            is_obligate_carrier=True,
        )
        layout = compute_layout(fam, PedigreeStyle())

        node_a = next(n for n in layout.nodes if n.individual_id == "A")
        node_b = next(n for n in layout.nodes if n.individual_id == "B")
        assert node_a.is_carrier is True
        assert node_b.is_obligate_carrier is True


class TestCarrierRendering:
    """Smoke test that carrier rendering doesn't crash."""

    def test_render_carrier_pedigree(self) -> None:
        carrier_file = Path(__file__).parent / "data" / "phenotype" / "carrier.ped"
        pedigree = parse_pedigree(str(carrier_file))
        from pypedigree.renderer import render_pedigree

        fig = render_pedigree(pedigree)
        assert fig is not None
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_render_carrier_to_file(self, tmp_path: Path) -> None:
        carrier_file = Path(__file__).parent / "data" / "phenotype" / "carrier.ped"
        pedigree = parse_pedigree(str(carrier_file))
        from pypedigree.renderer import render_pedigree

        output = tmp_path / "carrier.png"
        fig = render_pedigree(pedigree, output=str(output))
        assert output.exists()
        assert output.stat().st_size > 0
        import matplotlib.pyplot as plt

        plt.close(fig)


class TestCarrierAPI:
    """Test carrier support in PedigreeBuilder."""

    def test_builder_carrier(self) -> None:
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("A", sex="male", carrier=True)
        builder.add_individual("B", sex="female", obligate_carrier=True)
        ped = builder.build()
        fam = ped.families["FAM1"]
        assert fam.individuals["A"].is_carrier is True
        assert fam.individuals["B"].is_obligate_carrier is True

    def test_builder_render_carrier(self, tmp_path: Path) -> None:
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("father", sex="male", carrier=True)
        builder.add_individual("mother", sex="female")
        builder.add_individual(
            "child", sex="male", father="father", mother="mother", obligate_carrier=True
        )
        output = tmp_path / "api_carrier.png"
        fig = builder.render(output=str(output))
        assert output.exists()
        import matplotlib.pyplot as plt

        plt.close(fig)


# ---------------------------------------------------------------------------
# NSGC 2022 notation (default) and cr_ carrier columns
# ---------------------------------------------------------------------------

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402
from click.testing import CliRunner  # noqa: E402
from conftest import TEST_DATA_DIR  # noqa: E402

from pypedigree import PedigreeStyle, render_pedigree  # noqa: E402
from pypedigree.carriers import (  # noqa: E402
    LEGACY_CARRIER_KEY,
    carrier_state,
    resolve_carrier_scheme,
)
from pypedigree.errors import PedigreeError  # noqa: E402
from pypedigree.models import PhenotypeDisplayMode  # noqa: E402
from pypedigree.parser import parse_pedigree_from_string  # noqa: E402

MULTI = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tcr_CFTR\tcr_HEXA\n"
    "F\tgf\t0\t0\t1\t1\t1\t0\n"
    "F\tgm\t0\t0\t2\t1\t0\t1\n"
    "F\tfa\tgf\tgm\t1\t1\t1\t1\n"
    "F\tmo\t0\t0\t2\t2\t0\t0\n"
    "F\tc1\tfa\tmo\t1\t2\t2\t0\n"
    "F\tc2\tfa\tmo\t2\t1\t0\t-9\n"
)


def _hatched(ax):
    return [p for p in ax.patches if p.get_hatch()]


class TestCarrierColumns:
    def test_cr_prefix_declares_a_carrier_column(self) -> None:
        ped = parse_pedigree_from_string(MULTI)
        assert ped.phenotype_names == ["CFTR", "HEXA"]
        assert ped.phenotype_display_modes["CFTR"] == PhenotypeDisplayMode.CARRIER

    def test_carrier_codes_are_kept(self) -> None:
        fam = parse_pedigree_from_string(MULTI).get_family("F")
        assert carrier_state(fam.individuals["fa"], "CFTR") == 1
        assert carrier_state(fam.individuals["c1"], "CFTR") == 2
        assert (
            carrier_state(fam.individuals["c3" if "c3" in fam.individuals else "mo"], "CFTR") == 0
        )
        assert carrier_state(fam.individuals["c2"], "HEXA") is None

    def test_scheme_orders_sources_and_finds_obligates(self) -> None:
        ped = parse_pedigree(TEST_DATA_DIR / "phenotype/carrier.ped")
        people = list(ped.get_family("FAM1").individuals.values())
        scheme = resolve_carrier_scheme(people, ["CFTR"], "nsgc2022")
        assert scheme.keys == (LEGACY_CARRIER_KEY, "CFTR")
        assert LEGACY_CARRIER_KEY in scheme.obligate_keys
        legacy = resolve_carrier_scheme(people, [], "legacy")
        assert legacy.keys == ()
        assert legacy.has_legacy_plain and legacy.has_legacy_obligate


class TestNsgc2022Rendering:
    def test_default_notation_is_nsgc2022(self) -> None:
        assert PedigreeStyle().carrier_notation == "nsgc2022"

    def test_unknown_notation_is_rejected(self) -> None:
        ped = parse_pedigree(TEST_DATA_DIR / "phenotype/carrier.ped")
        with pytest.raises(PedigreeError, match="carrier_notation"):
            render_pedigree(ped, style=PedigreeStyle(carrier_notation="1995"))

    def test_every_carrier_gets_a_patterned_section(self) -> None:
        """carrier.ped has 3 carriers and 2 obligate carriers."""
        fig = render_pedigree(parse_pedigree(TEST_DATA_DIR / "phenotype/carrier.ped"))
        try:
            hatches = [p.get_hatch() for p in _hatched(fig.axes[0])]
            assert sorted(hatches) == sorted(["////"] * 3 + ["//"] * 2)
        finally:
            plt.close(fig)

    def test_legacy_notation_draws_no_pattern(self) -> None:
        style = PedigreeStyle(carrier_notation="legacy")
        fig = render_pedigree(parse_pedigree(TEST_DATA_DIR / "phenotype/carrier.ped"), style=style)
        try:
            assert _hatched(fig.axes[0]) == []
        finally:
            plt.close(fig)

    def test_affected_carrier_keeps_the_affected_section(self) -> None:
        """NSGC 2022: an ATM carrier with breast cancer shows both."""
        builder = PedigreeBuilder("F")
        builder.add_individual("x", sex="female", phenotype="affected", carrier=True)
        fig = builder.render()
        try:
            ax = fig.axes[0]
            assert len(_hatched(ax)) == 1
            filled = [
                p for p in ax.patches if not p.get_hatch() and p.get_facecolor()[:3] == (0, 0, 0)
            ]
            assert filled, "the affected section is filled"
        finally:
            plt.close(fig)

    def test_two_carrier_columns_get_two_patterns(self) -> None:
        fig = render_pedigree(parse_pedigree_from_string(MULTI))
        try:
            hatches = {p.get_hatch() for p in _hatched(fig.axes[0])}
            assert {"////", "\\\\\\\\", "//"} <= hatches
            labels = {t.get_text() for t in fig.axes[1].texts}
            assert {"Carrier: CFTR", "Obligate carrier: CFTR", "Carrier: HEXA"} <= labels
        finally:
            plt.close(fig)

    def test_carrier_column_is_a_band_next_to_phenotype_bands(self) -> None:
        text = (
            "FID\tIID\tFather\tMother\tSex\tb_Cancer\tcr_ATM\n"
            "F\ta\t0\t0\t1\t1\t1\n"
            "F\tb\t0\t0\t2\t2\t0\n"
        )
        fig = render_pedigree(parse_pedigree_from_string(text))
        try:
            assert len(_hatched(fig.axes[0])) == 1
            labels = {t.get_text() for t in fig.axes[1].texts}
            assert "Carrier: ATM" in labels and "Cancer" in labels
            assert "ATM" not in labels, "a carrier column is not a phenotype row"
        finally:
            plt.close(fig)

    def test_circular_layout_draws_the_sections(self) -> None:
        fig = render_pedigree(
            parse_pedigree(TEST_DATA_DIR / "phenotype/carrier.ped"), circular=True
        )
        try:
            assert len(_hatched(fig.axes[0])) == 5
        finally:
            plt.close(fig)

    def test_cli_option(self, tmp_path: Path) -> None:
        from pypedigree.cli import main

        out = tmp_path / "c.png"
        result = CliRunner().invoke(
            main,
            [
                str(TEST_DATA_DIR / "phenotype/carrier.ped"),
                "-o",
                str(out),
                "--carrier-notation",
                "legacy",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
