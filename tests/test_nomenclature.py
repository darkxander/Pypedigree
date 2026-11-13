"""Tests for the standard notation beyond sex, affection status and death.

Pregnancies and pregnancy losses, adoption, childlessness, assisted
reproduction, gender identity and the consultand (NSGC nomenclature,
Bennett et al. 2008 and 2022): parsing, the shared decisions of
``pypedigree.symbols``, the layout, and both renderers.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest
from matplotlib.axes import Axes
from matplotlib.patches import Polygon
from matplotlib.text import Text

from pypedigree import symbols
from pypedigree.api import PedigreeBuilder
from pypedigree.config import PedigreeStyle
from pypedigree.errors import PedigreeError, PedigreeParseError
from pypedigree.layout import compute_layout, compute_layout_circular
from pypedigree.models import (
    Adoption,
    Childlessness,
    EdgeType,
    Gender,
    Individual,
    LineStyle,
    NodeLayout,
    Notation,
    Pedigree,
    Phenotype,
    Pregnancy,
    ReproductiveRole,
    Sex,
)
from pypedigree.parser import build_notation, parse_pedigree, parse_pedigree_from_string
from pypedigree.renderer import _draw_label, _render_layout, render_pedigree
from pypedigree.svg_renderer import render_interactive_svg

FIXTURE = Path(__file__).parent / "data" / "nomenclature" / "nomenclature.ped"
HALF = PedigreeStyle().shape_size / 2

HEADER = "FID\tIID\tFather\tMother\tSex\tPhenotype\t{}"


def _one(column: str, value: str, sex: str = "0") -> Individual:
    """Parse a one-person file with one notation column."""
    text = HEADER.format(column) + f"\nF\ta\t0\t0\t{sex}\t1\t{value}\n"
    return parse_pedigree_from_string(text).families["F"].individuals["a"]


def _family_with_child(**child: str) -> Pedigree:
    builder = PedigreeBuilder("F")
    builder.add_individual("f", sex="male")
    builder.add_individual("m", sex="female")
    builder.add_individual("c", father="f", mother="m", **child)  # type: ignore[arg-type]
    return builder.build()


def _node(ped: Pedigree, sid: str) -> NodeLayout:
    layout = compute_layout(ped.families["F"])
    return next(n for n in layout.nodes if n.individual_id == sid)


def _svg_root(ped: Pedigree, **kwargs: object) -> ET.Element:
    svg = render_interactive_svg(ped, **kwargs)  # type: ignore[arg-type]
    return ET.fromstring(svg.replace(' xmlns="http://www.w3.org/2000/svg"', ""))


def _svg_node(root: ET.Element, sid: str) -> ET.Element:
    found = root.find(f".//g[@data-id='{sid}']")
    assert found is not None
    return found


class TestParsing:
    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            ("P", Pregnancy.ONGOING),
            ("sab", Pregnancy.MISCARRIAGE),
            ("Miscarriage", Pregnancy.MISCARRIAGE),
            ("TOP", Pregnancy.TERMINATION),
            ("ECT", Pregnancy.ECTOPIC),
            ("SB", Pregnancy.STILLBIRTH),
            ("stillbirth", Pregnancy.STILLBIRTH),
        ],
    )
    def test_pregnancy_codes(self, value: str, expected: Pregnancy) -> None:
        assert _one("Pregnancy", value).notation.pregnancy == expected

    @pytest.mark.parametrize(
        ("column", "value", "attr", "expected"),
        [
            ("Adopted", "in", "adoption", Adoption.IN),
            ("Adoption", "OUT", "adoption", Adoption.OUT),
            ("Adopted", "1", "adoption", Adoption.IN),
            ("Childless", "choice", "childlessness", Childlessness.BY_CHOICE),
            ("Childless", "2", "childlessness", Childlessness.INFERTILITY),
            ("Infertility", "infertile", "childlessness", Childlessness.INFERTILITY),
            ("ART", "D", "reproductive_role", ReproductiveRole.DONOR),
            ("ART", "surrogate", "reproductive_role", ReproductiveRole.SURROGATE),
            ("ART", "GC", "reproductive_role", ReproductiveRole.SURROGATE),
            ("Gender_identity", "woman", "gender", Gender.WOMAN),
            ("Identity", "NB", "gender", Gender.NONBINARY),
            ("gender_identity", "1", "gender", Gender.MAN),
        ],
    )
    def test_column_aliases_and_codes(
        self, column: str, value: str, attr: str, expected: object
    ) -> None:
        assert getattr(_one(column, value).notation, attr) == expected

    @pytest.mark.parametrize("missing", ["0", "-9", "NA", ".", "na"])
    def test_missing_values_note_nothing(self, missing: str) -> None:
        assert _one("Pregnancy", missing).notation == Notation()

    def test_consultand_column(self) -> None:
        assert _one("Consultand", "1").notation.is_consultand
        assert not _one("Consultand", "0").notation.is_consultand

    def test_unknown_value_names_the_line(self) -> None:
        with pytest.raises(PedigreeParseError, match=r"Line 2: unknown pregnancy value 'XX'"):
            _one("Pregnancy", "XX")

    def test_gender_column_is_the_sex_when_alone(self) -> None:
        """A file that calls its sex column "gender" keeps working."""
        text = "FID\tIID\tFather\tMother\tGender\tPhenotype\tASD\nF\ta\t0\t0\t2\t1\t2\n"
        ind = parse_pedigree_from_string(text).families["F"].individuals["a"]
        assert ind.sex == Sex.FEMALE
        assert ind.notation.gender is None

    def test_gender_column_is_the_identity_next_to_a_sex_column(self) -> None:
        text = "FID\tIID\tFather\tMother\tSex\tPhenotype\tGender\nF\ta\t0\t0\t1\t1\twoman\n"
        ind = parse_pedigree_from_string(text).families["F"].individuals["a"]
        assert ind.sex == Sex.MALE
        assert ind.notation.gender == Gender.WOMAN

    def test_notation_columns_are_not_phenotypes(self) -> None:
        ped = parse_pedigree(FIXTURE)
        assert ped.phenotype_names == []
        assert ped.annotation_names == ["Note"]

    def test_fixture(self) -> None:
        fam = parse_pedigree(FIXTURE).families["FAM1"]
        notation = {sid: ind.notation for sid, ind in fam.individuals.items()}
        assert notation["II-1"].childlessness == Childlessness.INFERTILITY
        assert notation["II-2"].is_consultand
        assert notation["II-4"].gender == Gender.WOMAN
        assert notation["III-2"].adoption == Adoption.IN
        assert notation["III-3"].adoption == Adoption.OUT
        assert notation["III-7"].pregnancy == Pregnancy.ONGOING

    def test_duplicate_partner_keeps_notation_but_not_consultand(self) -> None:
        text = (
            HEADER.format("Gender_identity\tConsultand")
            + "\nF\ta\t0\t0\t1\t1\twoman\t1"
            + "\nF\tb\t0\t0\t2\t1\t0\t0"
            + "\nF\tc\t0\t0\t2\t1\t0\t0"
            + "\nF\tk1\ta\tb\t1\t1\t0\t0"
            + "\nF\tk2\ta\tc\t1\t1\t0\t0\n"
        )
        fam = parse_pedigree_from_string(text).families["F"]
        dups = [i for i in fam.individuals.values() if i.duplicate_of == "a"]
        assert dups
        assert all(d.notation.gender == Gender.WOMAN for d in dups)
        assert not any(d.notation.is_consultand for d in dups)


class TestBuildNotation:
    def test_codes_and_enums(self) -> None:
        notation = build_notation(
            pregnancy="top", adopted=Adoption.OUT, childless="choice", art="s", gender="x"
        )
        assert notation == Notation(
            pregnancy=Pregnancy.TERMINATION,
            adoption=Adoption.OUT,
            childlessness=Childlessness.BY_CHOICE,
            reproductive_role=ReproductiveRole.SURROGATE,
            gender=Gender.NONBINARY,
        )

    def test_unknown_code(self) -> None:
        with pytest.raises(PedigreeError, match="Unknown art value"):
            build_notation(art="Z")

    def test_builder(self) -> None:
        ped = _family_with_child(pregnancy="SB", adopted="in")
        child = ped.families["F"].individuals["c"]
        assert child.notation.pregnancy == Pregnancy.STILLBIRTH
        assert child.notation.adoption == Adoption.IN

    def test_builder_rejects_unknown_code(self) -> None:
        builder = PedigreeBuilder("F")
        with pytest.raises(PedigreeError):
            builder.add_individual("a", gender="martian")


class TestSymbols:
    @pytest.mark.parametrize(
        ("sex", "gender", "shape", "note"),
        [
            (Sex.MALE, None, Sex.MALE, None),
            (Sex.MALE, Gender.MAN, Sex.MALE, None),  # cisgender: nothing written
            (Sex.MALE, Gender.WOMAN, Sex.FEMALE, "AMAB"),
            (Sex.FEMALE, Gender.MAN, Sex.MALE, "AFAB"),
            (Sex.FEMALE, Gender.NONBINARY, Sex.UNKNOWN, "AFAB"),
            (Sex.UNKNOWN, Gender.NONBINARY, Sex.UNKNOWN, None),  # sex at birth unknown
        ],
    )
    def test_gender_decides_shape_and_note(
        self, sex: Sex, gender: Gender | None, shape: Sex, note: str | None
    ) -> None:
        notation = Notation(gender=gender)
        assert symbols.symbol_sex(sex, notation) == shape
        assert symbols.sex_at_birth_note(sex, notation) == note

    def test_notes_order(self) -> None:
        ind = Individual(
            "a",
            "F",
            sex=Sex.FEMALE,
            notation=Notation(pregnancy=Pregnancy.STILLBIRTH, gender=Gender.MAN),
        )
        assert symbols.notes(ind) == ["SB", "AFAB"]

    def test_pregnancy_losses_are_triangles(self) -> None:
        losses = {p for p in Pregnancy if symbols.is_pregnancy_loss(Notation(pregnancy=p))}
        assert losses == {Pregnancy.MISCARRIAGE, Pregnancy.TERMINATION, Pregnancy.ECTOPIC}

    def test_stillbirth_is_struck_through(self) -> None:
        ind = Individual("a", "F", notation=Notation(pregnancy=Pregnancy.STILLBIRTH))
        assert symbols.is_struck_through(ind)
        assert not symbols.is_struck_through(Individual("b", "F"))

    def test_inner_letter(self) -> None:
        assert symbols.inner_letter(Notation(pregnancy=Pregnancy.ONGOING)) == "P"
        assert symbols.inner_letter(Notation(reproductive_role=ReproductiveRole.DONOR)) == "D"
        assert symbols.inner_letter(Notation(reproductive_role=ReproductiveRole.SURROGATE)) == "S"
        assert symbols.inner_letter(Notation()) is None

    def test_triangle_apex_meets_the_line_from_the_parents(self) -> None:
        apex, left, right = symbols.triangle_points(0.0, 0.0, HALF)
        assert apex == (0.0, HALF)
        assert left[1] == right[1] < 0
        flipped = symbols.triangle_points(0.0, 0.0, HALF, up=-1.0)
        assert flipped[0] == (0.0, -HALF)

    @pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, Sex.UNKNOWN])
    def test_brackets_enclose_the_symbol(self, sex: Sex) -> None:
        left, right = symbols.bracket_lines(0.0, 0.0, HALF, sex)
        reach = symbols.shape_half_width(sex, HALF)
        assert max(x for x, _ in left) < 0 and left[1][0] < -reach
        assert right[1][0] > reach
        assert left[1][1] > reach - 1e-9 and left[2][1] < -reach + 1e-9

    def test_childless_marker(self) -> None:
        segs = symbols.childless_segments(0.0, 0.0, Sex.MALE, HALF, Childlessness.BY_CHOICE)
        assert len(segs) == 2
        (stem_top, stem_bottom), (bar_l, bar_r) = segs
        assert stem_top == pytest.approx((0.0, -HALF))
        assert stem_bottom == pytest.approx((0.0, -HALF - symbols.CHILDLESS_STEM))
        assert bar_l[1] == pytest.approx(bar_r[1])
        infertile = symbols.childless_segments(0.0, 0.0, Sex.MALE, HALF, Childlessness.INFERTILITY)
        assert len(infertile) == 3
        assert symbols.childless_depth(Notation(childlessness=Childlessness.INFERTILITY)) > (
            symbols.childless_depth(Notation(childlessness=Childlessness.BY_CHOICE))
        )
        assert symbols.childless_depth(Notation()) == 0.0

    def test_childless_marker_follows_its_direction(self) -> None:
        segs = symbols.childless_segments(
            0.0, 0.0, Sex.FEMALE, HALF, Childlessness.BY_CHOICE, direction=0.0
        )
        (start, end), _ = segs
        assert start == pytest.approx((HALF, 0.0))
        assert end[0] > start[0] and end[1] == pytest.approx(0.0)


class TestLayout:
    def test_node_shape_follows_gender(self) -> None:
        node = _node(_family_with_child(sex="male", gender="woman"), "c")
        assert node.sex == Sex.FEMALE
        assert node.notes == ["AMAB"]

    def test_stillbirth_node_is_deceased(self) -> None:
        node = _node(_family_with_child(pregnancy="SB"), "c")
        assert node.is_deceased
        assert node.notes == ["SB"]

    @pytest.mark.parametrize("siblings", [0, 2])
    def test_adopted_in_line_is_dashed(self, siblings: int) -> None:
        builder = PedigreeBuilder("F")
        builder.add_individual("f", sex="male")
        builder.add_individual("m", sex="female")
        builder.add_individual("c", father="f", mother="m", adopted="in")
        for i in range(siblings):
            builder.add_individual(f"s{i}", father="f", mother="m")
        layout = compute_layout(builder.build().families["F"])
        x = next(n.x for n in layout.nodes if n.individual_id == "c")
        lines = [e for e in layout.edges if e.edge_type == EdgeType.INDIVIDUAL]
        assert [e.line_style for e in lines if e.start[0] == x] == [LineStyle.DASHED]
        assert all(e.line_style == LineStyle.SOLID for e in lines if e.start[0] != x)

    def test_adopted_out_line_is_solid(self) -> None:
        layout = compute_layout(_family_with_child(adopted="out").families["F"])
        lines = [e for e in layout.edges if e.edge_type == EdgeType.INDIVIDUAL]
        assert [e.line_style for e in lines] == [LineStyle.SOLID]

    def test_childless_edges_and_room_below(self) -> None:
        ped = _family_with_child(childless="infertility")
        layout = compute_layout(ped.families["F"])
        assert len([e for e in layout.edges if e.edge_type == EdgeType.CHILDLESS]) == 3
        node = next(n for n in layout.nodes if n.individual_id == "c")
        assert node.marker_depth == symbols.childless_depth(node.notation)
        plain = compute_layout(_family_with_child().families["F"])
        gap = abs(node.y)
        assert gap > abs(next(n.y for n in plain.nodes if n.individual_id == "c"))

    def test_circular_childless_marker_points_outward(self) -> None:
        ped = parse_pedigree(FIXTURE)
        layout = compute_layout_circular(ped.families["FAM1"])
        node = next(n for n in layout.nodes if n.individual_id == "II-1")
        assert node.label_angle is not None
        stem = next(e for e in layout.edges if e.edge_type == EdgeType.CHILDLESS)
        direction = math.atan2(stem.end[1] - stem.start[1], stem.end[0] - stem.start[0])
        assert math.cos(direction - node.label_angle) == pytest.approx(1.0)

    def test_existing_pedigrees_unchanged(self) -> None:
        """No notation, no change: nodes get default notation and no notes."""
        layout = compute_layout(_family_with_child().families["F"])
        assert all(n.notes == [] and n.marker_depth == 0.0 for n in layout.nodes)
        assert all(e.edge_type != EdgeType.CHILDLESS for e in layout.edges)


class TestPngRenderer:
    @staticmethod
    def _draw(ped: Pedigree, sid: str = "c") -> tuple[Axes, NodeLayout]:
        style = PedigreeStyle()
        layout = compute_layout(ped.families["F"])
        fig, ax = plt.subplots()
        _render_layout(ax, layout, style, show_labels=True)
        plt.close(fig)
        return ax, next(n for n in layout.nodes if n.individual_id == sid)

    @staticmethod
    def _texts(ax: Axes) -> dict[str, Text]:
        return {t.get_text(): t for t in ax.texts}

    def test_pregnancy_loss_is_a_triangle(self) -> None:
        ax, node = self._draw(_family_with_child(pregnancy="sab", phenotype="affected"))
        triangles = [p for p in ax.patches if isinstance(p, Polygon)]
        assert len(triangles) == 1
        assert triangles[0].get_facecolor()[:3] == (0.0, 0.0, 0.0)

    def test_termination_is_struck_through(self) -> None:
        ax, _ = self._draw(_family_with_child(pregnancy="TOP"))
        ax_plain, _ = self._draw(_family_with_child(pregnancy="SAB"))
        assert len(ax.lines) == len(ax_plain.lines) + 1

    def test_letter_inside(self) -> None:
        ax, node = self._draw(_family_with_child(pregnancy="P"))
        letter = self._texts(ax)["P"]
        assert letter.get_position() == (node.x, node.y)

    def test_letter_white_on_filled_symbol(self) -> None:
        ax, _ = self._draw(_family_with_child(art="D", phenotype="affected"))
        assert self._texts(ax)["D"].get_color() == "white"

    def test_consultand_arrow_has_no_p(self) -> None:
        ax, _ = self._draw(_family_with_child(consultand=True))  # type: ignore[arg-type]
        assert "P" not in self._texts(ax)
        ax_proband, _ = self._draw(_family_with_child(proband=True))  # type: ignore[arg-type]
        assert "P" in self._texts(ax_proband)

    def test_brackets(self) -> None:
        ax, _ = self._draw(_family_with_child(adopted="out"))
        ax_plain, _ = self._draw(_family_with_child())
        assert len(ax.lines) == len(ax_plain.lines) + 2

    def test_label_moves_below_the_childless_marker(self) -> None:
        style = PedigreeStyle()
        fig, ax = plt.subplots()
        for depth in (0.0, 0.4):
            node = NodeLayout("a", 0.0, 0.0, Sex.MALE, Phenotype.UNAFFECTED, label="a")
            node.marker_depth = depth
            _draw_label(ax, node, style)
        plt.close(fig)
        ys = [t.get_position()[1] for t in ax.texts]
        assert ys[0] - ys[1] == pytest.approx(0.4)

    def test_notes_below_label(self) -> None:
        ax, _ = self._draw(_family_with_child(sex="female", gender="nonbinary"))
        texts = self._texts(ax)
        assert texts["AFAB"].get_position()[1] < texts["c"].get_position()[1]

    def test_triangle_takes_first_band_color(self) -> None:
        builder = PedigreeBuilder("F")
        builder.add_individual("f", sex="male", phenotypes={"b_ASD": 1})
        builder.add_individual("m", sex="female", phenotypes={"b_ASD": 1})
        builder.add_individual(
            "c", father="f", mother="m", pregnancy="SAB", phenotypes={"b_ASD": 2}
        )
        style = PedigreeStyle(phenotype_colors={"ASD": "#e74c3c"})
        fig = render_pedigree(builder.build(), style=style)
        ax = fig.axes[0]
        triangles = [p for p in ax.patches if isinstance(p, Polygon) and p.get_closed()]
        colors = {matplotlib.colors.to_hex(p.get_facecolor()) for p in triangles}
        plt.close(fig)
        assert "#e74c3c" in colors

    def test_fixture_renders(self) -> None:
        fig = render_pedigree(parse_pedigree(FIXTURE))
        plt.close(fig)


class TestSvgRenderer:
    def test_pregnancy_loss_triangle_and_termination(self) -> None:
        root = _svg_root(_family_with_child(pregnancy="TOP"))
        node = _svg_node(root, "c")
        assert node.find("polygon[@class='shape']") is not None
        assert node.find("line[@class='termination']") is not None

    def test_letter_brackets_notes(self) -> None:
        ped = _family_with_child(sex="male", gender="woman", adopted="in", pregnancy="P")
        node = _svg_node(_svg_root(ped), "c")
        assert node.find("text[@class='symbol-letter']").text == "P"  # type: ignore[union-attr]
        assert len(node.findall("polyline[@class='adoption-bracket']")) == 2
        assert [t.text for t in node.findall("text[@class='note']")] == ["AMAB"]

    def test_adopted_in_line_dashed(self) -> None:
        svg = render_interactive_svg(_family_with_child(adopted="in"))
        assert "stroke-dasharray" in svg

    def test_consultand(self) -> None:
        node = _svg_node(_svg_root(_family_with_child(consultand=True)), "c")  # type: ignore[arg-type]
        assert node.find("line[@marker-end]") is not None
        assert all(t.text != "P" for t in node.findall("text"))

    def test_tooltip(self) -> None:
        ped = _family_with_child(sex="female", gender="man", childless="choice")
        title = _svg_node(_svg_root(ped), "c").find("title")
        assert title is not None and title.text is not None
        assert "Sex: Female" in title.text
        assert "Gender: Man" in title.text
        assert "No children by choice" in title.text

    def test_childless_edges(self) -> None:
        plain = render_interactive_svg(_family_with_child())
        marked = render_interactive_svg(_family_with_child(childless="infertility"))
        assert marked.count('class="edge"') == plain.count('class="edge"') + 3

    def test_fixture_renders(self) -> None:
        root = _svg_root(parse_pedigree(FIXTURE))
        assert root.findall(".//g[@class='node']")
