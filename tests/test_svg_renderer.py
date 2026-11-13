"""Tests for the interactive SVG renderer."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

from pypedigree.api import PedigreeBuilder
from pypedigree.config import PedigreeStyle
from pypedigree.models import Pedigree
from pypedigree.svg_renderer import render_interactive_svg


def _parse_svg(svg_str: str) -> ET.Element:
    """Parse SVG string, stripping namespace for easier querying."""
    # Remove namespace for easier xpath
    svg_str = svg_str.replace(' xmlns="http://www.w3.org/2000/svg"', "")
    return ET.fromstring(svg_str)


def _build_trio() -> Pedigree:
    """Build a simple trio pedigree."""
    builder = PedigreeBuilder("FAM1")
    builder.add_individual("father", sex="male", phenotype="affected")
    builder.add_individual("mother", sex="female")
    builder.add_individual("child", sex="male", father="father", mother="mother")
    return builder.build()


class TestSvgBasic:
    """Test basic SVG generation."""

    def test_produces_valid_svg(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        assert svg.startswith("<?xml")
        assert "<svg" in svg
        assert "</svg>" in svg

    def test_parseable_xml(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        root = _parse_svg(svg)
        assert root.tag == "svg"

    def test_empty_pedigree(self) -> None:
        ped = Pedigree()
        svg = render_interactive_svg(ped)
        assert "<svg" in svg


class TestSvgNodes:
    """Test node rendering in SVG."""

    def test_node_elements_present(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        root = _parse_svg(svg)
        nodes = root.findall(".//{http://www.w3.org/2000/svg}g[@class='node']")
        # Also try without namespace
        nodes = [g for g in root.iter("g") if g.get("class") == "node"]
        # At least 3 nodes (father, mother, child) + possible phantoms
        assert len(nodes) >= 3

    def test_node_data_attributes(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        root = _parse_svg(svg)
        nodes = [g for g in root.iter("g") if g.get("class") == "node"]
        for node in nodes:
            assert node.get("data-id") is not None

    def test_male_renders_rect(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        # Father should be a rect
        assert 'data-id="father"' in svg
        # Find the father node and check for rect
        root = _parse_svg(svg)
        father_node = None
        for g in root.iter("g"):
            if g.get("data-id") == "father":
                father_node = g
                break
        assert father_node is not None
        rects = list(father_node.iter("rect"))
        assert len(rects) > 0

    def test_female_renders_circle(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        root = _parse_svg(svg)
        mother_node = None
        for g in root.iter("g"):
            if g.get("data-id") == "mother":
                mother_node = g
                break
        assert mother_node is not None
        circles = list(mother_node.iter("circle"))
        assert len(circles) > 0


class TestSvgTooltips:
    """Test tooltip content."""

    def test_title_elements_present(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        root = _parse_svg(svg)
        titles = list(root.iter("title"))
        assert len(titles) >= 3

    def test_tooltip_contains_id(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        assert "ID: father" in svg
        assert "ID: mother" in svg


class TestSvgInteractivity:
    """Test JavaScript embedding."""

    def test_script_present(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        assert "<script>" in svg
        assert "addEventListener" in svg

    def test_css_styles_present(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        assert ".node.dimmed" in svg
        assert "cursor: pointer" in svg


class TestSvgEdges:
    """Test edge rendering."""

    def test_edges_present(self) -> None:
        ped = _build_trio()
        svg = render_interactive_svg(ped)
        root = _parse_svg(svg)
        edges = [el for el in root.iter("line") if el.get("class") == "edge"]
        # At least some edges (relationship + descent + individual lines)
        assert len(edges) > 0


class TestSvgCarrier:
    """Test carrier rendering in SVG."""

    def test_carrier_clip_path(self) -> None:
        """Legacy notation: half-filled symbol clipped to the shape."""
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("A", sex="male", carrier=True)
        ped = builder.build()
        svg = render_interactive_svg(ped, style=PedigreeStyle(carrier_notation="legacy"))
        assert "clipPath" in svg
        assert 'class="carrier-half-fill"' in svg

    def test_obligate_carrier_dot(self) -> None:
        """Legacy notation: central dot."""
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("B", sex="female", obligate_carrier=True)
        ped = builder.build()
        svg = render_interactive_svg(ped, style=PedigreeStyle(carrier_notation="legacy"))
        # Should have a small filled circle for the dot
        root = _parse_svg(svg)
        b_node = None
        for g in root.iter("g"):
            if g.get("data-id") == "B":
                b_node = g
                break
        assert b_node is not None
        circles = list(b_node.iter("circle"))
        # At least 2 circles: shape + dot
        assert len(circles) >= 2


class TestSvgFileOutput:
    """Test file writing."""

    def test_write_to_file(self, tmp_path: Path) -> None:
        ped = _build_trio()
        output = tmp_path / "test.svg"
        svg = render_interactive_svg(ped, output=str(output))
        assert output.exists()
        assert output.stat().st_size > 0
        content = output.read_text(encoding="utf-8")
        assert content == svg


class TestSvgMultiFamily:
    """Test multi-family rendering."""

    def test_two_families(self) -> None:
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("A", sex="male", family_id="FAM1")
        builder.add_individual("B", sex="female", family_id="FAM2")
        ped = builder.build()
        svg = render_interactive_svg(ped)
        assert 'data-id="A"' in svg
        assert 'data-id="B"' in svg


class TestSvgBuilderIntegration:
    """Test PedigreeBuilder.render_interactive_svg()."""

    def test_builder_method(self, tmp_path: Path) -> None:
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("father", sex="male")
        builder.add_individual("mother", sex="female")
        builder.add_individual("child", sex="female", father="father", mother="mother")
        output = tmp_path / "builder.svg"
        svg = builder.render_interactive_svg(output=str(output))
        assert output.exists()
        assert "<svg" in svg


class TestSvgCli:
    """Test CLI --interactive flag."""

    def test_interactive_cli(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from pypedigree.cli import main

        trio_file = Path(__file__).parent / "data" / "basic" / "simple_trio.ped"
        output = tmp_path / "interactive.svg"
        runner = CliRunner()
        result = runner.invoke(main, [str(trio_file), "-o", str(output), "--interactive"])
        assert result.exit_code == 0, result.output
        assert output.exists()
        content = output.read_text(encoding="utf-8")
        assert "<svg" in content
        assert "<script>" in content


class TestSvgPhenotypes:
    """The SVG draws the same phenotype encoding as the matplotlib renderer."""

    MULTI_PHENO = (
        "FID\tID\tFather\tMother\tSex\tb_ASD\tIQ\tc_Dx\to_Proband\tPGS_x\ta_Variant\n"
        "F\t1\t0\t0\t1\t2\t90\tASD\t2\t0.4\tDMXL2\n"
        "F\t2\t0\t0\t2\t1\t118\tControl\t1\t-0.2\t\n"
        "F\t3\t1\t2\t1\t-9\t92\tADHD\t1\t0.9\tCAMK2B\n"
    )

    def _pedigree(self) -> Pedigree:
        from pypedigree.parser import parse_pedigree_from_string

        return parse_pedigree_from_string(self.MULTI_PHENO)

    def _svg(self, **kwargs: object) -> str:
        return render_interactive_svg(self._pedigree(), **kwargs)  # type: ignore[arg-type]

    def test_every_individual_gets_a_tooltip(self) -> None:
        assert self._svg().count("<title>") == 3

    def test_tooltips_carry_numeric_and_categorical_values(self) -> None:
        svg = self._svg()
        assert "IQ: 90" in svg
        assert "Dx: ASD" in svg
        assert "PGS_x" in svg

    def test_bands_match_the_matplotlib_decisions(self) -> None:
        """Every band is filled with the color `band_cell` gives the PNG renderer."""
        from pypedigree.config import resolve_style
        from pypedigree.phenotypes import band_cell, resolve_phenotype_render, to_hex

        pedigree = self._pedigree()
        names = pedigree.phenotype_names
        style = resolve_style(PedigreeStyle(), names)
        render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, style)
        root = _parse_svg(self._svg())
        for g in root.iter("g"):
            sid = g.get("data-id")
            if sid is None:
                continue
            ind = pedigree.families["F"].individuals[sid]
            fills = [r.get("fill") for r in g.iter("rect") if r.get("class") == "band"]
            expected = [
                to_hex(band_cell(ind, n, style, render, pedigree.phenotype_display_modes).color)
                for n in render.band_names(names)
            ]
            assert fills == expected, sid

    def test_style_phenotype_colors_are_used(self) -> None:
        style = PedigreeStyle()
        style.phenotype_colors = {"ASD": "#123456"}
        assert 'fill="#123456"' in self._svg(style=style)

    def test_na_band_is_written(self) -> None:
        root = _parse_svg(self._svg())
        texts = [el.text for el in root.iter("text")]
        assert "NA" in texts

    def test_outline_column_colors_the_border(self) -> None:
        style = PedigreeStyle()
        style.phenotype_colors = {"Proband": "#ff8800"}
        root = _parse_svg(self._svg(style=style))
        borders = {
            g.get("data-id"): next(e for e in g if e.get("class") == "shape").get("stroke")
            for g in root.iter("g")
            if g.get("data-id")
        }
        assert borders["1"] == "#ff8800"
        assert borders["2"] != "#ff8800"

    def test_pgs_curve_beside_each_scored_individual(self) -> None:
        root = _parse_svg(self._svg())
        curves = [g for g in root.iter("g") if g.get("class") == "pgs"]
        assert len(curves) == 3
        assert all(c.find("polyline") is not None for c in curves)

    def test_annotations_are_drawn_under_the_label(self) -> None:
        root = _parse_svg(self._svg())
        notes = [el.text for el in root.iter("text") if el.get("class") == "annotation"]
        assert notes == ["DMXL2", "CAMK2B"]

    def test_legend_lists_every_column(self) -> None:
        root = _parse_svg(self._svg())
        legend = next(g for g in root.iter("g") if g.get("class") == "legend")
        texts = {el.text for el in legend.iter("text")}
        assert {"Legend", "ASD", "IQ", "Dx"} <= texts
        assert "Outline: Proband (colored border when affected)" in texts
        assert any(t and t.startswith("PGS distribution") for t in texts)

    def test_legend_sits_below_the_pedigree(self) -> None:
        root = _parse_svg(self._svg())
        legend = next(g for g in root.iter("g") if g.get("class") == "legend")
        legend_top = min(float(r.get("y") or 0) for r in legend.iter("rect"))
        lowest_node = max(
            float(e.get("cy") or e.get("y") or 0)
            for g in root.iter("g")
            if g.get("data-id")
            for e in g
            if e.get("class") == "shape"
        )
        assert legend_top > lowest_node

    def test_legend_can_be_hidden(self) -> None:
        assert 'class="legend"' not in self._svg(show_legend=False)

    def test_no_legend_without_phenotypes_or_carriers(self) -> None:
        assert 'class="legend"' not in render_interactive_svg(_build_trio())

    def test_text_placement_does_not_rely_on_dominant_baseline(self) -> None:
        """librsvg and older Inkscape ignore it, which misplaced captions."""
        assert "dominant-baseline" not in self._svg()


class TestSvgFixtures:
    """Every fixture renders to a well-formed SVG with unique element ids."""

    def test_all_fixtures(self) -> None:
        from collections import Counter

        from conftest import TEST_DATA_DIR

        from pypedigree.parser import parse_pedigree, parse_plink_ped

        files = [f for f in sorted(TEST_DATA_DIR.rglob("*.ped")) if "malformed" not in f.name]
        assert files
        for path in files:
            parse = parse_plink_ped if path.parent.name == "plink" else parse_pedigree
            root = _parse_svg(render_interactive_svg(parse(path)))
            ids = Counter(el.get("id") for el in root.iter() if el.get("id"))
            assert not [i for i, n in ids.items() if n > 1], path


class TestSvgCarrierNsgc2022:
    """Default notation: patterned sections defined once as SVG patterns."""

    def test_carrier_uses_a_pattern(self) -> None:
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("A", sex="male", carrier=True)
        svg = render_interactive_svg(builder.build())
        assert '<pattern id="carrier-diagonal-dense"' in svg
        assert 'fill="url(#carrier-diagonal-dense)"' in svg
        assert "carrier-half-fill" not in svg, "no legacy half-fill"

    def test_obligate_carrier_uses_the_sparse_pattern(self) -> None:
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("B", sex="female", obligate_carrier=True)
        svg = render_interactive_svg(builder.build())
        assert 'fill="url(#carrier-diagonal-sparse)"' in svg

    def test_affected_carrier_keeps_its_affected_section(self) -> None:
        builder = PedigreeBuilder("FAM1")
        builder.add_individual("C", sex="female", phenotype="affected", carrier=True)
        root = _parse_svg(render_interactive_svg(builder.build()))
        node = next(g for g in root.iter("g") if g.get("data-id") == "C")
        fills = [r.get("fill") for r in node.iter("rect")]
        assert "black" in fills
        assert "url(#carrier-diagonal-dense)" in fills


class TestSvgOrientation:
    def test_founders_are_drawn_above_their_descendants(self) -> None:
        """SVG's y axis points down: generation 0 must have the smallest y."""
        from conftest import TEST_DATA_DIR

        from pypedigree.parser import parse_pedigree

        pedigree = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
        root = _parse_svg(render_interactive_svg(pedigree))
        y: dict[str, float] = {}
        for g in root.iter("g"):
            sid = g.get("data-id")
            if sid is None:
                continue
            shape = next(e for e in g if e.get("class") == "shape")
            y[sid] = float(shape.get("cy") or shape.get("y"))
        assert y["gfather"] < y["father"] < y["child"]


class TestSvgMultiFamilyLayout:
    def _svg_root(self) -> ET.Element:
        from conftest import TEST_DATA_DIR

        from pypedigree.parser import parse_pedigree

        pedigree = parse_pedigree(TEST_DATA_DIR / "basic/multi_family.ped")
        return _parse_svg(render_interactive_svg(pedigree))

    def test_families_do_not_overlap(self) -> None:
        """Each family used to be drawn at the same coordinates as the others."""
        spans: dict[str, list[float]] = {}
        for g in self._svg_root().iter("g"):
            sid = g.get("data-id")
            if sid is None:
                continue
            shape = next(e for e in g if e.get("class") == "shape")
            spans.setdefault(sid[:2], []).append(float(shape.get("cx") or shape.get("x") or 0))
        f1, f2 = spans["f1"], spans["f2"]
        assert max(f1) < min(f2)

    def test_families_are_titled(self) -> None:
        titles = [
            el.text for el in self._svg_root().iter("text") if el.get("class") == "family-title"
        ]
        assert titles == ["Family: FAM1", "Family: FAM2"]


class TestSvgCliLegend:
    def test_no_legend_flag(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from pypedigree.cli import main

        ped = Path(__file__).parent / "data" / "phenotype" / "prefix_modes.ped"
        for flags, expected in (([], True), (["--no-legend"], False)):
            out = tmp_path / f"legend{len(flags)}.svg"
            result = CliRunner().invoke(main, [str(ped), "-o", str(out), "--interactive", *flags])
            assert result.exit_code == 0, result.output
            assert ('class="legend"' in out.read_text(encoding="utf-8")) is expected


def test_circular_with_interactive_warns(tmp_path: Path) -> None:
    from click.testing import CliRunner

    from pypedigree.cli import main

    trio = Path(__file__).parent / "data" / "basic" / "simple_trio.ped"
    out = tmp_path / "c.svg"
    result = CliRunner().invoke(main, [str(trio), "-o", str(out), "--interactive", "--circular"])
    assert result.exit_code == 0, result.output
    assert "--circular is not supported with --interactive" in result.output
    assert out.exists()


class TestSvgSize:
    """The default pixel size keeps the drawing's proportions."""

    @staticmethod
    def _sizes(ped: Pedigree) -> tuple[float, float, float, float]:
        root = _parse_svg(render_interactive_svg(ped))
        _, _, vb_w, vb_h = (float(v) for v in root.attrib["viewBox"].split())
        return float(root.attrib["width"]), float(root.attrib["height"]), vb_w, vb_h

    def test_wide_pedigree_is_not_letterboxed(self) -> None:
        builder = PedigreeBuilder("F")
        builder.add_individual("f", sex="male")
        builder.add_individual("m", sex="female")
        for i in range(40):
            builder.add_individual(f"c{i}", father="f", mother="m")
        width, height, vb_w, vb_h = self._sizes(builder.build())
        assert width <= 1600 and height <= 1200
        assert abs(width / height - vb_w / vb_h) < 0.02 * vb_w / vb_h

    def test_small_pedigree_is_at_least_600_px_wide(self) -> None:
        width, height, vb_w, vb_h = self._sizes(_build_trio())
        assert width == 600
        assert abs(width / height - vb_w / vb_h) < 0.02 * vb_w / vb_h
