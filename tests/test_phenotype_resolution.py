"""Tests for per-figure phenotype resolution and style resolution.

These cover the guarantees that keep bands, legend and spacing consistent:
style resolution is pure and idempotent, and every phenotype display
decision is taken once per figure rather than accumulated in module state.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest
from conftest import PHENOTYPE_DIR

from pypedigree import PedigreeStyle, parse_pedigree, parse_pedigree_from_string, render_pedigree
from pypedigree.config import apply_pgs_spacing, is_pgs_column, resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    diverging_color_at,
    diverging_color_for_value,
    normalize_to_unit,
    resolve_phenotype_render,
    to_rgb,
)

PREFIX_MODES = PHENOTYPE_DIR / "prefix_modes.ped"
PGS_PED = PHENOTYPE_DIR / "multi_phenotype_pgs.ped"


def _render_context(pedigree, style=None):
    """Resolve style + phenotype display the same way render_pedigree does."""
    resolved = resolve_style(style, pedigree.phenotype_names)
    return resolved, resolve_phenotype_render(
        pedigree, pedigree.phenotype_names, pedigree.phenotype_display_modes, resolved
    )


class TestIsPgsColumn:
    @pytest.mark.parametrize("name", ["PGS", "pgs", "PGS_autism", "pgs_IQ"])
    def test_recognised(self, name):
        assert is_pgs_column(name)

    @pytest.mark.parametrize("name", ["IQ", "SRS", "epilepsy", "myPGS"])
    def test_not_recognised(self, name):
        assert not is_pgs_column(name)


class TestResolveStyle:
    def test_none_gives_defaults(self):
        assert resolve_style(None).to_dict() == PedigreeStyle().to_dict()

    def test_assigns_palette_colors(self):
        resolved = resolve_style(None, ["epilepsy", "IQ"])
        assert set(resolved.phenotype_colors) == {"epilepsy", "IQ"}

    def test_skips_pgs_columns(self):
        resolved = resolve_style(None, ["IQ", "PGS_autism"])
        assert "PGS_autism" not in resolved.phenotype_colors
        assert "IQ" in resolved.phenotype_colors

    def test_user_colors_win(self):
        style = PedigreeStyle()
        style.phenotype_colors = {"IQ": "#123456"}
        resolved = resolve_style(style, ["IQ", "epilepsy"])
        assert resolved.phenotype_colors == {"IQ": "#123456"}

    def test_does_not_mutate_input(self):
        style = PedigreeStyle()
        resolve_style(style, ["IQ", "PGS"])
        assert style.phenotype_colors == {}
        assert style.couple_gap == PedigreeStyle().couple_gap

    def test_is_idempotent(self):
        once = resolve_style(None, ["IQ", "PGS"])
        twice = resolve_style(once, ["IQ", "PGS"])
        assert once.to_dict() == twice.to_dict()

    def test_copies_nested_containers(self):
        style = PedigreeStyle()
        style.gradient_range = {"IQ": [40, 160]}
        resolved = resolve_style(style, ["IQ"])
        resolved.gradient_range["IQ"][0] = 0
        assert style.gradient_range["IQ"] == [40, 160]


class TestApplyPgsSpacing:
    def test_widens_once(self):
        base = PedigreeStyle()
        widened = apply_pgs_spacing(base, ["PGS"])
        expected = base.couple_gap + base.shape_size * 0.55
        assert widened.couple_gap == pytest.approx(expected)
        assert widened.sibling_gap == pytest.approx(base.sibling_gap + base.shape_size * 0.55)

    def test_no_pgs_returns_input_unchanged(self):
        base = PedigreeStyle()
        assert apply_pgs_spacing(base, ["IQ"]) is base
        assert apply_pgs_spacing(base, None) is base

    def test_prefixed_pgs_columns_widen(self):
        base = PedigreeStyle()
        assert apply_pgs_spacing(base, ["PGS_autism"]).couple_gap > base.couple_gap

    def test_render_applies_spacing_only_once(self):
        """A caller that pre-resolves must not get double spacing."""
        pedigree = parse_pedigree(PGS_PED)
        pre_resolved = resolve_style(None, pedigree.phenotype_names)
        expected = pre_resolved.couple_gap + pre_resolved.shape_size * 0.55
        applied = apply_pgs_spacing(pre_resolved, pedigree.phenotype_names)
        assert applied.couple_gap == pytest.approx(expected)


class TestEffectiveModeResolution:
    def test_prefixes_map_to_modes(self):
        pedigree = parse_pedigree(PREFIX_MODES)
        _, ctx = _render_context(pedigree)
        assert ctx.effective("Epilepsy") == "binary"
        assert ctx.effective("Proband") == "outline"
        assert ctx.effective("Diagnosis") == "categorical"

    def test_auto_continuous_column_is_gradient(self):
        """An unprefixed continuous column must not be treated as binary."""
        pedigree = parse_pedigree(PREFIX_MODES)
        _, ctx = _render_context(pedigree)
        assert ctx.effective("IQ") == "gradient"

    def test_auto_binary_column_is_binary(self):
        text = "FID\tID\tFather\tMother\tSex\tflag\nF\ta\t0\t0\t1\t2\nF\tb\t0\t0\t2\t1\n"
        pedigree = parse_pedigree_from_string(text)
        _, ctx = _render_context(pedigree)
        assert ctx.effective("flag") == "binary"

    def test_auto_mode_is_decided_per_column_not_per_individual(self):
        """A column with both 0/1/2-like and out-of-range values is one mode."""
        text = "FID\tID\tFather\tMother\tSex\tseizures\nF\ta\t0\t0\t1\t1\nF\tb\t0\t0\t2\t17\n"
        pedigree = parse_pedigree_from_string(text)
        _, ctx = _render_context(pedigree)
        assert ctx.effective("seizures") == "gradient"

    def test_pgs_columns_are_not_bands(self):
        pedigree = parse_pedigree(PGS_PED)
        _, ctx = _render_context(pedigree)
        assert ctx.effective("PGS") == "pgs"

    def test_explicit_gradient_range_implies_gradient(self):
        text = "FID\tID\tFather\tMother\tSex\tSRS\nF\ta\t0\t0\t1\t1\nF\tb\t0\t0\t2\t2\n"
        pedigree = parse_pedigree_from_string(text)
        style = PedigreeStyle()
        style.gradient_range = {"SRS": [0, 10]}
        _, ctx = _render_context(pedigree, style)
        assert ctx.effective("SRS") == "gradient"


class TestCategoricalColors:
    def test_categories_get_distinct_colors(self):
        pedigree = parse_pedigree(PREFIX_MODES)
        _, ctx = _render_context(pedigree)
        colors = ctx.categorical_colors["Diagnosis"]
        assert set(colors) == {"ASD", "Control", "ADHD"}
        assert len(set(colors.values())) == len(colors)

    def test_style_base_color_used_for_first_category(self):
        pedigree = parse_pedigree(PREFIX_MODES)
        style = PedigreeStyle()
        style.phenotype_colors = {"Diagnosis": "#123456"}
        _, ctx = _render_context(pedigree, style)
        assert ctx.categorical_colors["Diagnosis"]["ASD"] == "#123456"

    def test_categories_avoid_other_column_colors(self):
        pedigree = parse_pedigree(PREFIX_MODES)
        resolved, ctx = _render_context(pedigree)
        others = {c for n, c in resolved.phenotype_colors.items() if n != "Diagnosis"}
        assigned = set(ctx.categorical_colors["Diagnosis"].values())
        base = resolved.phenotype_colors["Diagnosis"]
        assert not (assigned - {base}) & others

    def test_resolution_is_reproducible(self):
        pedigree = parse_pedigree(PREFIX_MODES)
        _, first = _render_context(pedigree)
        _, second = _render_context(pedigree)
        assert first.categorical_colors == second.categorical_colors

    def test_style_change_is_honored_on_second_render(self):
        """Colors must not be cached across renders in the same process."""
        pedigree = parse_pedigree(PREFIX_MODES)
        seen = []
        for hex_color in ("#111111", "#00ff00"):
            style = PedigreeStyle()
            style.phenotype_colors = {"Diagnosis": hex_color}
            fig = render_pedigree(pedigree, style=style)
            plt.close(fig)
            _, ctx = _render_context(pedigree, style)
            seen.append(ctx.categorical_colors["Diagnosis"]["ASD"])
        assert seen == ["#111111", "#00ff00"]

    def test_no_category_leak_between_pedigrees(self):
        """A pedigree's legend must list only its own categories."""
        first = parse_pedigree(PREFIX_MODES)
        header = PREFIX_MODES.read_text().splitlines()[0]
        second = parse_pedigree_from_string(f"{header}\nFAM2\tp1\t0\t0\t1\t1\t100\t1\tSCZ\n")

        fig = render_pedigree(first)
        plt.close(fig)
        fig = render_pedigree(second)
        plt.close(fig)

        _, ctx = _render_context(second)
        assert set(ctx.categorical_colors["Diagnosis"]) == {"SCZ"}


class TestApiCliParity:
    def test_render_pedigree_bands_without_explicit_style(self):
        """render_pedigree must draw bands the CLI would draw."""
        pedigree = parse_pedigree(PREFIX_MODES)
        resolved, _ = _render_context(pedigree)
        assert resolved.phenotype_colors, "bands require auto-assigned colors"

        fig = render_pedigree(pedigree)
        # Legend axes are added on top of the family axes when bands are drawn.
        assert len(fig.axes) > 1
        plt.close(fig)

    def test_builder_render_is_repeatable(self):
        """Rendering twice must not compound PGS spacing on the shared style."""
        from pypedigree import PedigreeBuilder

        style = PedigreeStyle()
        builder = PedigreeBuilder("FAM1").set_style(style)
        builder.add_individual("a", sex="male", phenotypes={"PGS": 0.5})
        builder.add_individual("b", sex="female", phenotypes={"PGS": -0.5})

        plt.close(builder.render())
        plt.close(builder.render())
        assert style.couple_gap == PedigreeStyle().couple_gap
        assert style.phenotype_colors == {}


class TestStyleColorValidation:
    """Invalid colors are reported as PedigreeError, not matplotlib tracebacks."""

    @pytest.fixture
    def pedigree(self):
        return parse_pedigree(PREFIX_MODES)

    def test_invalid_phenotype_color_raises_pedigree_error(self, pedigree):
        style = PedigreeStyle()
        style.phenotype_colors = {"Epilepsy": "#zzzzzz"}
        with pytest.raises(PedigreeError, match=r"phenotype_colors\['Epilepsy'\]"):
            render_pedigree(pedigree, style=style)

    def test_invalid_outline_color_raises(self, pedigree):
        style = PedigreeStyle(outline_color="#12345")
        with pytest.raises(PedigreeError, match="outline_color"):
            render_pedigree(pedigree, style=style)

    def test_all_bad_colors_reported_at_once(self, pedigree):
        style = PedigreeStyle(outline_color="nope", fill_color="alsonope")
        style.phenotype_colors = {"Epilepsy": "#zzzzzz"}
        with pytest.raises(PedigreeError) as excinfo:
            render_pedigree(pedigree, style=style)
        message = str(excinfo.value)
        assert "outline_color" in message
        assert "fill_color" in message
        assert "phenotype_colors['Epilepsy']" in message

    def test_invalid_diverging_gradient_color_raises(self, pedigree):
        style = PedigreeStyle()
        style.gradient_diverging = {"IQ": {"min": 40, "max": 160, "color_low": "#nothex"}}
        with pytest.raises(PedigreeError, match="gradient_diverging"):
            render_pedigree(pedigree, style=style)

    def test_named_colors_are_accepted(self, pedigree):
        style = PedigreeStyle(outline_color="grey", line_color="darkblue")
        plt.close(render_pedigree(pedigree, style=style))

    def test_error_is_raised_before_any_drawing(self, pedigree):
        """Validation happens up front, so no half-drawn figure is left behind."""
        before = len(plt.get_fignums())
        style = PedigreeStyle(outline_color="nope")
        with pytest.raises(PedigreeError):
            render_pedigree(pedigree, style=style)
        assert len(plt.get_fignums()) == before


class TestVersion:
    def test_version_matches_installed_metadata(self):
        from importlib.metadata import version

        import pypedigree

        assert pypedigree.__version__ == version("pypedigree")

    def test_version_is_read_from_package_metadata(self):
        """pyproject.toml is the single source of the version."""
        from pathlib import Path

        source = (Path(__file__).parent.parent / "src/pypedigree/__init__.py").read_text()
        assert '_metadata_version("pypedigree")' in source

    def test_version_matches_pyproject(self):
        import re
        from pathlib import Path

        import pypedigree

        pyproject = (Path(__file__).parent.parent / "pyproject.toml").read_text()
        declared = re.search(r'^version = "([^"]+)"', pyproject, re.M).group(1)
        assert pypedigree.__version__ == declared


class TestDivergingGradient:
    """Bands and legend share one implementation, so they cannot drift apart."""

    CONF = {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}

    def test_midpoint_is_white(self):
        assert diverging_color_at(self.CONF, 0.5) == pytest.approx((1.0, 1.0, 1.0))

    def test_extremes_are_full_color(self):
        assert diverging_color_at(self.CONF, 0.0) == pytest.approx(to_rgb("#e74c3c"))
        assert diverging_color_at(self.CONF, 1.0) == pytest.approx(to_rgb("#2ecc71"))

    def test_value_and_position_agree(self):
        """The band (value-based) and legend (position-based) paths match."""
        for value, t in [(40, 0.0), (70, 0.25), (100, 0.5), (130, 0.75), (160, 1.0)]:
            assert diverging_color_for_value(self.CONF, value) == pytest.approx(
                diverging_color_at(self.CONF, t)
            )

    def test_values_outside_range_are_clamped(self):
        assert diverging_color_for_value(self.CONF, 0) == pytest.approx(
            diverging_color_for_value(self.CONF, 40)
        )
        assert diverging_color_for_value(self.CONF, 500) == pytest.approx(
            diverging_color_for_value(self.CONF, 160)
        )

    def test_degenerate_range_does_not_divide_by_zero(self):
        conf = {"min": 100, "mid": 100, "max": 100}
        assert len(diverging_color_for_value(conf, 100)) == 3

    def test_defaults_are_applied(self):
        """An entry with no colors falls back to the documented red/green."""
        assert diverging_color_at({}, 0.0) == pytest.approx(to_rgb("#e74c3c"))
        assert diverging_color_at({}, 1.0) == pytest.approx(to_rgb("#2ecc71"))


class TestNormalizeToUnit:
    def test_maps_range_onto_unit(self):
        assert normalize_to_unit(50, 0, 100) == pytest.approx(0.5)

    def test_clamps_outliers(self):
        assert normalize_to_unit(-10, 0, 100) == 0.0
        assert normalize_to_unit(200, 0, 100) == 1.0

    def test_degenerate_range_returns_midpoint(self):
        assert normalize_to_unit(5, 5, 5) == 0.5
