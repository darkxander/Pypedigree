"""Tests for the legend panel.

The legend is sized from its font and its row count, never from the figure,
and each row type is drawn by its own function. These tests cover the
geometry contract the row drawers rely on, plus the row composition of a
rendered figure.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest
from conftest import PHENOTYPE_DIR

from pypedigree import PedigreeStyle, parse_pedigree, parse_pedigree_from_string, render_pedigree
from pypedigree.phenotypes import (
    LEG_GRAD_W,
    LEG_X0,
    legend_metrics,
    pgs_legend_description,
)

PREFIX_MODES = PHENOTYPE_DIR / "prefix_modes.ped"
CARRIER = PHENOTYPE_DIR / "carrier.ped"
TEST_ALL = PHENOTYPE_DIR / "test_all.ped"
SYNESTHESIA = PHENOTYPE_DIR.parent / "complex" / "synesthesia.ped"


class TestLegendMetrics:
    def test_single_row_is_not_degenerate(self):
        m = legend_metrics(1, PedigreeStyle())
        assert m.n_rows == 1
        assert m.row_h > 0
        assert m.total_h > m.row_h

    def test_zero_rows_is_clamped_to_one(self):
        assert legend_metrics(0, PedigreeStyle()).n_rows == 1

    def test_rows_keep_their_height(self):
        """The panel grows with its rows; rows never shrink to fit."""
        heights = {legend_metrics(n, PedigreeStyle()).row_h for n in (1, 4, 50)}
        assert len(heights) == 1

    def test_total_height_grows_with_rows(self):
        totals = [legend_metrics(n, PedigreeStyle()).total_h for n in (1, 4, 10)]
        assert totals == sorted(totals)

    def test_binary_headers_get_their_own_room(self):
        plain = legend_metrics(3, PedigreeStyle())
        headed = legend_metrics(3, PedigreeStyle(), binary_headers=True)
        assert headed.total_h > plain.total_h

    def test_font_size_follows_style(self):
        small = legend_metrics(3, PedigreeStyle(font_size=6.0))
        large = legend_metrics(3, PedigreeStyle(font_size=14.0))
        assert large.fs > small.fs

    def test_units_scale_with_the_font(self):
        """Swatches and text grow together, so captions keep their room."""
        small = legend_metrics(3, PedigreeStyle(font_size=8.0))
        large = legend_metrics(3, PedigreeStyle(font_size=16.0))
        assert large.unit_x == pytest.approx(2 * small.unit_x)
        assert large.unit_y == pytest.approx(2 * small.unit_y)
        assert large.fs == pytest.approx(2 * small.fs)

    def test_swatches_do_not_overlap(self):
        m = legend_metrics(3, PedigreeStyle())
        for i in range(2):
            assert m.swatch_x(i) + m.sw_w <= m.swatch_x(i + 1)

    def test_row_label_clears_the_last_swatch(self):
        m = legend_metrics(3, PedigreeStyle())
        assert m.text_x > m.swatch_x(2) + m.sw_w

    def test_icon_label_clears_the_widest_icon(self):
        m = legend_metrics(3, PedigreeStyle())
        assert m.icon_text_x > LEG_X0 + m.sw_w

    def test_gradient_bar_fits_before_the_row_label(self):
        """The high-value caption is right-aligned on the bar's right edge."""
        m = legend_metrics(3, PedigreeStyle())
        assert m.text_x >= LEG_X0 + LEG_GRAD_W

    def test_rows_descend_and_do_not_overlap(self):
        m = legend_metrics(4, PedigreeStyle())
        centers = [m.row_center(i, m.total_h) for i in range(4)]
        assert centers == sorted(centers, reverse=True)
        for a, b in zip(centers[:-1], centers[1:], strict=True):
            assert a - b == pytest.approx(m.row_h)

    def test_metrics_are_immutable(self):
        m = legend_metrics(2, PedigreeStyle())
        with pytest.raises(AttributeError):
            m.row_h = 99.0


class TestPgsLegendDescription:
    def test_single_plain_pgs_column(self):
        text = pgs_legend_description(["PGS"])
        assert text == "PGS distribution (vertical) with individual score"

    def test_prefixed_columns_list_their_traits(self):
        text = pgs_legend_description(["PGS_autism", "PGS_IQ"])
        assert text == "PGS distribution (autism, IQ)"

    def test_unprefixed_name_is_kept_whole(self):
        assert pgs_legend_description(["PGSautism"]) == "PGS distribution (PGSautism)"

    def test_single_prefixed_column_still_lists_its_trait(self):
        assert pgs_legend_description(["PGS_ADHD"]) == "PGS distribution (ADHD)"


class TestLegendPhysicalSize:
    """The panel is sized by its font, not by the figure it sits in.

    It used to take a fixed share of the figure, with swatches sized on the
    figure width and text in points: a wide pedigree got huge swatches and
    tiny text, a circular one a panel stretched over the whole height, and
    a larger font made the captions overlap.
    """

    @staticmethod
    def _inches(fig, ax):
        return ax.get_position().height * fig.get_size_inches()[1]

    def test_same_legend_height_whatever_the_figure(self):
        """Same rows, same font: same panel, for a 51 cm-wide hierarchical
        figure and a square circular one."""
        wide = render_pedigree(parse_pedigree(SYNESTHESIA))
        circular = render_pedigree(parse_pedigree(SYNESTHESIA), circular=True)
        try:
            wide_h = self._inches(wide, wide.axes[1])
            assert wide_h == pytest.approx(self._inches(circular, circular.axes[1]))
            assert wide_h < 2.0
        finally:
            plt.close(wide)
            plt.close(circular)

    def test_legend_adds_room_instead_of_squeezing_the_pedigree(self):
        with_legend = render_pedigree(parse_pedigree(SYNESTHESIA))
        without = render_pedigree(parse_pedigree(SYNESTHESIA), show_legend=False)
        try:
            main_with = self._inches(with_legend, with_legend.axes[0])
            main_without = self._inches(without, without.axes[0])
            assert main_with == pytest.approx(main_without)
            legend_h = self._inches(with_legend, with_legend.axes[1])
            assert with_legend.get_size_inches()[1] == pytest.approx(
                without.get_size_inches()[1] + legend_h
            )
        finally:
            plt.close(with_legend)
            plt.close(without)

    @pytest.mark.parametrize("font_size", [8.0, 12.0, 20.0])
    def test_column_headers_never_overlap(self, font_size):
        fig = render_pedigree(parse_pedigree(TEST_ALL), style=PedigreeStyle(font_size=font_size))
        try:
            legend = fig.axes[1]
            headers = [
                t for t in legend.texts if t.get_text() in ("Affected", "Unaffected", "Unknown")
            ]
            boxes = sorted((t.get_window_extent() for t in headers), key=lambda b: b.x0)
            assert len(boxes) == 3
            for left, right in zip(boxes, boxes[1:], strict=False):
                assert left.x1 <= right.x0
        finally:
            plt.close(fig)

    def test_border_encloses_every_entry(self):
        fig = render_pedigree(parse_pedigree(TEST_ALL), style=PedigreeStyle(font_size=14.0))
        try:
            legend = fig.axes[1]
            border = max(legend.patches, key=lambda p: p.get_window_extent().width)
            frame = border.get_window_extent()
            for text in legend.texts:
                box = text.get_window_extent()
                assert frame.x0 <= box.x0 and box.x1 <= frame.x1, text.get_text()
        finally:
            plt.close(fig)


class TestLegendComposition:
    """The legend axes exists exactly when there is something to explain."""

    def test_legend_axes_added_for_phenotypes(self):
        fig = render_pedigree(parse_pedigree(PREFIX_MODES))
        assert len(fig.axes) == 2
        plt.close(fig)

    def test_no_legend_axes_when_disabled(self):
        fig = render_pedigree(parse_pedigree(PREFIX_MODES), show_legend=False)
        assert len(fig.axes) == 1
        plt.close(fig)

    def test_no_legend_axes_without_phenotypes_or_carriers(self):
        text = "F\ta\t0\t0\t1\t1\nF\tb\t0\t0\t2\t2\nF\tc\ta\tb\t1\t2\n"
        fig = render_pedigree(parse_pedigree_from_string(text))
        assert len(fig.axes) == 1
        plt.close(fig)

    def test_carriers_alone_get_a_legend(self):
        fig = render_pedigree(parse_pedigree(CARRIER))
        assert len(fig.axes) == 2
        plt.close(fig)

    def test_legend_grows_with_more_phenotype_columns(self):
        """More rows reserve more of the figure for the legend."""
        one = render_pedigree(parse_pedigree(PREFIX_MODES))
        many = render_pedigree(parse_pedigree(TEST_ALL))
        one_h = one.axes[1].get_position().height
        many_h = many.axes[1].get_position().height
        assert many_h >= one_h
        plt.close(one)
        plt.close(many)

    def test_legend_renders_for_every_row_type(self, tmp_path):
        """Binary, gradient, diverging, categorical, outline, PGS and carrier."""
        text = (
            "FID\tID\tFather\tMother\tSex\tb_Bin\tg_Grad\tgi_Inv\t"
            "o_Out\tc_Cat\tIQ\tPGS_x\tcarrier\n"
            "F\ta\t0\t0\t1\t2\t90\t70\t1\tASD\t105\t0.4\t1\n"
            "F\tb\t0\t0\t2\t1\t70\t95\t0\tControl\t118\t-0.2\t2\n"
            "F\tc\ta\tb\t1\t2\t80\t60\t1\tADHD\t92\t0.9\t0\n"
        )
        pedigree = parse_pedigree_from_string(text)
        style = PedigreeStyle()
        style.gradient_diverging = {"IQ": {"min": 40, "mid": 100, "max": 160}}
        out = tmp_path / "all_row_types.png"
        fig = render_pedigree(pedigree, style=style, output=out)
        assert len(fig.axes) == 2
        plt.close(fig)
        assert out.stat().st_size > 0
