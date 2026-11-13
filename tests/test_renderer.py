"""Tests for pypedigree.renderer."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest
from conftest import TEST_DATA_DIR

from pypedigree.config import PedigreeStyle
from pypedigree.parser import parse_pedigree
from pypedigree.renderer import render_family, render_pedigree


@pytest.fixture
def style():
    return PedigreeStyle()


class TestRenderSmoke:
    """Smoke tests: verify rendering doesn't crash."""

    def test_simple_trio(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_two_gen_siblings(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/two_gen_siblings.ped")
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_three_gen(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_multi_family(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/multi_family.ped")
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_single_parent_with_phantom(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/single_parent.ped")
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_affected_mixed(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/affected_mixed.ped")
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)


class TestShapeTypes:
    def test_has_patches(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        # Should have at least 3 patches (one per individual)
        assert len(ax.patches) >= 3
        plt.close(fig)

    def test_affected_has_filled_shape(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        # At least one patch should be filled (black facecolor)
        filled = [p for p in ax.patches if p.get_facecolor()[:3] == (0.0, 0.0, 0.0)]
        assert len(filled) >= 1  # child is affected
        plt.close(fig)


class TestSaveOutput:
    def test_save_png(self, style, tmp_path):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        out = tmp_path / "test.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)

    def test_save_svg(self, style, tmp_path):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        out = tmp_path / "test.svg"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)

    def test_save_pdf(self, style, tmp_path):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        out = tmp_path / "test.pdf"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)


class TestRenderOnExistingAxes:
    def test_render_on_provided_axes(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fig, ax = plt.subplots(figsize=(10, 8))
        returned_fig = render_pedigree(ped, style=style, ax=ax)
        assert returned_fig is fig
        assert len(ax.patches) > 0
        plt.close(fig)


class TestRenderFamily:
    def test_render_single_family(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        fig, layout = render_family(fam, style=style)
        assert fig is not None
        assert len(layout.nodes) == 3
        plt.close(fig)


class TestPhantomRendering:
    def test_phantom_parent_dashed(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/single_parent.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        # At least one patch should have dashed linestyle
        dashed = [p for p in ax.patches if p.get_linestyle() == "--"]
        assert len(dashed) >= 1
        plt.close(fig)


class TestDeceasedRendering:
    def test_deceased_marker_drawn(self, style):
        """Deceased individuals should have extra line elements drawn (the diagonal)."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_deceased.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        # The deceased marker is drawn as ax.plot() which creates Line2D objects.
        # A trio without deceased has only connection lines.
        # With one deceased individual, there should be extra line(s).
        lines = ax.get_lines()
        assert len(lines) >= 4  # 3 connection lines + at least 1 deceased marker
        plt.close(fig)

    def test_deceased_renders_without_error(self, style, tmp_path):
        """Smoke test: extended file with deceased renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_deceased.ped")
        out = tmp_path / "deceased.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)


class TestProbandRendering:
    def test_proband_marker_drawn(self, style):
        """Proband individuals should have a 'P' text annotation."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_proband.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        # The proband marker draws a text "P" — check ax.texts
        text_contents = [t.get_text() for t in ax.texts]
        assert "P" in text_contents
        plt.close(fig)

    def test_proband_renders_without_error(self, style, tmp_path):
        """Smoke test: extended file with proband renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_proband.ped")
        out = tmp_path / "proband.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)


class TestExtendedFullRendering:
    def test_full_extended_renders(self, style, tmp_path):
        """Full 10-column file renders without error."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_full.ped")
        out = tmp_path / "extended_full.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)

    def test_deceased_and_proband_together(self, style):
        """File with both deceased and proband renders correctly."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_full.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        # Should have "P" for proband (child1)
        text_contents = [t.get_text() for t in ax.texts]
        assert "P" in text_contents
        # Should have extra lines for deceased (gfather)
        lines = ax.get_lines()
        assert len(lines) >= 1  # at least the deceased marker
        plt.close(fig)


class TestConsanguinityRendering:
    def test_consanguinous_renders_without_error(self, style, tmp_path):
        """Pedigree with consanguinous union renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        out = tmp_path / "consanguinity.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)

    def test_consanguinous_has_extra_lines(self, style):
        """Consanguinous relationship should produce extra line elements (double line)."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_sibling_marriage.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        lines = ax.get_lines()
        # A double relationship line produces 2 Line2D objects instead of 1
        # Normal trio has connection lines; with double there are extra lines
        assert len(lines) >= 4
        plt.close(fig)

    def test_uncle_niece_renders(self, style, tmp_path):
        """Uncle-niece consanguinity renders to PNG without error."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_uncle_niece.ped")
        out = tmp_path / "uncle_niece.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)


class TestBrokenRelationshipRendering:
    def test_broken_relationship_renders_without_error(self, style, tmp_path):
        """Programmatic broken union renders to PNG."""
        from pypedigree.models import Family, Individual, Pedigree, Sex, Union

        fam = Family(family_id="FAM1")
        fam.individuals["f"] = Individual(sample_id="f", family_id="FAM1", sex=Sex.MALE)
        fam.individuals["m"] = Individual(sample_id="m", family_id="FAM1", sex=Sex.FEMALE)
        fam.individuals["c"] = Individual(
            sample_id="c",
            family_id="FAM1",
            sex=Sex.MALE,
            father_id="f",
            mother_id="m",
        )
        fam.unions = [
            Union(
                partner_ids=("f", "m"),
                children_ids=["c"],
                is_broken=True,
            )
        ]
        ped = Pedigree(families={"FAM1": fam})
        out = tmp_path / "broken.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 0
        plt.close(fig)


class TestTwinRendering:
    def test_mz_twins_render_without_error(self, style, tmp_path):
        """MZ twin pedigree renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        out = tmp_path / "twins_mz.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_dz_twins_render_without_error(self, style, tmp_path):
        """DZ twin pedigree renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_dz.ped")
        out = tmp_path / "twins_dz.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_unknown_twins_render_with_question_mark(self, style):
        """Unknown zygosity twins render with '?' text annotation."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_unknown.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        text_contents = [t.get_text() for t in ax.texts]
        assert "?" in text_contents
        plt.close(fig)

    def test_triplets_render_without_error(self, style, tmp_path):
        """Triplet pedigree renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_triplets.ped")
        out = tmp_path / "triplets.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_twins_with_siblings_render(self, style, tmp_path):
        """Twins mixed with non-twin siblings render correctly."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_with_siblings.ped")
        out = tmp_path / "twins_siblings.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_mz_twins_have_extra_lines_for_bar(self, style):
        """MZ twins should produce extra Line2D elements for the horizontal bar."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        lines = ax.get_lines()
        # Descent line + 2 V-lines + 1 MZ bar = at least 4 lines
        # (plus relationship line = 5 total minimum)
        assert len(lines) >= 4
        plt.close(fig)


class TestDuplicateArcRendering:
    def test_duplicate_arc_renders_without_error(self, style, tmp_path):
        """Pedigree with duplicate individual renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        out = tmp_path / "multi_partner.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_male_multi_partner_renders(self, style, tmp_path):
        """Male with multiple partners renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_male.ped")
        out = tmp_path / "multi_male.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_complex_multi_partner_renders(self, style, tmp_path):
        """Complex multi-gen with duplicate renders to PNG."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_complex.ped")
        out = tmp_path / "multi_complex.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_duplicate_arc_creates_patch(self, style):
        """Duplicate arc should add a FancyArrowPatch to the axes."""
        import matplotlib.patches as mpatches

        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        fancy_patches = [p for p in ax.patches if isinstance(p, mpatches.FancyArrowPatch)]
        assert len(fancy_patches) >= 1
        # Check it uses dotted linestyle
        arc_patch = fancy_patches[0]
        assert arc_patch.get_linestyle() == ":"
        plt.close(fig)

    def test_duplicate_labels_show_original_name(self, style):
        """Both original and duplicate should show the same name in labels."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fig = render_pedigree(ped, style=style)
        ax = fig.axes[0]
        text_contents = [t.get_text() for t in ax.texts]
        # "mother" should appear twice (original + duplicate)
        assert text_contents.count("mother") == 2
        plt.close(fig)


class TestLegendNoOverlap:
    """Verify that legend text elements do not overlap."""

    def _get_legend_text_bboxes(self, fig):
        """Extract rendered bounding boxes for all text in the legend axes."""
        renderer = fig.canvas.get_renderer()
        # Legend axes is the last axes added (via fig.add_axes)
        legend_ax = None
        for ax in fig.axes:
            # Legend axes has xlim (0, 20) and is independent from subplot axes
            xl = ax.get_xlim()
            if abs(xl[1] - 20.0) < 0.1 and abs(xl[0]) < 0.1:
                legend_ax = ax
                break
        if legend_ax is None:
            return []
        bboxes = []
        for txt in legend_ax.texts:
            if not txt.get_text().strip():
                continue
            bb = txt.get_window_extent(renderer=renderer)
            bboxes.append((txt.get_text(), bb))
        return bboxes

    def _check_no_overlap(self, bboxes):
        """Check that no two text bounding boxes overlap significantly."""
        overlaps = []
        for i in range(len(bboxes)):
            for j in range(i + 1, len(bboxes)):
                name_i, bb_i = bboxes[i]
                name_j, bb_j = bboxes[j]
                # Check intersection
                x_overlap = max(0, min(bb_i.x1, bb_j.x1) - max(bb_i.x0, bb_j.x0))
                y_overlap = max(0, min(bb_i.y1, bb_j.y1) - max(bb_i.y0, bb_j.y0))
                overlap_area = x_overlap * y_overlap
                min_area = min(bb_i.width * bb_i.height, bb_j.width * bb_j.height)
                if min_area > 0 and overlap_area / min_area > 0.15:
                    overlaps.append((name_i, name_j, overlap_area / min_area))
        return overlaps

    def test_legend_no_text_overlap_multi_phenotype(self, style):
        """Legend text should not overlap with 3 band phenotypes."""
        ped = parse_pedigree(TEST_DATA_DIR / "phenotype/multi_phenotype_inline.ped")
        fig = render_pedigree(ped, style=style)
        fig.canvas.draw()
        bboxes = self._get_legend_text_bboxes(fig)
        overlaps = self._check_no_overlap(bboxes)
        assert not overlaps, f"Legend text overlaps detected: {overlaps}"
        plt.close(fig)

    def test_legend_no_text_overlap_prefix_modes(self, style):
        """Legend text should not overlap with mixed display modes."""
        ped = parse_pedigree(TEST_DATA_DIR / "phenotype/prefix_modes.ped")
        fig = render_pedigree(ped, style=style)
        fig.canvas.draw()
        bboxes = self._get_legend_text_bboxes(fig)
        overlaps = self._check_no_overlap(bboxes)
        assert not overlaps, f"Legend text overlaps detected: {overlaps}"
        plt.close(fig)


class TestEdgeCaseRendering:
    def test_empty_pedigree_renders(self, style):
        """An empty Pedigree (no families) should render without error."""
        from pypedigree.models import Pedigree

        ped = Pedigree(families={})
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_render_family_with_provided_axes(self, style):
        """render_family() with an explicit ax= argument should work."""
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        fig, ax = plt.subplots(figsize=(8, 6))
        returned_fig, layout = render_family(fam, style=style, ax=ax)
        assert returned_fig is fig
        assert len(layout.nodes) == 3
        assert len(ax.patches) > 0
        plt.close(fig)
