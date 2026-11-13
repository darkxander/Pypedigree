"""End-to-end integration tests for pypedigree."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest
from conftest import TEST_DATA_DIR

from pypedigree import (
    PedigreeStyle,
    Phenotype,
    Sex,
    compute_layout,
    parse_pedigree,
    parse_pedigree_from_string,
    render_pedigree,
)


@pytest.fixture
def style():
    return PedigreeStyle()


class TestEndToEnd:
    def test_string_to_figure(self, style):
        text = "FAM1\tf\t0\t0\t1\t1\nFAM1\tm\t0\t0\t2\t1\nFAM1\tc\tf\tm\t2\t2\n"
        ped = parse_pedigree_from_string(text)
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_file_to_saved_png(self, style, tmp_path):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
        out = tmp_path / "three_gen.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_all_valid_files_render(self, style):
        valid_files = [
            "simple_trio.ped",
            "two_gen_siblings.ped",
            "three_gen.ped",
            "single_parent.ped",
            "multi_family.ped",
            "comma_separated.csv",
            "space_separated.ped",
            "with_header.ped",
            "affected_mixed.ped",
            "consanguinity_first_cousins.ped",
            "consanguinity_uncle_niece.ped",
            "consanguinity_sibling_marriage.ped",
            "twins_mz.ped",
            "twins_dz.ped",
            "twins_unknown.ped",
            "twins_triplets.ped",
            "twins_with_siblings.ped",
            "twins_one_sided.ped",
            "dob_ordering.ped",
            "multiple_partners.ped",
            "multiple_partners_male.ped",
            "multiple_partners_complex.ped",
            "synesthesia.ped",
        ]
        for fname in valid_files:
            fpath = TEST_DATA_DIR / fname
            if fpath.exists():
                ped = parse_pedigree(fpath)
                fig = render_pedigree(ped, style=style)
                assert fig is not None, f"Failed to render {fname}"
                plt.close(fig)

    def test_consanguinity_end_to_end(self, style, tmp_path):
        """Full pipeline: parse consanguinous .ped -> detect -> layout -> render."""
        from pypedigree.models import EdgeType, LineStyle

        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        fam = ped.get_family("FAM1")

        # Verify detection
        cousin_union = fam.get_union_for_parents("cousin1", "cousin2")
        assert cousin_union.is_consanguinous is True

        # Verify layout
        layout = compute_layout(fam, style)
        double_edges = [
            e
            for e in layout.edges
            if e.edge_type == EdgeType.RELATIONSHIP and e.line_style == LineStyle.DOUBLE
        ]
        assert len(double_edges) == 1

        # Verify render
        out = tmp_path / "consanguinity_e2e.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_broken_relationship_api_end_to_end(self, style, tmp_path):
        """Full pipeline: programmatic broken union -> layout -> render."""
        from pypedigree.models import (
            EdgeType,
            Family,
            Individual,
            LineStyle,
            Pedigree,
            Union,
        )

        fam = Family(family_id="FAM1")
        fam.individuals["f"] = Individual(sample_id="f", family_id="FAM1", sex=Sex.MALE)
        fam.individuals["m"] = Individual(sample_id="m", family_id="FAM1", sex=Sex.FEMALE)
        fam.individuals["c1"] = Individual(
            sample_id="c1",
            family_id="FAM1",
            sex=Sex.MALE,
            father_id="f",
            mother_id="m",
        )
        fam.unions = [
            Union(
                partner_ids=("f", "m"),
                children_ids=["c1"],
                is_broken=True,
            )
        ]
        ped = Pedigree(families={"FAM1": fam})

        layout = compute_layout(fam, style)
        broken_edges = [
            e
            for e in layout.edges
            if e.edge_type == EdgeType.RELATIONSHIP and e.line_style == LineStyle.BROKEN
        ]
        assert len(broken_edges) == 1

        out = tmp_path / "broken_e2e.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        plt.close(fig)

    def test_programmatic_api(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")

        # Inspect data model
        assert fam.individuals["child"].phenotype == Phenotype.AFFECTED
        assert fam.individuals["father"].sex == Sex.MALE

        # Compute layout separately
        layout = compute_layout(fam, style)
        assert len(layout.nodes) == 3
        assert len([e for e in layout.edges if e.edge_type.value == "relationship"]) == 1

        # Render
        fig = render_pedigree(ped, style=style)
        assert fig is not None
        plt.close(fig)

    def test_twins_end_to_end(self, style, tmp_path):
        """Full pipeline: parse twins -> layout with connectors -> render."""
        from pypedigree.models import TwinType

        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")

        # Verify parsing
        assert len(fam.individuals["twin1"].twin_info) == 1

        # Verify layout
        layout = compute_layout(fam, style)
        assert len(layout.twin_connectors) == 1
        tc = layout.twin_connectors[0]
        assert tc.twin_type == TwinType.MONOZYGOTIC
        assert len(tc.child_positions) == 2

        # Verify render
        out = tmp_path / "twins_e2e.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_twins_with_siblings_end_to_end(self, style, tmp_path):
        """Full pipeline: twins + DOB-ordered siblings -> layout -> render."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_with_siblings.ped")
        fam = ped.get_family("FAM1")

        # Verify DOB + twin ordering in layout
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        assert node_map["older"].x < node_map["twin1"].x
        assert node_map["twin2"].x < node_map["younger"].x

        # Verify twin connector + sibship line coexistence
        assert len(layout.twin_connectors) == 1
        from pypedigree.models import EdgeType

        sibship_edges = [e for e in layout.edges if e.edge_type == EdgeType.SIBSHIP]
        assert len(sibship_edges) == 1

        # Verify render
        out = tmp_path / "twins_siblings_e2e.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_multiple_partners_end_to_end(self, style, tmp_path):
        """Full pipeline: parse multiple partners -> detect duplicates -> layout -> render."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")

        # Verify detection
        duplicates = [ind for ind in fam.individuals.values() if ind.is_duplicate]
        assert len(duplicates) == 1
        assert duplicates[0].duplicate_of == "mother"

        # Verify layout
        layout = compute_layout(fam, style)
        assert len(layout.duplicate_arcs) == 1
        dup_nodes = [n for n in layout.nodes if n.is_duplicate]
        assert len(dup_nodes) == 1
        assert dup_nodes[0].label == "mother"

        # Verify render
        out = tmp_path / "multi_partner_e2e.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_complex_multi_partner_end_to_end(self, style, tmp_path):
        """Full pipeline: multi-gen with multiple partners -> layout -> render."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_complex.ped")
        fam = ped.get_family("FAM1")

        # Verify detection
        duplicates = [ind for ind in fam.individuals.values() if ind.is_duplicate]
        assert len(duplicates) == 1
        assert duplicates[0].duplicate_of == "father"

        # Verify layout
        layout = compute_layout(fam, style)
        assert len(layout.duplicate_arcs) == 1
        # Original and duplicate should have the same label
        node_labels = {n.individual_id: n.label for n in layout.nodes}
        dup_id = [sid for sid in node_labels if sid.startswith("_dup_")][0]
        assert node_labels[dup_id] == "father"

        # Verify render
        out = tmp_path / "multi_complex_e2e.png"
        fig = render_pedigree(ped, style=style, output=str(out))
        assert out.exists()
        assert out.stat().st_size > 1000
        plt.close(fig)

    def test_synesthesia_pedigree_parsing(self, style):
        """Parse real-world synesthesia pedigree (45 individuals, couple 44+45)."""
        ped = parse_pedigree(TEST_DATA_DIR / "complex/synesthesia.ped")
        fam = ped.get_family("FAM_SYN01")
        assert fam is not None

        # 45 individuals: couple 44+45 with many children who marry out
        assert len(fam.individuals) == 45

        # Parents 44 and 45 are founders
        assert fam.individuals["44"].is_founder
        assert fam.individuals["45"].is_founder

        # Header detected: b_ASD parsed as binary phenotype
        # Primary phenotype comes from first phenotype column (b_ASD)
        assert fam.individuals["1"].phenotype == Phenotype.AFFECTED  # ASD=2
        assert fam.individuals["3"].phenotype == Phenotype.UNAFFECTED  # ASD=1

        # Multiple unions from children of 44+45 marrying different spouses
        assert len(fam.unions) >= 10

    def test_synesthesia_pedigree_layout(self, style):
        """Layout for synesthesia pedigree: no overlaps, 45 nodes positioned."""
        ped = parse_pedigree(TEST_DATA_DIR / "complex/synesthesia.ped")
        fam = ped.get_family("FAM_SYN01")
        layout = compute_layout(fam, style)

        assert len(layout.nodes) == 45

        # No overlapping nodes on the same generation row
        positions = [(n.x, n.y) for n in layout.nodes]
        for i, (x1, y1) in enumerate(positions):
            for x2, y2 in positions[i + 1 :]:
                if y1 == y2:
                    assert abs(x1 - x2) >= style.min_node_distance * 0.9

    def test_synesthesia_pedigree_render(self, style, tmp_path):
        """Full render of synesthesia pedigree to PNG and interactive SVG."""
        ped = parse_pedigree(TEST_DATA_DIR / "complex/synesthesia.ped")

        out_png = tmp_path / "synesthesia.png"
        fig = render_pedigree(ped, style=style, output=str(out_png))
        assert out_png.exists()
        assert out_png.stat().st_size > 5000
        plt.close(fig)

        from pypedigree.svg_renderer import render_interactive_svg

        out_svg = tmp_path / "synesthesia.svg"
        svg = render_interactive_svg(ped, style=style, output=str(out_svg))
        assert out_svg.exists()
        assert "<svg" in svg
        assert 'data-id="1"' in svg

    def test_multi_format_output(self, style, tmp_path):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        for ext in ["png", "svg", "pdf"]:
            out = tmp_path / f"test.{ext}"
            fig = render_pedigree(ped, style=style, output=str(out))
            assert out.exists(), f"Failed to save {ext}"
            assert out.stat().st_size > 0
            plt.close(fig)
