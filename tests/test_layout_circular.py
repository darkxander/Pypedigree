"""Tests for the circular (radial) layout.

`compute_layout_circular` places founders at the centre and each later
generation on its own concentric ring, in the left-to-right order of the
hierarchical layout. These tests pin the geometric invariants that the
renderer relies on.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import pytest
from conftest import (
    BASIC_DIR,
    COMPLEX_DIR,
    CONSANG_DIR,
    PARTNERS_DIR,
    PHENOTYPE_DIR,
    TWINS_DIR,
)

from pypedigree import (
    EdgeType,
    LineStyle,
    PedigreeStyle,
    compute_layout,
    compute_layout_circular,
    parse_pedigree,
    parse_pedigree_from_string,
    render_pedigree,
)
from pypedigree.layout import _arc_waypoints, _assign_generations
from pypedigree.models import Family, Pedigree, Sex
from pypedigree.renderer import _render_layout, _shape_boundary_radius, _shape_lower_extent

SIMPLE_TRIO = BASIC_DIR / "simple_trio.ped"
THREE_GEN = BASIC_DIR / "three_gen.ped"
TWO_GEN_SIBLINGS = BASIC_DIR / "two_gen_siblings.ped"
CONSANGUINITY = CONSANG_DIR / "consanguinity.ped"
SYNESTHESIA = COMPLEX_DIR / "synesthesia.ped"
TWINS_MZ = TWINS_DIR / "twins_mz.ped"
MULTI_PARTNER = PARTNERS_DIR / "multiple_partners.ped"

# Every parseable fixture, so an invariant break shows up wherever it happens.
ALL_FIXTURES = [
    BASIC_DIR / "simple_trio.ped",
    BASIC_DIR / "single_parent.ped",
    BASIC_DIR / "three_gen.ped",
    BASIC_DIR / "two_gen_siblings.ped",
    BASIC_DIR / "affected_mixed.ped",
    CONSANG_DIR / "consanguinity.ped",
    CONSANG_DIR / "consanguinity_first_cousins.ped",
    CONSANG_DIR / "consanguinity_uncle_niece.ped",
    PARTNERS_DIR / "multiple_partners.ped",
    PARTNERS_DIR / "multiple_partners_complex.ped",
    PARTNERS_DIR / "multiple_relationships.ped",
    PHENOTYPE_DIR / "carrier.ped",
    PHENOTYPE_DIR / "carrier_genes.ped",
    PHENOTYPE_DIR / "test_all.ped",
    TWINS_DIR / "twins_mz.ped",
    TWINS_DIR / "twins_triplets.ped",
    TWINS_DIR / "twins_with_siblings.ped",
    COMPLEX_DIR / "synesthesia.ped",
]


def first_family(path) -> Family:
    """Load a fixture and return its first family."""
    return next(iter(parse_pedigree(path).families.values()))


def radius(node) -> float:
    """Distance of a node from the centre of the diagram."""
    return math.hypot(node.x, node.y)


def radii_by_generation(family: Family, layout) -> dict[int, set[float]]:
    """Map each generation to the set of radii its nodes were placed at."""
    generations = _assign_generations(family)
    grouped: dict[int, set[float]] = {}
    for node in layout.nodes:
        gen = generations.get(node.individual_id)
        grouped.setdefault(gen, set()).add(round(radius(node), 6))
    return grouped


class TestCircularLayoutBasics:
    def test_empty_family_gives_empty_layout(self):
        layout = compute_layout_circular(Family(family_id="EMPTY"))
        assert layout.family_id == "EMPTY"
        assert layout.nodes == []
        assert layout.edges == []

    def test_family_id_is_preserved(self):
        layout = compute_layout_circular(first_family(SIMPLE_TRIO))
        assert layout.family_id == "FAM1"

    def test_default_style_is_used_when_none(self):
        family = first_family(SIMPLE_TRIO)
        assert compute_layout_circular(family) == compute_layout_circular(family, PedigreeStyle())

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_every_individual_gets_one_node(self, path):
        family = first_family(path)
        layout = compute_layout_circular(family)
        assert len(layout.nodes) == len(family.individuals)
        assert {n.individual_id for n in layout.nodes} == set(family.individuals)

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_all_coordinates_are_finite(self, path):
        layout = compute_layout_circular(first_family(path))
        for node in layout.nodes:
            assert math.isfinite(node.x)
            assert math.isfinite(node.y)
            assert math.isfinite(node.rotation)

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_layout_is_deterministic(self, path):
        family = first_family(path)
        first = compute_layout_circular(family)
        second = compute_layout_circular(family)
        assert [(n.individual_id, n.x, n.y, n.rotation) for n in first.nodes] == [
            (n.individual_id, n.x, n.y, n.rotation) for n in second.nodes
        ]


class TestCircularRings:
    """Generations map onto concentric rings."""

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_one_radius_per_generation(self, path):
        family = first_family(path)
        grouped = radii_by_generation(family, compute_layout_circular(family))
        for gen, values in grouped.items():
            assert len(values) == 1, f"generation {gen} spread over radii {sorted(values)}"

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_radius_increases_with_generation(self, path):
        family = first_family(path)
        grouped = radii_by_generation(family, compute_layout_circular(family))
        radii = [next(iter(grouped[g])) for g in sorted(grouped)]
        assert radii == sorted(radii)
        assert len(set(radii)) == len(radii), "two generations share a ring"

    def test_founders_sit_near_the_centre(self):
        family = first_family(THREE_GEN)
        layout = compute_layout_circular(family)
        generations = _assign_generations(family)
        style = PedigreeStyle()
        for node in layout.nodes:
            if generations.get(node.individual_id) == 0:
                # Offset from the exact centre only by the couple gap.
                assert radius(node) < style.base_radius / 2

    def test_first_descendant_ring_uses_base_radius(self):
        style = PedigreeStyle(base_radius=4.0)
        family = first_family(SIMPLE_TRIO)
        grouped = radii_by_generation(family, compute_layout_circular(family, style))
        assert next(iter(grouped[1])) == pytest.approx(4.0)

    def test_later_rings_are_spaced_by_radial_gap(self):
        style = PedigreeStyle(base_radius=4.0, radial_gap=3.0)
        family = first_family(THREE_GEN)
        grouped = radii_by_generation(family, compute_layout_circular(family, style))
        assert next(iter(grouped[2])) - next(iter(grouped[1])) == pytest.approx(3.0)

    def test_dense_ring_expands_base_radius(self):
        """A ring too crowded for base_radius is pushed outward, not overlapped."""
        family = first_family(SYNESTHESIA)
        style = PedigreeStyle(base_radius=3.0)
        grouped = radii_by_generation(family, compute_layout_circular(family, style))
        assert next(iter(grouped[1])) > style.base_radius

    def test_sparse_ring_keeps_base_radius(self):
        family = first_family(SIMPLE_TRIO)
        style = PedigreeStyle(base_radius=3.0)
        grouped = radii_by_generation(family, compute_layout_circular(family, style))
        assert next(iter(grouped[1])) == pytest.approx(3.0)


class TestCircularSpacing:
    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_no_two_nodes_overlap(self, path):
        style = PedigreeStyle()
        layout = compute_layout_circular(first_family(path), style)
        points = [(n.x, n.y) for n in layout.nodes]
        for i, a in enumerate(points):
            for b in points[i + 1 :]:
                assert math.dist(a, b) >= style.shape_size, (
                    f"nodes closer than one shape width in {path.stem}"
                )

    def test_larger_shapes_still_do_not_overlap(self):
        style = PedigreeStyle(shape_size=1.6)
        layout = compute_layout_circular(first_family(SYNESTHESIA), style)
        points = [(n.x, n.y) for n in layout.nodes]
        closest = min(math.dist(a, b) for i, a in enumerate(points) for b in points[i + 1 :])
        assert closest >= style.shape_size


class TestCircularNodes:
    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_shapes_stay_upright(self, path):
        """A square turned to face outward would read as the unknown-sex diamond."""
        layout = compute_layout_circular(first_family(path))
        assert all(n.rotation == 0.0 for n in layout.nodes)

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_labels_point_radially_outward(self, path):
        """Away from the lines, which reach each node from the centre."""
        style = PedigreeStyle()
        layout = compute_layout_circular(first_family(path), style)
        for node in layout.nodes:
            if radius(node) <= style.shape_size:
                assert node.label_angle is None
            else:
                direction = math.atan2(node.y, node.x)
                assert math.cos(node.label_angle - direction) == pytest.approx(1.0)

    def test_pgs_curve_sides_differ_within_a_couple(self):
        """Partners put their PGS curve on opposite sides, not on top of each other."""
        family = first_family(PHENOTYPE_DIR / "test_all.ped")
        layout = compute_layout_circular(family)
        sides = {n.individual_id: n.pgs_side for n in layout.nodes}
        for union in family.unions:
            p1, p2 = union.partner_ids
            if p1 in sides and p2 in sides:
                assert sides[p1] != sides[p2]

    def test_node_metadata_survives_the_layout(self):
        family = first_family(BASIC_DIR / "affected_mixed.ped")
        linear = {n.individual_id: n for n in compute_layout(family).nodes}
        circular = {n.individual_id: n for n in compute_layout_circular(family).nodes}
        for sid, node in circular.items():
            assert node.sex == linear[sid].sex
            assert node.phenotype == linear[sid].phenotype
            assert node.label == linear[sid].label

    def test_extended_markers_survive_the_layout(self):
        family = first_family(BASIC_DIR.parent / "extended/extended_markers.ped")
        linear = {n.individual_id: n for n in compute_layout(family).nodes}
        circular = {n.individual_id: n for n in compute_layout_circular(family).nodes}
        for sid, node in circular.items():
            assert node.is_deceased == linear[sid].is_deceased
            assert node.is_proband == linear[sid].is_proband

    def test_phenotype_values_survive_the_layout(self):
        family = first_family(PHENOTYPE_DIR / "test_all.ped")
        circular = compute_layout_circular(family)
        assert any(n.extra_phenotypes for n in circular.nodes)
        for node in circular.nodes:
            source = family.individuals[node.individual_id]
            assert node.extra_phenotypes == source.extra_phenotypes


class TestCircularEdges:
    def test_every_union_gets_a_relationship_edge(self):
        family = first_family(TWO_GEN_SIBLINGS)
        layout = compute_layout_circular(family)
        relationships = [e for e in layout.edges if e.edge_type is EdgeType.RELATIONSHIP]
        assert len(relationships) == len(family.unions)

    def test_single_child_union_gets_a_descent_edge(self):
        """Typed the same way the hierarchical layout types it."""
        family = first_family(THREE_GEN)
        circular = compute_layout_circular(family)
        linear = compute_layout(family)
        assert sum(1 for e in circular.edges if e.edge_type is EdgeType.DESCENT) == sum(
            1 for e in linear.edges if e.edge_type is EdgeType.DESCENT
        )

    def test_multi_child_union_gets_descent_sibship_and_stubs(self):
        family = first_family(TWO_GEN_SIBLINGS)
        layout = compute_layout_circular(family)
        union = family.unions[0]
        assert len(union.children_ids) > 1
        assert sum(1 for e in layout.edges if e.edge_type is EdgeType.DESCENT) == 1
        assert sum(1 for e in layout.edges if e.edge_type is EdgeType.SIBSHIP) == 1
        assert sum(1 for e in layout.edges if e.edge_type is EdgeType.INDIVIDUAL) == len(
            union.children_ids
        )

    def test_sibship_edge_is_an_arc_of_constant_radius(self):
        layout = compute_layout_circular(first_family(TWO_GEN_SIBLINGS))
        sibships = [e for e in layout.edges if e.edge_type is EdgeType.SIBSHIP]
        assert sibships, "expected a sibship arc for a multi-child union"
        arc = sibships[0]
        assert arc.waypoints, "sibship arc must be drawn as a polyline"
        points = [arc.start, *arc.waypoints, arc.end]
        radii = [math.hypot(x, y) for x, y in points]
        assert max(radii) - min(radii) < 1e-6

    def test_consanguinous_union_uses_a_double_line(self):
        layout = compute_layout_circular(first_family(CONSANGUINITY))
        styles = {e.line_style for e in layout.edges if e.edge_type is EdgeType.RELATIONSHIP}
        assert LineStyle.DOUBLE in styles

    def test_broken_union_uses_a_broken_line(self):
        text = (
            "FID\tID\tFather\tMother\tSex\tPheno\n"
            "F\tdad\t0\t0\t1\t1\n"
            "F\tmum\t0\t0\t2\t1\n"
            "F\tkid\tdad\tmum\t1\t2\n"
        )
        family = next(iter(parse_pedigree_from_string(text).families.values()))
        family.unions[0].is_broken = True
        layout = compute_layout_circular(family)
        relationship = next(e for e in layout.edges if e.edge_type is EdgeType.RELATIONSHIP)
        assert relationship.line_style is LineStyle.BROKEN

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_all_edge_endpoints_are_finite(self, path):
        layout = compute_layout_circular(first_family(path))
        for edge in layout.edges:
            for x, y in [edge.start, edge.end, *edge.waypoints]:
                assert math.isfinite(x)
                assert math.isfinite(y)

    def test_relationship_edges_join_actual_node_positions(self):
        layout = compute_layout_circular(first_family(TWO_GEN_SIBLINGS))
        positions = {(round(n.x, 9), round(n.y, 9)) for n in layout.nodes}
        for edge in layout.edges:
            if edge.edge_type is EdgeType.RELATIONSHIP:
                assert (round(edge.start[0], 9), round(edge.start[1], 9)) in positions
                assert (round(edge.end[0], 9), round(edge.end[1], 9)) in positions


class TestArcWaypoints:
    def test_points_lie_on_the_requested_radius(self):
        points = _arc_waypoints(0.0, math.pi / 2, 5.0, n_points=8)
        for x, y in points:
            assert math.hypot(x, y) == pytest.approx(5.0)

    def test_endpoints_match_the_requested_angles(self):
        points = _arc_waypoints(0.0, math.pi / 2, 5.0, n_points=8)
        assert points[0] == pytest.approx((5.0, 0.0))
        assert points[-1] == pytest.approx((0.0, 5.0), abs=1e-9)

    def test_point_count_is_honoured(self):
        assert len(_arc_waypoints(0.0, 1.0, 2.0, n_points=12)) == 12

    def test_zero_span_arc_is_degenerate_but_valid(self):
        points = _arc_waypoints(1.0, 1.0, 2.0, n_points=4)
        for x, y in points:
            assert math.hypot(x, y) == pytest.approx(2.0)


class TestCircularRendering:
    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_renders_without_error(self, path, tmp_path):
        pedigree = parse_pedigree(path)
        out = tmp_path / f"{path.stem}.png"
        fig = render_pedigree(pedigree, circular=True, output=out)
        plt.close(fig)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_figure_is_square(self):
        fig = render_pedigree(parse_pedigree(THREE_GEN), circular=True)
        width, height = fig.get_size_inches()
        assert width == pytest.approx(height)
        plt.close(fig)

    def test_axes_use_equal_aspect(self):
        fig = render_pedigree(parse_pedigree(THREE_GEN), circular=True)
        assert fig.axes[0].get_aspect() == 1.0
        plt.close(fig)

    def test_shapes_are_enlarged_relative_to_linear(self):
        """Circular mode scales shapes up; the caller's style is untouched."""
        style = PedigreeStyle(shape_size=0.8)
        plt.close(render_pedigree(parse_pedigree(THREE_GEN), style=style, circular=True))
        assert style.shape_size == 0.8

    def test_renders_with_phenotype_bands(self, tmp_path):
        out = tmp_path / "pheno_circular.png"
        fig = render_pedigree(
            parse_pedigree(PHENOTYPE_DIR / "test_all.ped"), circular=True, output=out
        )
        plt.close(fig)
        assert out.stat().st_size > 0

    def test_renders_multiple_families(self, tmp_path):
        out = tmp_path / "multi_circular.png"
        fig = render_pedigree(
            parse_pedigree(BASIC_DIR / "multi_family.ped"), circular=True, output=out
        )
        plt.close(fig)
        assert out.stat().st_size > 0

    def test_empty_pedigree_renders(self):
        fig = render_pedigree(Pedigree(families={}), circular=True)
        plt.close(fig)

    def test_large_pedigree_completes(self):
        fig = render_pedigree(parse_pedigree(SYNESTHESIA), circular=True)
        plt.close(fig)


class TestRotatedNodeLabels:
    """Labels and proband markers of rotated nodes stay readable.

    The circular layout keeps shapes upright, but the renderer still honours
    ``NodeLayout.rotation`` for callers who set it. A square at 45° reaches
    √2 × half below its centre, so a label placed at a fixed offset was
    hidden behind its own shape, and the proband marker rotated with the
    shape. These tests pin the fix on the circular layout with every node
    turned to face outward, as 0.4.0 development versions drew it.
    """

    @pytest.fixture
    def rendered(self):
        """Render the large fixture's circular layout, shapes turned outward."""
        from dataclasses import replace

        family = first_family(SYNESTHESIA)
        style = PedigreeStyle()
        base = compute_layout_circular(family, style)
        layout = replace(
            base,
            nodes=[
                replace(
                    n,
                    rotation=(
                        math.atan2(n.y, n.x) - math.pi / 2 if radius(n) > style.shape_size else 0.0
                    ),
                    label_angle=None,
                )
                for n in base.nodes
            ],
        )
        fig, ax = plt.subplots()
        _render_layout(ax, layout, style, show_labels=True)
        yield layout, style, ax
        plt.close(fig)

    @staticmethod
    def _label_of(ax, node):
        """The Text artist carrying a node's ID label."""
        matches = [t for t in ax.texts if t.get_text() == node.label and t.get_va() == "top"]
        return min(matches, key=lambda t: abs(t.get_position()[0] - node.x))

    def test_fixture_has_rotated_squares(self, rendered):
        layout, _, _ = rendered
        assert any(n.rotation != 0.0 and n.sex is Sex.MALE for n in layout.nodes)

    def test_every_label_sits_below_its_shape(self, rendered):
        layout, style, ax = rendered
        half = style.shape_size / 2
        for node in layout.nodes:
            if not node.label:
                continue
            _, label_top = self._label_of(ax, node).get_position()
            lowest = node.y - _shape_lower_extent(node.sex, half, node.rotation)
            assert label_top < lowest, f"label of {node.individual_id} starts inside its shape"

    def test_labels_are_not_pushed_further_than_needed(self, rendered):
        """Nodes already clear of their shape keep the standard offset."""
        layout, style, ax = rendered
        half = style.shape_size / 2
        standard = half + style.label_offset
        for node in layout.nodes:
            if not node.label:
                continue
            _, label_top = self._label_of(ax, node).get_position()
            needed = _shape_lower_extent(node.sex, half, node.rotation) + style.label_offset
            assert node.y - label_top <= max(standard, needed) + 1e-9

    def test_rotated_labels_have_a_white_backing(self, rendered):
        layout, _, ax = rendered
        for node in layout.nodes:
            if node.label and node.rotation != 0.0:
                assert self._label_of(ax, node).get_bbox_patch() is not None

    def test_proband_letter_stays_upright(self, rendered):
        _, _, ax = rendered
        letters = [t for t in ax.texts if t.get_text() == "P"]
        assert letters, "fixture must contain a proband"
        assert all(t.get_rotation() == 0.0 for t in letters)

    def test_proband_arrow_points_at_the_rotated_outline(self, rendered):
        layout, style, ax = rendered
        half = style.shape_size / 2
        proband = next(n for n in layout.nodes if n.is_proband)
        assert proband.rotation != 0.0, "the check is only meaningful on a rotated node"
        arrows = [t for t in ax.texts if hasattr(t, "xy") and t.get_text() == ""]
        tip = min(arrows, key=lambda a: math.dist(a.xy, (proband.x, proband.y))).xy
        direction = math.atan2(tip[1] - proband.y, tip[0] - proband.x)
        assert math.degrees(direction) % 360 == pytest.approx(225)
        expected = _shape_boundary_radius(proband.sex, half, proband.rotation, direction)
        assert math.dist(tip, (proband.x, proband.y)) == pytest.approx(expected)

    def test_no_text_is_upside_down(self):
        """Band values and NA marks turn with the shape, but never past 90°."""
        from dataclasses import replace

        from pypedigree.config import resolve_style
        from pypedigree.phenotypes import resolve_phenotype_render

        pedigree = parse_pedigree(PHENOTYPE_DIR / "test_all.ped")
        names, modes = pedigree.phenotype_names, pedigree.phenotype_display_modes
        style = resolve_style(PedigreeStyle(), names)
        base = compute_layout_circular(next(iter(pedigree.families.values())), style)
        layout = replace(
            base,
            nodes=[
                replace(n, rotation=math.atan2(n.y, n.x) - math.pi / 2, label_angle=None)
                for n in base.nodes
            ],
        )
        fig, ax = plt.subplots()
        try:
            render = resolve_phenotype_render(pedigree, names, modes, style)
            _render_layout(ax, layout, style, False, names, modes, render)
            turned = [text for text in ax.texts if text.get_rotation() != 0.0]
            assert turned, "the check needs band text on rotated nodes"
            for text in turned:
                angle = (text.get_rotation() + 180.0) % 360.0 - 180.0  # in [-180, 180)
                assert -90.0 <= angle <= 90.0, f"{text.get_text()!r} is upside down"
        finally:
            plt.close(fig)

    def test_unrotated_proband_marker_is_unchanged(self):
        """The hierarchical layout keeps the arrow on the bounding-box corner."""
        family = first_family(BASIC_DIR.parent / "extended/extended_proband.ped")
        style = PedigreeStyle()
        layout = compute_layout(family, style)
        fig, ax = plt.subplots()
        _render_layout(ax, layout, style, show_labels=True)
        half = style.shape_size / 2
        proband = next(n for n in layout.nodes if n.is_proband)
        arrow = next(t for t in ax.texts if hasattr(t, "xy") and t.get_text() == "")
        assert arrow.xy == pytest.approx((proband.x - half, proband.y - half))
        plt.close(fig)


class TestShapeGeometry:
    """The outline helpers behind label and proband placement."""

    H = 0.5

    def test_upright_square_reaches_half_below(self):
        assert _shape_lower_extent(Sex.MALE, self.H, 0.0) == pytest.approx(self.H)

    def test_square_at_45_degrees_reaches_root_two_half(self):
        extent = _shape_lower_extent(Sex.MALE, self.H, math.radians(45))
        assert extent == pytest.approx(self.H * math.sqrt(2))

    @pytest.mark.parametrize("deg", [0, 17, 45, 90, 133, 270])
    def test_circle_extent_ignores_rotation(self, deg):
        assert _shape_lower_extent(Sex.FEMALE, self.H, math.radians(deg)) == self.H

    def test_diamond_extent_shrinks_when_rotated_onto_its_side(self):
        upright = _shape_lower_extent(Sex.UNKNOWN, self.H, 0.0)
        tilted = _shape_lower_extent(Sex.UNKNOWN, self.H, math.radians(45))
        assert tilted < upright

    def test_square_boundary_reaches_the_corner_on_the_diagonal(self):
        r = _shape_boundary_radius(Sex.MALE, self.H, 0.0, math.radians(225))
        assert r == pytest.approx(self.H * math.sqrt(2))

    def test_square_boundary_reaches_the_side_on_an_axis(self):
        r = _shape_boundary_radius(Sex.MALE, self.H, 0.0, math.radians(180))
        assert r == pytest.approx(self.H)

    def test_boundary_follows_rotation(self):
        """Rotating the square by 45° puts a side, not a corner, on the diagonal."""
        r = _shape_boundary_radius(Sex.MALE, self.H, math.radians(45), math.radians(225))
        assert r == pytest.approx(self.H)

    def test_circle_boundary_is_its_radius(self):
        r = _shape_boundary_radius(Sex.FEMALE, self.H, 1.0, 2.0)
        assert r == pytest.approx(self.H)


class TestRadialLabels:
    """Labels of the circular layout sit outside their shape, radially."""

    def test_labels_clear_their_shape_outward(self):
        family = first_family(SYNESTHESIA)
        style = PedigreeStyle()
        layout = compute_layout_circular(family, style)
        fig, ax = plt.subplots()
        try:
            _render_layout(ax, layout, style, show_labels=True)
            half = style.shape_size / 2
            texts = {t.get_text(): t for t in ax.texts}
            for node in layout.nodes:
                if node.label_angle is None or not node.label:
                    continue
                text = texts[node.label]
                reach = _shape_boundary_radius(node.sex, half, 0.0, node.label_angle)
                assert math.dist(text.get_position(), (node.x, node.y)) >= reach
                assert text.get_bbox_patch() is not None, "white backing against lines"
        finally:
            plt.close(fig)


class TestCircularTwinsAndDuplicates:
    """Twin connectors and duplicate arcs, as in the hierarchical layout."""

    def test_twins_hang_from_one_point_of_the_sibship_arc(self):
        family = first_family(TWINS_MZ)
        layout = compute_layout_circular(family)
        assert len(layout.twin_connectors) == len(compute_layout(family).twin_connectors) == 1
        connector = layout.twin_connectors[0]
        assert len(connector.child_positions) == 2

    def test_twins_get_no_individual_line_of_their_own(self):
        family = first_family(TWINS_DIR / "twins_with_siblings.ped")
        layout = compute_layout_circular(family)
        twins = {m for n in layout.twin_connectors for m in n.child_positions}
        assert twins
        node_at = {(round(n.x, 6), round(n.y, 6)): n.individual_id for n in layout.nodes}
        stubs = {
            node_at.get((round(e.end[0], 6), round(e.end[1], 6)))
            for e in layout.edges
            if e.edge_type is EdgeType.INDIVIDUAL
        }
        assert stubs == {"older", "younger"}

    def test_only_twins_need_no_sibship_arc(self):
        """A union whose children are all triplets has one connection point."""
        layout = compute_layout_circular(first_family(TWINS_DIR / "twins_triplets.ped"))
        assert not [e for e in layout.edges if e.edge_type is EdgeType.SIBSHIP]
        assert len(layout.twin_connectors) == 1

    def test_duplicates_are_joined_to_their_original(self):
        family = first_family(PARTNERS_DIR / "multiple_partners_complex.ped")
        layout = compute_layout_circular(family)
        assert len(layout.duplicate_arcs) == 1
        dup = next(n for n in layout.nodes if n.is_duplicate)
        original = next(n for n in layout.nodes if n.individual_id == "father")
        assert radius(dup) == pytest.approx(radius(original)), "same generation, same ring"

    @pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.stem)
    def test_sibship_arcs_of_one_ring_never_overlap(self, path):
        """Inherited from the hierarchical order: sibships stay apart."""
        layout = compute_layout_circular(first_family(path))
        arcs: dict[float, list[tuple[float, float]]] = {}
        for e in layout.edges:
            if e.edge_type is not EdgeType.SIBSHIP:
                continue
            points = [e.start, *e.waypoints, e.end]
            angles = [math.atan2(y, x) for x, y in points]
            unwrapped = [angles[0]]
            for a in angles[1:]:
                while a - unwrapped[-1] > math.pi:
                    a -= 2 * math.pi
                while a - unwrapped[-1] < -math.pi:
                    a += 2 * math.pi
                unwrapped.append(a)
            arcs.setdefault(round(math.hypot(*e.start), 6), []).append(
                (min(unwrapped), max(unwrapped))
            )
        for spans in arcs.values():
            spans.sort()
            for (_, hi), (lo, _) in zip(spans, spans[1:], strict=False):
                assert lo >= hi - 1e-9
