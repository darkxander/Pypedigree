"""Tests for pypedigree.layout."""

import math

import pytest
from conftest import TEST_DATA_DIR

from pypedigree.config import PedigreeStyle
from pypedigree.layout import _assign_generations, compute_layout
from pypedigree.models import EdgeType, Family
from pypedigree.parser import parse_pedigree


@pytest.fixture
def style():
    return PedigreeStyle()


@pytest.fixture
def trio_family():
    ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
    return ped.get_family("FAM1")


@pytest.fixture
def siblings_family():
    ped = parse_pedigree(TEST_DATA_DIR / "basic/two_gen_siblings.ped")
    return ped.get_family("FAM1")


@pytest.fixture
def three_gen_family():
    ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
    return ped.get_family("FAM1")


class TestGenerationAssignment:
    def test_trio_generations(self, trio_family):
        gens = _assign_generations(trio_family)
        assert gens["father"] == 0
        assert gens["mother"] == 0
        assert gens["child"] == 1

    def test_partners_same_generation(self, trio_family):
        gens = _assign_generations(trio_family)
        assert gens["father"] == gens["mother"]

    def test_three_gen(self, three_gen_family):
        gens = _assign_generations(three_gen_family)
        assert gens["gfather"] == 0
        assert gens["gmother"] == 0
        assert gens["father"] == 1
        assert gens["mother"] == 1
        assert gens["child"] == 2

    def test_siblings_same_generation(self, siblings_family):
        gens = _assign_generations(siblings_family)
        assert gens["son1"] == gens["daughter1"] == gens["son2"]
        assert gens["son1"] == 1

    def test_children_below_parents(self, three_gen_family):
        gens = _assign_generations(three_gen_family)
        for ind in three_gen_family.individuals.values():
            if ind.father_id and ind.father_id in gens:
                assert gens[ind.sample_id] > gens[ind.father_id]
            if ind.mother_id and ind.mother_id in gens:
                assert gens[ind.sample_id] > gens[ind.mother_id]

    def test_empty_family(self, style):
        fam = Family(family_id="EMPTY")
        layout = compute_layout(fam, style)
        assert layout.nodes == []
        assert layout.edges == []


class TestHorizontalPositioning:
    def test_partners_side_by_side(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        father = node_map["father"]
        mother = node_map["mother"]
        # Partners should be separated by couple_gap
        assert abs(abs(father.x - mother.x) - style.couple_gap) < 0.01

    def test_child_centered_under_parents(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        father = node_map["father"]
        mother = node_map["mother"]
        child = node_map["child"]
        parent_center = (father.x + mother.x) / 2
        assert abs(child.x - parent_center) < 0.01

    def test_siblings_spread_horizontally(self, siblings_family, style):
        layout = compute_layout(siblings_family, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        children_x = [node_map[sid].x for sid in ["son1", "daughter1", "son2"]]
        # Children should be spread out
        for i in range(1, len(children_x)):
            assert children_x[i] > children_x[i - 1]

    def test_no_overlaps(self, siblings_family, style):
        layout = compute_layout(siblings_family, style)
        # Group nodes by y (generation)
        gen_groups: dict[float, list[float]] = {}
        for node in layout.nodes:
            gen_groups.setdefault(node.y, []).append(node.x)
        for _y, xs in gen_groups.items():
            xs_sorted = sorted(xs)
            for i in range(1, len(xs_sorted)):
                assert xs_sorted[i] - xs_sorted[i - 1] >= style.min_node_distance - 0.01

    def test_three_gen_vertical_order(self, three_gen_family, style):
        layout = compute_layout(three_gen_family, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        # Grandparents above parents above child (y decreases with generation)
        assert node_map["gfather"].y > node_map["father"].y
        assert node_map["father"].y > node_map["child"].y


class TestConnectionLines:
    def test_trio_has_relationship_line(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        rel_edges = [e for e in layout.edges if e.edge_type == EdgeType.RELATIONSHIP]
        assert len(rel_edges) == 1

    def test_trio_has_descent_line(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        desc_edges = [e for e in layout.edges if e.edge_type == EdgeType.DESCENT]
        assert len(desc_edges) == 1

    def test_trio_has_individual_line(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        ind_edges = [e for e in layout.edges if e.edge_type == EdgeType.INDIVIDUAL]
        # Single child connects directly from descent
        assert len(ind_edges) == 1

    def test_siblings_have_sibship_line(self, siblings_family, style):
        layout = compute_layout(siblings_family, style)
        sib_edges = [e for e in layout.edges if e.edge_type == EdgeType.SIBSHIP]
        assert len(sib_edges) == 1

    def test_siblings_have_individual_lines(self, siblings_family, style):
        layout = compute_layout(siblings_family, style)
        ind_edges = [e for e in layout.edges if e.edge_type == EdgeType.INDIVIDUAL]
        # 3 children = 3 individual lines
        assert len(ind_edges) == 3

    def test_three_gen_has_two_relationship_lines(self, three_gen_family, style):
        layout = compute_layout(three_gen_family, style)
        rel_edges = [e for e in layout.edges if e.edge_type == EdgeType.RELATIONSHIP]
        assert len(rel_edges) == 2

    def test_relationship_line_is_horizontal(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        rel_edge = [e for e in layout.edges if e.edge_type == EdgeType.RELATIONSHIP][0]
        assert rel_edge.start[1] == rel_edge.end[1]  # Same y

    def test_descent_line_is_vertical(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        desc_edge = [e for e in layout.edges if e.edge_type == EdgeType.DESCENT][0]
        assert desc_edge.start[0] == desc_edge.end[0]  # Same x
        assert desc_edge.start[1] > desc_edge.end[1]  # Goes downward


class TestLabels:
    def test_real_individuals_have_labels(self, trio_family, style):
        layout = compute_layout(trio_family, style)
        for node in layout.nodes:
            if not node.individual_id.startswith("_phantom_"):
                assert node.label != ""

    def test_phantom_has_empty_label(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "basic/single_parent.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        phantom_nodes = [n for n in layout.nodes if n.is_phantom]
        assert len(phantom_nodes) == 1
        assert phantom_nodes[0].label == ""


class TestNodeAttributes:
    def test_deceased_propagated_to_node(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_deceased.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        assert node_map["father"].is_deceased is True
        assert node_map["mother"].is_deceased is False
        assert node_map["child"].is_deceased is False

    def test_proband_propagated_to_node(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_proband.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        assert node_map["child"].is_proband is True
        assert node_map["father"].is_proband is False
        assert node_map["mother"].is_proband is False

    def test_affected_propagated_to_node(self, style):
        ped = parse_pedigree(TEST_DATA_DIR / "extended/extended_full.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        from pypedigree.models import Phenotype

        assert node_map["child1"].phenotype == Phenotype.AFFECTED
        assert node_map["child2"].phenotype == Phenotype.UNAFFECTED


class TestConsanguinityLineStyle:
    def test_consanguinous_union_produces_double_line(self, style):
        """Consanguinous union should produce DOUBLE line style on relationship edge."""
        from pypedigree.models import LineStyle

        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_sibling_marriage.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        rel_edges = [e for e in layout.edges if e.edge_type == EdgeType.RELATIONSHIP]
        double_edges = [e for e in rel_edges if e.line_style == LineStyle.DOUBLE]
        assert len(double_edges) == 1

    def test_non_consanguinous_union_has_solid_line(self, style):
        """Normal unions should have SOLID line style."""
        from pypedigree.models import LineStyle

        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        rel_edges = [e for e in layout.edges if e.edge_type == EdgeType.RELATIONSHIP]
        for e in rel_edges:
            assert e.line_style == LineStyle.SOLID

    def test_first_cousin_has_one_double_line(self, style):
        """In first-cousin pedigree, only the cousin union should be DOUBLE."""
        from pypedigree.models import LineStyle

        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        rel_edges = [e for e in layout.edges if e.edge_type == EdgeType.RELATIONSHIP]
        double_edges = [e for e in rel_edges if e.line_style == LineStyle.DOUBLE]
        solid_edges = [e for e in rel_edges if e.line_style == LineStyle.SOLID]
        assert len(double_edges) == 1  # cousin1 + cousin2
        assert len(solid_edges) >= 2  # gfather+gmother, father+mother, uncle+aunt


class TestBrokenRelationshipLineStyle:
    def test_broken_union_produces_broken_line(self, style):
        """Union with is_broken=True should produce BROKEN line style."""
        from pypedigree.models import Individual, LineStyle, Sex, Union

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
        layout = compute_layout(fam, style)
        rel_edges = [e for e in layout.edges if e.edge_type == EdgeType.RELATIONSHIP]
        assert len(rel_edges) == 1
        assert rel_edges[0].line_style == LineStyle.BROKEN


class TestCoupleGrouping:
    def test_first_cousin_couples_grouped(self, style):
        """In first-cousin pedigree, couples should be adjacent blocks."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}

        # father and mother should be adjacent (couple block)
        assert abs(abs(node_map["father"].x - node_map["mother"].x) - style.couple_gap) < 0.01

        # uncle and aunt should be adjacent (couple block)
        assert abs(abs(node_map["uncle"].x - node_map["aunt"].x) - style.couple_gap) < 0.01

        # The two couple blocks should NOT interleave:
        # all of {father, mother} should be left of or right of all of {uncle, aunt}
        fm_max = max(node_map["father"].x, node_map["mother"].x)
        ua_min = min(node_map["uncle"].x, node_map["aunt"].x)
        assert fm_max < ua_min or ua_min > fm_max  # non-overlapping

    def test_first_cousins_drawn_as_a_couple(self, style):
        """Consanguinous cousins stand side by side, a couple gap apart."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}

        gap = abs(node_map["cousin1"].x - node_map["cousin2"].x)
        assert gap == pytest.approx(style.couple_gap)

    def test_sibling_marriage_partners_adjacent(self, style):
        """In sibling marriage, son and daughter should remain properly spaced."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_sibling_marriage.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}

        # Son and daughter should be at couple_gap distance
        assert abs(abs(node_map["son"].x - node_map["daughter"].x) - style.couple_gap) < 0.01


class TestDobOrdering:
    def test_children_ordered_by_dob(self, style):
        """Children with DOB should be ordered oldest (left) to youngest (right)."""
        ped = parse_pedigree(TEST_DATA_DIR / "extended/dob_ordering.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}

        # oldest (1990) < middle (1995) < youngest (1999)
        assert node_map["oldest"].x < node_map["middle"].x
        assert node_map["middle"].x < node_map["youngest"].x

    def test_twins_with_siblings_ordered_by_dob(self, style):
        """In twins_with_siblings.ped: older < twin1,twin2 < younger."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_with_siblings.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}

        assert node_map["older"].x < node_map["twin1"].x
        assert node_map["twin1"].x < node_map["twin2"].x
        assert node_map["twin2"].x < node_map["younger"].x


class TestTwinAdjacency:
    def test_mz_twins_adjacent(self, style):
        """MZ twins should be placed adjacent to each other."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}

        gap = abs(node_map["twin1"].x - node_map["twin2"].x)
        # Twins should be within sibling_gap distance
        assert gap <= style.sibling_gap + 0.01

    def test_triplets_all_adjacent(self, style):
        """Triplets should all be placed adjacent."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_triplets.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}

        xs = sorted([node_map["trip1"].x, node_map["trip2"].x, node_map["trip3"].x])
        # Each pair should be within sibling_gap
        for i in range(1, len(xs)):
            assert xs[i] - xs[i - 1] <= style.sibling_gap + 0.01


class TestTwinConnectors:
    def test_mz_twins_have_connector(self, style):
        """MZ twins produce a TwinConnector in the layout."""
        from pypedigree.models import TwinType

        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.twin_connectors) == 1
        tc = layout.twin_connectors[0]
        assert tc.twin_type == TwinType.MONOZYGOTIC
        assert len(tc.child_positions) == 2

    def test_dz_twins_have_connector(self, style):
        """DZ twins produce a TwinConnector."""
        from pypedigree.models import TwinType

        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_dz.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.twin_connectors) == 1
        assert layout.twin_connectors[0].twin_type == TwinType.DIZYGOTIC

    def test_unknown_twins_have_connector(self, style):
        """Unknown zygosity twins produce a TwinConnector."""
        from pypedigree.models import TwinType

        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_unknown.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.twin_connectors) == 1
        assert layout.twin_connectors[0].twin_type == TwinType.UNKNOWN

    def test_triplets_have_one_connector_with_three_children(self, style):
        """Triplets produce one TwinConnector with 3 child positions."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_triplets.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.twin_connectors) == 1
        tc = layout.twin_connectors[0]
        assert len(tc.child_positions) == 3

    def test_no_twin_connectors_for_normal_pedigree(self, style):
        """Non-twin pedigrees have no twin connectors."""
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert layout.twin_connectors == []

    def test_twins_only_no_sibship_line(self, style):
        """When all children are twins, no sibship line should be produced."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        sibship_edges = [e for e in layout.edges if e.edge_type == EdgeType.SIBSHIP]
        assert len(sibship_edges) == 0

    def test_twins_no_individual_lines(self, style):
        """Twin members should not have individual lines (only V-connectors)."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        individual_edges = [e for e in layout.edges if e.edge_type == EdgeType.INDIVIDUAL]
        assert len(individual_edges) == 0

    def test_twins_with_siblings_has_sibship_line(self, style):
        """Twins mixed with non-twins should have a sibship line."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_with_siblings.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        sibship_edges = [e for e in layout.edges if e.edge_type == EdgeType.SIBSHIP]
        assert len(sibship_edges) == 1

    def test_twins_with_siblings_non_twin_individual_lines(self, style):
        """Non-twin siblings should still have individual lines."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_with_siblings.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        individual_edges = [e for e in layout.edges if e.edge_type == EdgeType.INDIVIDUAL]
        # older and younger have individual lines; twin1/twin2 don't
        assert len(individual_edges) == 2

    def test_convergence_point_centered(self, style):
        """Twin convergence point should be centered between twins."""
        ped = parse_pedigree(TEST_DATA_DIR / "twins/twins_mz.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        tc = layout.twin_connectors[0]
        twin_center = (node_map["twin1"].x + node_map["twin2"].x) / 2
        assert abs(tc.convergence_point[0] - twin_center) < 0.01


def _nodes_inside_twin_groups(family, layout) -> list[tuple[str, list[str]]]:
    """(node, twin group) for every non-twin drawn strictly between twins of a group."""
    from pypedigree.layout import _find_twin_groups

    found = []
    for group, _tt in _find_twin_groups(family):
        twins = [n for n in layout.nodes if n.individual_id in group]
        if len(twins) < 2:
            continue
        lo, hi = min(n.x for n in twins), max(n.x for n in twins)
        row = round(twins[0].y, 6)
        for node in layout.nodes:
            if (
                node.individual_id not in group
                and round(node.y, 6) == row
                and lo + 1e-6 < node.x < hi - 1e-6
            ):
                found.append((node.individual_id, sorted(group)))
    return found


class TestTwinSpouses:
    """A married twin's partner is drawn outside the twin group.

    The partner used to be placed on the side given by the partnership's
    partner order, which could put them between two twins; the twin
    connector then spanned across the spouse.
    """

    # Columns: FID ID father mother sex phenotype DOB twins deceased proband
    CASES = {
        "left twin married": """
            F f 0 0 1 1 0 0 0 0
            F m 0 0 2 1 0 0 0 0
            F t1 f m 1 1 0 DZ:t2 0 0
            F t2 f m 2 1 0 DZ:t1 0 0
            F s 0 0 2 1 0 0 0 0
            F k t1 s 1 1 0 0 0 0
        """,
        "right twin married, wife listed second": """
            F f 0 0 1 1 0 0 0 0
            F m 0 0 2 1 0 0 0 0
            F t1 f m 1 1 0 MZ:t2 0 0
            F t2 f m 2 1 0 MZ:t1 0 0
            F h 0 0 1 1 0 0 0 0
            F k h t2 2 1 0 0 0 0
        """,
        "both twins married": """
            F f 0 0 1 1 0 0 0 0
            F m 0 0 2 1 0 0 0 0
            F t1 f m 2 1 0 DZ:t2 0 0
            F t2 f m 1 1 0 DZ:t1 0 0
            F h 0 0 1 1 0 0 0 0
            F w 0 0 2 1 0 0 0 0
            F k1 h t1 1 1 0 0 0 0
            F k2 t2 w 2 1 0 0 0 0
        """,
        "end triplets married": """
            F f 0 0 1 1 0 0 0 0
            F m 0 0 2 1 0 0 0 0
            F a f m 1 1 0 TW:b,TW:c 0 0
            F b f m 2 1 0 TW:a,TW:c 0 0
            F c f m 2 1 0 TW:a,TW:b 0 0
            F w 0 0 2 1 0 0 0 0
            F h 0 0 1 1 0 0 0 0
            F k1 a w 1 1 0 0 0 0
            F k2 h c 2 1 0 0 0 0
        """,
        "married twins among siblings": """
            F f 0 0 1 1 0 0 0 0
            F m 0 0 2 1 0 0 0 0
            F o f m 2 1 1990 0 0 0
            F t1 f m 1 1 1992 DZ:t2 0 0
            F t2 f m 2 1 1992 DZ:t1 0 0
            F y f m 1 1 1995 0 0 0
            F s1 0 0 2 1 0 0 0 0
            F s2 0 0 1 1 0 0 0 0
            F s3 0 0 2 1 0 0 0 0
            F k1 t1 s1 1 1 0 0 0 0
            F k2 s2 t2 2 1 0 0 0 0
            F k3 y s3 1 1 0 0 0 0
        """,
    }

    @pytest.mark.parametrize("case", list(CASES))
    def test_no_spouse_between_twins(self, style, case):
        family = _family_from_rows(self.CASES[case])
        layout = compute_layout(family, style)
        assert _nodes_inside_twin_groups(family, layout) == []
        assert _overlapping_sibships(layout) == []

    def test_spouse_on_outer_side(self, style):
        family = _family_from_rows(self.CASES["left twin married"])
        x = {n.individual_id: n.x for n in compute_layout(family, style).nodes}
        assert x["s"] < x["t1"] < x["t2"]
        assert min(x["s"], x["t1"]) < x["k"] < max(x["s"], x["t1"])

    def test_no_spouse_between_twins_in_fixtures(self, style):
        for family in TestOrthogonalDescent._all_fixture_families():
            inside = _nodes_inside_twin_groups(family, compute_layout(family, style))
            assert not inside, f"{family.family_id}: drawn between twins {inside}"


class TestDuplicateNodes:
    def test_duplicate_has_separate_position(self, style):
        """Duplicate node should have its own x position, not same as original."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        original = node_map["mother"]
        dup = [n for n in layout.nodes if n.is_duplicate][0]
        assert dup.x != original.x

    def test_duplicate_label_shows_original_name(self, style):
        """Duplicate node should show the original individual's name as label."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        dup_nodes = [n for n in layout.nodes if n.is_duplicate]
        assert len(dup_nodes) == 1
        assert dup_nodes[0].label == "mother"

    def test_duplicate_preserves_sex(self, style):
        """Duplicate node should have same sex as original."""
        from pypedigree.models import Sex

        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        dup_nodes = [n for n in layout.nodes if n.is_duplicate]
        assert dup_nodes[0].sex == Sex.FEMALE

    def test_no_duplicate_arcs_for_simple_pedigree(self, style):
        """Pedigree without duplicates should have no duplicate arcs."""
        ped = parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert layout.duplicate_arcs == []


class TestDuplicateArcs:
    def test_one_arc_for_one_duplicate(self, style):
        """One duplicate should produce one DuplicateArc."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.duplicate_arcs) == 1

    def test_arc_connects_original_to_duplicate(self, style):
        """Arc from_pos should match original, to_pos should match duplicate."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        node_map = {n.individual_id: n for n in layout.nodes}
        arc = layout.duplicate_arcs[0]
        original = node_map["mother"]
        dup = [n for n in layout.nodes if n.is_duplicate][0]
        assert arc.from_pos == (original.x, original.y)
        assert arc.to_pos == (dup.x, dup.y)

    def test_complex_has_arc(self, style):
        """Complex multi-gen pedigree with duplicate should have an arc."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_complex.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.duplicate_arcs) == 1

    def test_male_multi_partner_has_arc(self, style):
        """Male with multiple partners should produce a duplicate arc."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_male.ped")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.duplicate_arcs) == 1


class TestEdgeCaseLayouts:
    def test_single_individual_layout(self, style):
        """A single founder with no relationships should produce 1 node, 0 edges."""
        from pypedigree.parser import parse_pedigree_from_string

        ped = parse_pedigree_from_string("FAM1\tfounder\t0\t0\t1\t1\n")
        fam = ped.get_family("FAM1")
        layout = compute_layout(fam, style)
        assert len(layout.nodes) == 1
        assert len(layout.edges) == 0
        assert layout.nodes[0].individual_id == "founder"

    def test_empty_family_layout(self, style):
        """An empty family should produce an empty layout."""
        fam = Family(family_id="EMPTY")
        layout = compute_layout(fam, style)
        assert len(layout.nodes) == 0
        assert len(layout.edges) == 0


class TestOrthogonalDescent:
    """Pedigree nomenclature uses horizontal and vertical lines only.

    A single child that could not be centred under its parents used to be
    joined by one diagonal line. It is now reached the way a sibship is: a
    short horizontal segment at sibship level, then a vertical line down.
    """

    @staticmethod
    def _all_fixture_families():
        import glob
        import warnings

        families = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for path in sorted(glob.glob(str(TEST_DATA_DIR / "*" / "*.ped"))):
                try:
                    families += list(parse_pedigree(path).families.values())
                except Exception:
                    continue
        from conftest import (
            generate_deep_pedigree,
            generate_realistic_pedigree,
            generate_wide_pedigree,
        )

        families += [generate_deep_pedigree(), generate_wide_pedigree()]
        families += [generate_realistic_pedigree(n, seed) for n in (100, 200) for seed in (1, 42)]
        return families

    def test_every_hierarchical_edge_is_horizontal_or_vertical(self, style):
        for family in self._all_fixture_families():
            for edge in compute_layout(family, style).edges:
                (x0, y0), (x1, y1) = edge.start, edge.end
                assert x0 == pytest.approx(x1) or y0 == pytest.approx(y1), (
                    f"{family.family_id}: diagonal {edge.edge_type.value} edge "
                    f"from {edge.start} to {edge.end}"
                )

    def test_every_descent_line_reaches_its_children(self, style):
        """The bottom of each descent line touches a sibship line, a child's
        vertical line, or a twin convergence point — never empty space."""
        for family in self._all_fixture_families():
            layout = compute_layout(family, style)
            sibships = [e for e in layout.edges if e.edge_type is EdgeType.SIBSHIP]
            stubs = [e for e in layout.edges if e.edge_type is EdgeType.INDIVIDUAL]
            for descent in (e for e in layout.edges if e.edge_type is EdgeType.DESCENT):
                x, y = descent.end
                on_sibship = any(
                    s.start[1] == pytest.approx(y)
                    and min(s.start[0], s.end[0]) - 1e-9 <= x <= max(s.start[0], s.end[0]) + 1e-9
                    for s in sibships
                )
                on_stub = any(e.start == pytest.approx((x, y)) for e in stubs)
                on_twins = any(
                    t.convergence_point == pytest.approx((x, y)) for t in layout.twin_connectors
                )
                assert on_sibship or on_stub or on_twins, (
                    f"{family.family_id}: descent line ends at {descent.end}, touching nothing"
                )

    def test_offset_single_child_gets_a_horizontal_bridge(self, style, trio_family):
        """An only child placed off the parents' midpoint.

        The layout centres an only child under their parents unless the
        child is drawn next to a partner from another branch (a loop), so
        positions are set by hand.
        """
        from pypedigree.layout import _compute_edges

        father, mother, child = (
            trio_family.unions[0].partner_ids[0],
            trio_family.unions[0].partner_ids[1],
            trio_family.unions[0].children_ids[0],
        )
        x = {father: 0.0, mother: style.couple_gap, child: -1.0}
        generations = {father: 0, mother: 0, child: 1}
        edges = _compute_edges(trio_family, x, generations, style)
        stub = next(e for e in edges if e.edge_type is EdgeType.INDIVIDUAL)
        bridge = [e for e in edges if e.edge_type is EdgeType.SIBSHIP]
        assert stub.start[0] == pytest.approx(stub.end[0]), "the child's own line is vertical"
        assert len(bridge) == 1, "the offset child must hang from a horizontal segment"
        assert sorted((bridge[0].start[0], bridge[0].end[0])) == pytest.approx(
            [-1.0, style.couple_gap / 2]
        )

    def test_centred_single_child_needs_no_bridge(self, style, trio_family):
        layout = compute_layout(trio_family, style)
        assert not [e for e in layout.edges if e.edge_type is EdgeType.SIBSHIP]

    def test_offset_twin_group_is_connected_to_the_descent_line(self, style):
        """A twin pair as only children, placed off the parents' midpoint.

        No fixture produces this, so positions are set by hand: before the
        fix the descent line stopped beside the twins' convergence point.
        """
        from pypedigree.layout import _compute_edges, _compute_twin_connectors, _find_twin_groups
        from pypedigree.parser import parse_pedigree_from_string

        text = (
            "FID\tID\tFather\tMother\tSex\tPheno\tDOB\tTwins\n"
            "F\tdad\t0\t0\t1\t1\t0\t0\n"
            "F\tmum\t0\t0\t2\t1\t0\t0\n"
            "F\tt1\tdad\tmum\t1\t1\t0\tMZ:t2\n"
            "F\tt2\tdad\tmum\t1\t1\t0\tMZ:t1\n"
        )
        family = next(iter(parse_pedigree_from_string(text).families.values()))
        generations = _assign_generations(family)
        x = {"dad": 0.0, "mum": 1.5, "t1": 3.0, "t2": 4.5}  # twins centred at 3.75, parents at 0.75
        groups = _find_twin_groups(family)
        assert groups, "fixture must form a twin group"
        edges = _compute_edges(family, x, generations, style, groups)
        connector = _compute_twin_connectors(family, x, generations, style, groups)[0]
        cx, cy = connector.convergence_point
        descent = next(e for e in edges if e.edge_type is EdgeType.DESCENT)
        bridge = [
            e for e in edges if e.edge_type is EdgeType.SIBSHIP and e.start[1] == pytest.approx(cy)
        ]
        assert bridge, "a horizontal segment must join the descent line to the twins"
        lo, hi = sorted((bridge[0].start[0], bridge[0].end[0]))
        assert lo == pytest.approx(min(descent.end[0], cx))
        assert hi == pytest.approx(max(descent.end[0], cx))


def _sibship_spans(layout):
    """Sibship segments grouped by height: {y: [(x_left, x_right), ...]}."""
    levels: dict[float, list[tuple[float, float]]] = {}
    for edge in layout.edges:
        if edge.edge_type is EdgeType.SIBSHIP:
            xs = sorted((edge.start[0], edge.end[0]))
            levels.setdefault(round(edge.start[1], 6), []).append((xs[0], xs[1]))
    return levels


def _overlapping_sibships(layout):
    overlaps = []
    for y, spans in _sibship_spans(layout).items():
        spans.sort()
        for i, (_lo, hi) in enumerate(spans):
            for lo2, hi2 in spans[i + 1 :]:
                if lo2 >= hi - 1e-6:
                    break
                overlaps.append((y, (_lo, hi), (lo2, hi2)))
    return overlaps


class TestSibshipSeparation:
    """Sibship lines of different unions never share a stretch of line.

    When they did, the children of two couples hung from one continuous
    line and could not be told apart. This happened in one small tree
    pedigree in ten and in nearly every pedigree of 50 individuals or more.
    """

    @staticmethod
    def _families():
        from conftest import (
            generate_clinical_pedigree,
            generate_deep_pedigree,
            generate_realistic_pedigree,
            generate_wide_pedigree,
        )

        families = TestOrthogonalDescent._all_fixture_families()
        families += [generate_clinical_pedigree(n, s) for n in (8, 12, 15) for s in range(10)]
        families += [generate_clinical_pedigree(n, s) for n in (50, 100, 200) for s in range(3)]
        families += [
            generate_clinical_pedigree(n, s, cousin_marriages=k)
            for n in (30, 120)
            for s in range(4)
            for k in (1, 3)
        ]
        families += [generate_realistic_pedigree(150, s) for s in (7, 99)]
        families += [generate_deep_pedigree(), generate_wide_pedigree()]
        return families

    def test_two_cousins_sibships_stay_apart(self, style):
        """An 11-person tree that drew two sibships as one line."""
        from pypedigree.parser import parse_pedigree_from_string

        text = (
            "F\tf0\t0\t0\t1\t1\n"
            "F\tf1\t0\t0\t2\t1\n"
            "F\ti2\tf0\tf1\t1\t1\n"
            "F\ti3\tf0\tf1\t1\t1\n"
            "F\ts4\t0\t0\t2\t1\n"
            "F\ts5\t0\t0\t2\t1\n"
            "F\ti6\ti3\ts5\t2\t1\n"
            "F\ti7\ti3\ts5\t1\t1\n"
            "F\ts8\t0\t0\t1\t1\n"
            "F\ts9\t0\t0\t2\t1\n"
            "F\ti10\ti2\ts4\t1\t1\n"
        )
        family = parse_pedigree_from_string(text).get_family("F")
        assert _overlapping_sibships(compute_layout(family, style)) == []

    def test_no_two_sibship_lines_overlap(self, style):
        for family in self._families():
            overlaps = _overlapping_sibships(compute_layout(family, style))
            assert not overlaps, f"{family.family_id}: overlapping sibship lines {overlaps[:3]}"

    def test_nodes_on_a_row_keep_their_distance(self, style):
        gap = min(style.sibling_gap, style.couple_gap)
        for family in self._families():
            rows: dict[float, list[float]] = {}
            for node in compute_layout(family, style).nodes:
                rows.setdefault(round(node.y, 6), []).append(node.x)
            for xs in rows.values():
                xs.sort()
                for a, b in zip(xs, xs[1:], strict=False):
                    assert b - a >= gap - 1e-9, f"{family.family_id}: nodes {b - a:.2f} apart"

    def test_tree_couples_are_centred_over_their_children(self, style):
        """Without loops, each couple sits over the middle of its children."""
        from conftest import generate_clinical_pedigree

        for family in [generate_clinical_pedigree(n, s) for n in (20, 100) for s in range(3)]:
            layout = compute_layout(family, style)
            x = {n.individual_id: n.x for n in layout.nodes}
            for union in family.unions:
                kids = [x[c] for c in union.children_ids]
                centre = (x[union.partner_ids[0]] + x[union.partner_ids[1]]) / 2
                assert centre == pytest.approx((min(kids) + max(kids)) / 2)

    def test_loops_stay_adjacent_when_lines_allow(self, style):
        """First cousins are drawn side by side, without a duplicate."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        layout = compute_layout(ped.get_family("FAM1"), style)
        assert not [n for n in layout.nodes if n.is_duplicate]

    def test_tangled_loops_are_broken_with_duplicates(self, style):
        """Where adjacency would tangle lines, a partner is drawn twice."""
        from conftest import generate_realistic_pedigree

        family = generate_realistic_pedigree(100, 42)
        layout = compute_layout(family, style)
        duplicates = [n for n in layout.nodes if n.is_duplicate]
        assert duplicates
        assert len(layout.duplicate_arcs) == len(duplicates)
        labels = {n.individual_id: n.label for n in layout.nodes}
        for node in duplicates:
            assert node.label in family.individuals, "a duplicate is labelled as its original"
        assert len(labels) == len(family.individuals) + len(duplicates)

    def test_loop_breaking_leaves_the_input_family_alone(self, style):
        import copy

        from conftest import generate_realistic_pedigree

        family = generate_realistic_pedigree(100, 42)
        before = copy.deepcopy(family)
        compute_layout(family, style)
        assert family == before

    def test_duplicate_arcs_stay_between_two_rows(self, style):
        from conftest import generate_realistic_pedigree

        for arc in compute_layout(generate_realistic_pedigree(150, 7), style).duplicate_arcs:
            length = math.dist(arc.from_pos, arc.to_pos)
            bulge = abs(arc.curvature) * length / 2
            assert bulge <= 0.4 * style.generation_gap + 1e-9


class TestCrossGenerationPartners:
    def test_duplicate_stands_on_the_original_generation(self, style):
        """multiple_partners_complex: the father's second partnership is
        drawn on his own generation, not among the founders."""
        ped = parse_pedigree(TEST_DATA_DIR / "multi_partner/multiple_partners_complex.ped")
        layout = compute_layout(ped.get_family("FAM1"), style)
        y = {n.individual_id: n.y for n in layout.nodes}
        dup = next(n.individual_id for n in layout.nodes if n.is_duplicate)
        assert y[dup] == pytest.approx(y["father"])
        assert y["child3"] == pytest.approx(y["child1"])

    def test_line_reaches_a_child_drawn_a_generation_lower(self, style):
        """consanguinity_uncle_niece: the uncle stands on his niece's row;
        the line from his parents must reach him there."""
        ped = parse_pedigree(TEST_DATA_DIR / "consanguinity/consanguinity_uncle_niece.ped")
        layout = compute_layout(ped.get_family("FAM1"), style)
        uncle = next(n for n in layout.nodes if n.individual_id == "uncle")
        stub = next(
            e
            for e in layout.edges
            if e.edge_type is EdgeType.INDIVIDUAL and e.end[0] == pytest.approx(uncle.x)
        )
        assert stub.end[1] == pytest.approx(uncle.y + style.shape_size / 2)


def _added_duplicates(family, style) -> int:
    before = sum(ind.is_duplicate for ind in family.individuals.values())
    return sum(n.is_duplicate for n in compute_layout(family, style).nodes) - before


def _family_from_rows(rows: str):
    from pypedigree.parser import parse_pedigree_from_string

    text = "\n".join("\t".join(r.split()) for r in rows.strip().splitlines()) + "\n"
    return parse_pedigree_from_string(text).get_family("F")


class TestLoopDuplicates:
    """A loop is drawn without duplicating anyone whenever its branches can meet.

    Duplicates are allowed by the nomenclature but make a pedigree harder to
    read, so they are a last resort: the branches of a loop are ordered to
    meet (each ancestor on one side placed last among their siblings, on the
    other side first), and when lines still overlap, one loop per tangle is
    broken before laying out again.
    """

    def test_cousins_from_non_adjacent_branches_meet(self, style):
        # p1's and p3's children marry; p2's family sat between them.
        family = _family_from_rows(
            """
            F g1 0 0 1 1
            F g2 0 0 2 1
            F p1 g1 g2 1 1
            F s1 0 0 2 1
            F p2 g1 g2 2 1
            F s2 0 0 1 1
            F p3 g1 g2 1 1
            F s3 0 0 2 1
            F a1 p1 s1 1 1
            F a2 p1 s1 2 1
            F a3 p1 s1 2 1
            F b1 s2 p2 2 1
            F b2 s2 p2 1 1
            F d1 p3 s3 1 1
            F d2 p3 s3 2 1
            F c1 a1 d2 1 1
            """
        )
        assert _added_duplicates(family, style) == 0
        assert _overlapping_sibships(compute_layout(family, style)) == []

    def test_two_marriages_between_two_sibships_need_one_duplicate(self, style):
        # Both couples cannot sit at the single boundary between the two
        # sibships, so one partner must be drawn twice, but only one.
        family = _family_from_rows(
            """
            F g1 0 0 1 1
            F g2 0 0 2 1
            F p1 g1 g2 1 1
            F s1 0 0 2 1
            F p2 g1 g2 2 1
            F s2 0 0 1 1
            F a1 p1 s1 1 1
            F a2 p1 s1 2 1
            F b1 s2 p2 2 1
            F b2 s2 p2 1 1
            F c1 a1 b1 1 1
            F c2 b2 a2 2 1
            """
        )
        assert _added_duplicates(family, style) == 1
        assert _overlapping_sibships(compute_layout(family, style)) == []

    def test_loop_duplicate_keeps_phenotype_values(self, style):
        """The partner drawn twice to break a loop gets the original's bands."""
        from pypedigree.parser import parse_pedigree_from_string

        rows = [
            "g1 0 0 1",
            "g2 0 0 2",
            "p1 g1 g2 1",
            "s1 0 0 2",
            "p2 g1 g2 2",
            "s2 0 0 1",
            "a1 p1 s1 1",
            "a2 p1 s1 2",
            "b1 s2 p2 2",
            "b2 s2 p2 1",
            "c1 a1 b1 1",
            "c2 b2 a2 2",
        ]
        lines = ["FID\tIID\tFather\tMother\tSex\tg_IQ"]
        lines += ["F\t" + "\t".join(r.split()) + f"\t{90 + i}" for i, r in enumerate(rows)]
        family = parse_pedigree_from_string("\n".join(lines) + "\n").get_family("F")
        nodes = compute_layout(family, style).nodes
        dups = [n for n in nodes if n.is_duplicate]
        assert len(dups) == 1
        original = next(n for n in nodes if n.individual_id == dups[0].label)
        assert dups[0].extra_phenotypes_raw == original.extra_phenotypes_raw != {}
        assert dups[0].extra_phenotypes == original.extra_phenotypes

    def test_one_cousin_marriage_in_a_clinical_tree_needs_no_duplicate(self, style):
        """Was duplicated in 7 of these 18 pedigrees, however distant the cousins."""
        from conftest import generate_clinical_pedigree

        for n in (30, 60, 120):
            for seed in range(6):
                family = generate_clinical_pedigree(n, seed, cousin_marriages=1)
                assert _added_duplicates(family, style) == 0, (n, seed)

    def test_several_cousin_marriages_stay_below_the_recorded_count(self, style):
        """Regression guard: 57 duplicates over these 24 pedigrees before the
        loop-aware ordering, 35 after. Lower is fine; higher is a regression."""
        from conftest import generate_clinical_pedigree

        total = sum(
            _added_duplicates(generate_clinical_pedigree(n, seed, cousin_marriages=3), style)
            for n in (30, 60, 120, 200)
            for seed in range(6)
        )
        assert total <= 35
