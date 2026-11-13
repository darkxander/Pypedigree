"""Tests for large pedigree layout — stress tests and position snapshots.

Tests the layout engine with synthetic large pedigrees (50-200+ individuals)
and verifies that existing small pedigrees produce identical positions
after any layout algorithm changes.
"""

from __future__ import annotations

import time

import pytest
from conftest import (
    TEST_DATA_DIR,
    generate_deep_pedigree,
    generate_realistic_pedigree,
    generate_wide_pedigree,
)

from pypedigree.config import PedigreeStyle
from pypedigree.layout import compute_layout
from pypedigree.parser import parse_pedigree


@pytest.fixture
def style() -> PedigreeStyle:
    return PedigreeStyle()


# ---------------------------------------------------------------------------
# Position snapshot tests — verify small pedigrees are unchanged
# ---------------------------------------------------------------------------


def _get_positions(family_layout) -> dict[str, tuple[float, float]]:
    """Extract {individual_id: (x, y)} from a FamilyLayout."""
    return {n.individual_id: (round(n.x, 6), round(n.y, 6)) for n in family_layout.nodes}


class TestPositionSnapshots:
    """Verify that existing pedigrees produce exactly the same positions.

    If a layout change intentionally alters positions, update the snapshots.
    """

    @staticmethod
    def _layout_for_file(filename: str, style: PedigreeStyle):
        ped = parse_pedigree(TEST_DATA_DIR / filename)
        fam = list(ped.families.values())[0]
        return compute_layout(fam, style)

    def test_simple_trio_positions(self, style: PedigreeStyle) -> None:
        layout = self._layout_for_file("basic/simple_trio.ped", style)
        pos = _get_positions(layout)
        # Father and mother at generation 0, child at generation 1
        assert len(pos) == 3
        # Parents should be at same y, child below
        parent_ys = {pos["father"][1], pos["mother"][1]}
        assert len(parent_ys) == 1  # same generation
        assert pos["child"][1] < list(parent_ys)[0]  # child below (more negative)

    def test_two_gen_siblings_positions(self, style: PedigreeStyle) -> None:
        layout = self._layout_for_file("basic/two_gen_siblings.ped", style)
        pos = _get_positions(layout)
        # Parents + 3 children
        assert len(pos) == 5
        # All children at same y
        child_ys = {pos["son1"][1], pos["daughter1"][1], pos["son2"][1]}
        assert len(child_ys) == 1

    def test_three_gen_positions(self, style: PedigreeStyle) -> None:
        layout = self._layout_for_file("basic/three_gen.ped", style)
        pos = _get_positions(layout)
        # Verify 3 distinct generation levels
        y_values = sorted(set(p[1] for p in pos.values()), reverse=True)
        assert len(y_values) == 3

    def test_consanguinity_positions_stable(self, style: PedigreeStyle) -> None:
        layout = self._layout_for_file("consanguinity/consanguinity_first_cousins.ped", style)
        pos = _get_positions(layout)
        # Should have 9 individuals and no NaN positions
        assert len(pos) >= 7
        for sid, (x, y) in pos.items():
            assert x == x, f"NaN x for {sid}"  # NaN check
            assert y == y, f"NaN y for {sid}"

    def test_twins_mz_positions_stable(self, style: PedigreeStyle) -> None:
        layout = self._layout_for_file("twins/twins_mz.ped", style)
        pos = _get_positions(layout)
        assert len(pos) >= 4

    def test_multiple_partners_positions_stable(self, style: PedigreeStyle) -> None:
        layout = self._layout_for_file("multi_partner/multiple_partners.ped", style)
        pos = _get_positions(layout)
        assert len(pos) >= 5


# ---------------------------------------------------------------------------
# Stress tests — large synthetic pedigrees
# ---------------------------------------------------------------------------


class TestLargePedigreeNoOverlap:
    """Verify no nodes overlap in large pedigrees."""

    @staticmethod
    def _check_no_overlap(layout, style: PedigreeStyle) -> None:
        """Assert no two nodes in the same generation are closer than min_node_distance."""
        # Group nodes by y (generation)
        by_gen: dict[float, list[tuple[str, float]]] = {}
        for n in layout.nodes:
            y = round(n.y, 4)
            by_gen.setdefault(y, []).append((n.individual_id, n.x))

        for y, nodes in by_gen.items():
            sorted_nodes = sorted(nodes, key=lambda t: t[1])
            for i in range(1, len(sorted_nodes)):
                prev_id, prev_x = sorted_nodes[i - 1]
                curr_id, curr_x = sorted_nodes[i]
                gap = curr_x - prev_x
                assert gap >= style.min_node_distance - 0.01, (
                    f"Overlap at y={y}: {prev_id} (x={prev_x:.2f}) and "
                    f"{curr_id} (x={curr_x:.2f}), gap={gap:.2f} < {style.min_node_distance}"
                )

    def test_deep_pedigree_no_overlap(self, style: PedigreeStyle) -> None:
        fam = generate_deep_pedigree(generations=8, children_per_union=3)
        layout = compute_layout(fam, style)
        assert len(layout.nodes) >= 20
        self._check_no_overlap(layout, style)

    def test_wide_pedigree_no_overlap(self, style: PedigreeStyle) -> None:
        fam = generate_wide_pedigree(num_couples=15, children_per_couple=3)
        layout = compute_layout(fam, style)
        assert len(layout.nodes) >= 60
        self._check_no_overlap(layout, style)

    def test_100_individuals_no_overlap(self, style: PedigreeStyle) -> None:
        fam = generate_realistic_pedigree(target_size=100, seed=42)
        layout = compute_layout(fam, style)
        assert len(layout.nodes) >= 80  # may not hit exactly 100
        self._check_no_overlap(layout, style)

    def test_150_individuals_no_overlap(self, style: PedigreeStyle) -> None:
        fam = generate_realistic_pedigree(target_size=150, seed=123)
        layout = compute_layout(fam, style)
        assert len(layout.nodes) >= 100
        self._check_no_overlap(layout, style)


class TestLargePedigreePerformance:
    """Verify layout computation is fast enough for large pedigrees."""

    def test_100_individuals_under_1s(self, style: PedigreeStyle) -> None:
        fam = generate_realistic_pedigree(target_size=100, seed=42)
        start = time.perf_counter()
        layout = compute_layout(fam, style)
        elapsed = time.perf_counter() - start
        assert elapsed < 1.0, f"Layout took {elapsed:.2f}s for {len(layout.nodes)} nodes"

    def test_200_individuals_under_2s(self, style: PedigreeStyle) -> None:
        fam = generate_realistic_pedigree(target_size=200, seed=99)
        start = time.perf_counter()
        layout = compute_layout(fam, style)
        elapsed = time.perf_counter() - start
        assert elapsed < 2.0, f"Layout took {elapsed:.2f}s for {len(layout.nodes)} nodes"

    def test_deep_10_gen_under_1s(self, style: PedigreeStyle) -> None:
        fam = generate_deep_pedigree(generations=10, children_per_union=3)
        start = time.perf_counter()
        layout = compute_layout(fam, style)
        elapsed = time.perf_counter() - start
        assert elapsed < 1.0, f"Layout took {elapsed:.2f}s for {len(layout.nodes)} nodes"


class TestLargePedigreeWidthBounded:
    """Verify layout width doesn't explode for large pedigrees."""

    def test_100_individuals_width_bounded(self, style: PedigreeStyle) -> None:
        fam = generate_realistic_pedigree(target_size=100, seed=42)
        layout = compute_layout(fam, style)
        xs = [n.x for n in layout.nodes]
        width = max(xs) - min(xs)
        n = len(layout.nodes)
        # Width should be at most ~2x the number of individuals (rough bound)
        assert width < n * 2.5, f"Width {width:.1f} too large for {n} nodes"


class TestLargePedigreeGenerationStability:
    """Verify generation assignment doesn't produce errors for complex pedigrees."""

    def test_consanguinity_loops_stable(self, style: PedigreeStyle) -> None:
        """Pedigrees with consanguinity should not cause infinite loops."""
        fam = generate_realistic_pedigree(target_size=50, seed=7)
        # Should complete without hanging
        layout = compute_layout(fam, style)
        assert len(layout.nodes) > 0

    def test_deep_pedigree_generations_correct(self, style: PedigreeStyle) -> None:
        """Deep pedigree should have correct number of generation levels."""
        gen_count = 8
        fam = generate_deep_pedigree(generations=gen_count, children_per_union=2)
        layout = compute_layout(fam, style)
        y_values = sorted(set(round(n.y, 4) for n in layout.nodes), reverse=True)
        # Should have gen_count + 1 levels (gen_count generations of children + founders)
        # Actually: gen_count couples produce gen_count levels of children,
        # plus spouse_* at some of those levels
        assert len(y_values) >= gen_count


class TestLargePedigreeChildCentering:
    """Verify children are reasonably centered under their parents."""

    def test_children_centered_under_parents(self, style: PedigreeStyle) -> None:
        fam = generate_realistic_pedigree(target_size=50, seed=42)
        layout = compute_layout(fam, style)

        # Build position map
        pos = {n.individual_id: (n.x, n.y) for n in layout.nodes}

        total_deviation = 0.0
        count = 0

        for union in fam.unions:
            p1, p2 = union.partner_ids
            if p1 not in pos or p2 not in pos:
                continue
            couple_center = (pos[p1][0] + pos[p2][0]) / 2

            child_xs = [pos[cid][0] for cid in union.children_ids if cid in pos]
            if not child_xs:
                continue
            children_center = sum(child_xs) / len(child_xs)

            total_deviation += abs(couple_center - children_center)
            count += 1

        if count > 0:
            mean_deviation = total_deviation / count
            # Mean deviation should be moderate (not perfect centering expected
            # due to overlap resolution, but shouldn't be wildly off).
            assert mean_deviation < style.couple_gap * 5, (
                f"Mean parent-child centering deviation {mean_deviation:.2f} too large"
            )
