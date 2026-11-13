"""Visual regression tests for pypedigree rendering.

Compares rendered output against golden reference images in
tests/visual/golden/, by counting the pixels that changed. A mean
difference over the whole image, used before, let real regressions
through: a descent line drawn diagonally, or a wrong legend, changes a few
hundred pixels of a mostly white image and moved the mean by less than
0.001.

To regenerate golden images (e.g. after an intentional visual change):
    pytest tests/visual/test_visual_regression.py --regen-golden
"""

from __future__ import annotations

import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pypedigree.config import PedigreeStyle
from pypedigree.parser import apply_phenotypes_to_pedigree, parse_pedigree, parse_phenotype_file
from pypedigree.renderer import render_pedigree

GOLDEN_DIR = Path(__file__).parent / "golden"
DATA_DIR = Path(__file__).parent.parent / "data"

# A pixel counts as changed when one of its channels moves by more than this
# (on a 0-1 scale): far above anti-aliasing jitter, far below any drawn line.
MAX_CHANNEL_DIFF = 0.1

# A render fails when more than this fraction of its pixels changed. Renders
# are reproducible to the pixel across platforms and the supported matplotlib
# versions (3.7 to 3.10 measured: 0 pixels); the smallest real regression
# seen, a legend change, moved 597 pixels (0.1%).
MAX_CHANGED_FRACTION = 1e-4


@pytest.fixture
def regen_golden(request: pytest.FixtureRequest) -> bool:
    return request.config.getoption("--regen-golden", default=False)


def _render_to_array(fig: plt.Figure, dpi: int = 100) -> np.ndarray:
    """Render a matplotlib figure to a numpy RGBA array."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    img = plt.imread(buf)
    plt.close(fig)
    return img


def _changed_pixels(actual: np.ndarray, expected: np.ndarray) -> np.ndarray | None:
    """Mask of the pixels that differ, or None when the image sizes differ.

    A size change is a layout change (renders are cropped to their content),
    so it is not tolerated.
    """
    if actual.shape != expected.shape:
        return None
    diff = np.abs(actual.astype(float) - expected.astype(float)).max(axis=2)
    return np.asarray(diff > MAX_CHANNEL_DIFF)


def _save_diff(actual: np.ndarray, expected: np.ndarray, mask: np.ndarray, path: Path) -> None:
    """Save the expected image greyed out, changed pixels in red."""
    grey = expected[..., :3].mean(axis=2, keepdims=True)
    out = np.concatenate([np.repeat(0.6 + 0.4 * grey, 3, axis=2), np.ones(grey.shape)], axis=2)
    out[mask] = (1.0, 0.0, 0.0, 1.0)
    plt.imsave(str(path), out)


def _golden_path(name: str) -> Path:
    return GOLDEN_DIR / f"{name}.png"


def _assert_visual_match(
    fig: plt.Figure,
    name: str,
    regen: bool,
    dpi: int = 100,
    max_changed: float = MAX_CHANGED_FRACTION,
) -> None:
    """Assert that a rendered figure matches its golden reference."""
    actual = _render_to_array(fig, dpi=dpi)
    golden = _golden_path(name)

    if regen:
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        plt.imsave(str(golden), actual)
        return

    if not golden.exists():
        # First run: generate golden automatically
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        plt.imsave(str(golden), actual)
        pytest.skip(f"Golden image generated: {golden}. Re-run to compare.")

    expected = plt.imread(str(golden))
    mask = _changed_pixels(actual, expected)
    fail_path = golden.with_suffix(".actual.png")
    if mask is None:
        plt.imsave(str(fail_path), actual)
        pytest.fail(
            f"Visual regression for '{name}': size {actual.shape[1]}x{actual.shape[0]} "
            f"instead of {expected.shape[1]}x{expected.shape[0]}. Actual saved to {fail_path}"
        )
    changed = int(mask.sum())
    if changed > max_changed * mask.size:
        plt.imsave(str(fail_path), actual)
        diff_path = golden.with_suffix(".diff.png")
        _save_diff(actual, expected, mask, diff_path)
        rows, cols = np.nonzero(mask)
        pytest.fail(
            f"Visual regression for '{name}': {changed} pixels changed "
            f"({changed / mask.size:.3%}, limit {max_changed:.3%}), in the box "
            f"x={cols.min()}-{cols.max()}, y={rows.min()}-{rows.max()}. "
            f"Actual saved to {fail_path}, changes in red in {diff_path}"
        )


@pytest.fixture
def style() -> PedigreeStyle:
    return PedigreeStyle()


class TestVisualSimple:
    """Visual regression tests for basic pedigree rendering."""

    def test_simple_trio(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "basic/simple_trio.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "simple_trio", regen_golden)

    def test_three_gen(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "basic/three_gen.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "three_gen", regen_golden)

    def test_two_gen_siblings(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "basic/two_gen_siblings.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "two_gen_siblings", regen_golden)


class TestVisualExtended:
    """Visual regression for extended features (deceased, proband, twins)."""

    def test_deceased(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "extended/extended_deceased.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "deceased", regen_golden)

    def test_proband(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "extended/extended_proband.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "proband", regen_golden)

    def test_twins_mz(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "twins/twins_mz.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "twins_mz", regen_golden)

    def test_twins_dz(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "twins/twins_dz.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "twins_dz", regen_golden)

    def test_consanguinity(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "consanguinity/consanguinity_first_cousins.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "consanguinity", regen_golden)

    def test_multigenerational_loops(self, style: PedigreeStyle, regen_golden: bool) -> None:
        """Five generations, two loops drawn side by side, bands and carriers."""
        ped = parse_pedigree(DATA_DIR / "complex/multigenerational.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "multigenerational", regen_golden)

    def test_nomenclature(self, style: PedigreeStyle, regen_golden: bool) -> None:
        """Pregnancies, losses, adoption, childlessness, donor, gender, consultand."""
        ped = parse_pedigree(DATA_DIR / "nomenclature/nomenclature.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "nomenclature", regen_golden)


class TestVisualMultiPhenotype:
    """Visual regression for multi-phenotype rendering."""

    def test_inline_phenotypes(self, regen_golden: bool) -> None:
        style = PedigreeStyle()
        style.phenotype_colors = {
            "epilepsy": "#e74c3c",
            "IQ": "#3498db",
            "SRS": "#e67e22",
        }
        ped = parse_pedigree(DATA_DIR / "phenotype/multi_phenotype_inline.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "multi_phenotype_inline", regen_golden)

    def test_inline_phenotypes_with_pgs(self, regen_golden: bool) -> None:
        style = PedigreeStyle()
        style.phenotype_colors = {
            "epilepsy": "#e74c3c",
            "IQ": "#3498db",
            "SRS": "#e67e22",
        }
        # render_pedigree widens the gaps for the PGS curves itself.
        ped = parse_pedigree(DATA_DIR / "phenotype/multi_phenotype_pgs.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "multi_phenotype_pgs", regen_golden)

    def test_pheno_file(self, regen_golden: bool) -> None:
        style = PedigreeStyle()
        style.phenotype_colors = {
            "epilepsy": "#e74c3c",
            "IQ": "#3498db",
            "SRS": "#e67e22",
        }
        ped = parse_pedigree(DATA_DIR / "basic/simple_trio.ped")
        pheno_data, pheno_names = parse_phenotype_file(DATA_DIR / "phenotype/test_pheno.tsv")
        apply_phenotypes_to_pedigree(ped, pheno_data, pheno_names)
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "pheno_file", regen_golden)

    def test_no_legend(self, regen_golden: bool) -> None:
        style = PedigreeStyle()
        style.phenotype_colors = {
            "epilepsy": "#e74c3c",
            "IQ": "#3498db",
            "SRS": "#e67e22",
        }
        ped = parse_pedigree(DATA_DIR / "phenotype/multi_phenotype_inline.ped")
        fig = render_pedigree(ped, style=style, show_legend=False)
        _assert_visual_match(fig, "multi_phenotype_no_legend", regen_golden)

    def test_no_labels(self, regen_golden: bool) -> None:
        style = PedigreeStyle()
        style.phenotype_colors = {
            "epilepsy": "#e74c3c",
            "IQ": "#3498db",
            "SRS": "#e67e22",
        }
        ped = parse_pedigree(DATA_DIR / "phenotype/multi_phenotype_inline.ped")
        fig = render_pedigree(ped, style=style, show_labels=False)
        _assert_visual_match(fig, "multi_phenotype_no_labels", regen_golden)


class TestVisualLayout:
    """Layout features added in 0.4.0: sibships, duplicates, loops."""

    def test_separate_sibships(self, style: PedigreeStyle, regen_golden: bool) -> None:
        """The 11-person tree whose two sibships used to share one line."""
        from pypedigree.parser import parse_pedigree_from_string

        text = "".join(
            f"F\t{i}\t{f}\t{m}\t{s}\t1\n"
            for i, f, m, s in [
                ("f0", 0, 0, 1),
                ("f1", 0, 0, 2),
                ("i2", "f0", "f1", 1),
                ("i3", "f0", "f1", 1),
                ("s4", 0, 0, 2),
                ("s5", 0, 0, 2),
                ("i6", "i3", "s5", 2),
                ("i7", "i3", "s5", 1),
                ("s8", 0, 0, 1),
                ("s9", 0, 0, 2),
                ("i10", "i2", "s4", 1),
            ]
        )
        fig = render_pedigree(parse_pedigree_from_string(text), style=style)
        _assert_visual_match(fig, "separate_sibships", regen_golden)

    def test_second_partnership(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "multi_partner/multiple_partners_complex.ped")
        fig = render_pedigree(ped, style=style)
        _assert_visual_match(fig, "multiple_partners_complex", regen_golden)

    def test_uncle_niece(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "consanguinity/consanguinity_uncle_niece.ped")
        _assert_visual_match(render_pedigree(ped, style=style), "uncle_niece", regen_golden)


class TestVisualCarriers:
    """NSGC 2022 carrier notation, and the legacy notation."""

    def test_carriers(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "phenotype/carrier.ped")
        _assert_visual_match(render_pedigree(ped, style=style), "carrier", regen_golden)

    def test_carrier_genes(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "phenotype/carrier_genes.ped")
        _assert_visual_match(render_pedigree(ped, style=style), "carrier_genes", regen_golden)

    def test_legacy_notation(self, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "phenotype/carrier.ped")
        style = PedigreeStyle(carrier_notation="legacy")
        _assert_visual_match(render_pedigree(ped, style=style), "carrier_legacy", regen_golden)


class TestVisualCircular:
    """Circular layout: rings, twins, duplicates, bands and labels."""

    def test_three_gen(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "basic/three_gen.ped")
        fig = render_pedigree(ped, style=style, circular=True)
        _assert_visual_match(fig, "circular_three_gen", regen_golden)

    def test_twins(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "twins/twins_with_siblings.ped")
        fig = render_pedigree(ped, style=style, circular=True)
        _assert_visual_match(fig, "circular_twins", regen_golden)

    def test_second_partnership(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "multi_partner/multiple_partners_complex.ped")
        fig = render_pedigree(ped, style=style, circular=True)
        _assert_visual_match(fig, "circular_multiple_partners", regen_golden)

    def test_bands_and_pgs(self, style: PedigreeStyle, regen_golden: bool) -> None:
        ped = parse_pedigree(DATA_DIR / "phenotype/test_all.ped")
        fig = render_pedigree(ped, style=style, circular=True)
        _assert_visual_match(fig, "circular_test_all", regen_golden)
