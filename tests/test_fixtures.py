"""Tests for the test fixtures themselves.

Two things are asserted here. First, the generated fixtures on disk match
what `generate_fixtures.py` produces, so the claim that they are synthetic
stays verifiable and a hand-edit cannot drift away unnoticed. Second, no
fixture carries an identifier shaped like a cohort accession, which is the
kind of value that must never reach a public repository.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import TEST_DATA_DIR

GENERATOR = TEST_DATA_DIR / "generate_fixtures.py"

# Files written by the generator, relative to tests/data/.
GENERATED = [
    "complex/synesthesia.ped",
    "complex/large_family.ped",
    "complex/multi_phenotype_pgs.ped",
    "phenotype/test_all.ped",
    "phenotype/annotations.ped",
    "complex/multigenerational.ped",
]

# A letter followed by digit groups joined by hyphens: the shape of a study
# accession or sample barcode (e.g. "X0000-00-000").
ACCESSION_RE = re.compile(r"[A-Za-z][0-9]{3,}-[0-9]{2,}-[0-9]{2,}")


class TestGeneratedFixturesAreReproducible:
    def test_generator_exists_and_is_documented(self):
        assert GENERATOR.exists()
        assert (TEST_DATA_DIR / "README.md").exists()

    def test_regenerating_changes_nothing(self, tmp_path):
        """Running the generator must reproduce the committed files exactly.

        This is what makes "these fixtures are synthetic" checkable rather
        than a claim in a comment.
        """
        before = {rel: (TEST_DATA_DIR / rel).read_bytes() for rel in GENERATED}

        result = subprocess.run(
            [sys.executable, str(GENERATOR)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr

        for rel, original in before.items():
            assert (TEST_DATA_DIR / rel).read_bytes() == original, (
                f"{rel} differs from generator output — regenerate it, or revert the hand-edit"
            )

    @pytest.mark.parametrize("rel", GENERATED)
    def test_generated_fixture_has_a_paired_image(self, rel):
        assert (TEST_DATA_DIR / rel).with_suffix(".png").exists()


class TestNoIdentifiersInFixtures:
    """No fixture may carry a value shaped like a study accession."""

    @pytest.mark.parametrize(
        "path",
        sorted(TEST_DATA_DIR.rglob("*.ped"))
        + sorted(TEST_DATA_DIR.rglob("*.fam"))
        + sorted(TEST_DATA_DIR.rglob("*.tsv"))
        + sorted(TEST_DATA_DIR.rglob("*.csv")),
        ids=lambda p: str(p.relative_to(TEST_DATA_DIR)),
    )
    def test_no_accession_shaped_identifier(self, path: Path):
        matches = ACCESSION_RE.findall(path.read_text(errors="replace"))
        assert not matches, f"{path.name} contains accession-shaped value(s): {matches}"

    def test_family_ids_are_from_the_synthetic_set(self):
        """Family IDs are invented labels, not numbers carried over from data."""
        allowed = {"FID", "FAM1", "FAM2", "FAM_SYN01"}
        seen: set[str] = set()
        for path in sorted(TEST_DATA_DIR.rglob("*.ped")):
            lines = path.read_text(errors="replace").splitlines()
            for line in lines[1:]:
                if line.strip():
                    seen.add(line.split("\t")[0].split(",")[0].split()[0])
        assert seen <= allowed, f"unexpected family IDs: {sorted(seen - allowed)}"

    def test_no_dates_of_birth_are_populated(self):
        """A real DOB is directly identifying; fixtures leave the column at 0."""
        for path in sorted(TEST_DATA_DIR.rglob("*.ped")):
            lines = path.read_text(errors="replace").splitlines()
            if not lines:
                continue
            header = lines[0].split("\t")
            if "DOB" not in header:
                continue
            i = header.index("DOB")
            for line in lines[1:]:
                fields = line.split("\t")
                if len(fields) > i and line.strip():
                    assert fields[i] in ("0", "-9", ""), (
                        f"{path.name} has a populated DOB: {fields[i]}"
                    )
