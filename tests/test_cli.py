"""Tests for pypedigree.cli."""

import shutil

import pytest
from click.testing import CliRunner
from conftest import TEST_DATA_DIR

from pypedigree.cli import main


@pytest.fixture
def runner():
    return CliRunner()


class TestCLIBasic:
    def test_simple_file_produces_output(self, runner, tmp_path):
        out = tmp_path / "result.png"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert out.stat().st_size > 0

    def test_default_output_name(self, runner, tmp_path):
        src = TEST_DATA_DIR / "basic/simple_trio.ped"
        dst = tmp_path / "simple_trio.ped"
        shutil.copy(src, dst)
        result = runner.invoke(main, [str(dst)])
        assert result.exit_code == 0, result.output
        expected_out = tmp_path / "simple_trio.png"
        assert expected_out.exists()

    def test_svg_format(self, runner, tmp_path):
        out = tmp_path / "result.svg"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_pdf_format_via_flag(self, runner, tmp_path):
        out = tmp_path / "result.pdf"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
                "--format",
                "pdf",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()


class TestCLIOptions:
    def test_family_filter(self, runner, tmp_path):
        out = tmp_path / "fam1.png"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/multi_family.ped"),
                "-o",
                str(out),
                "--family",
                "FAM1",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_no_labels(self, runner, tmp_path):
        out = tmp_path / "nolabels.png"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
                "--no-labels",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_dpi_flag(self, runner, tmp_path):
        out = tmp_path / "highres.png"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
                "--dpi",
                "300",
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()


class TestCLIErrors:
    def test_nonexistent_file(self, runner):
        result = runner.invoke(main, ["/nonexistent/file.ped"])
        assert result.exit_code != 0

    def test_invalid_family_id(self, runner, tmp_path):
        out = tmp_path / "bad.png"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
                "--family",
                "NONEXISTENT",
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower()

    def test_version_flag(self, runner):
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "pypedigree" in result.output.lower()

    def test_help_flag(self, runner):
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "pedigree" in result.output.lower()


class TestCLIEdgeCases:
    def test_style_file_overrides(self, runner, tmp_path):
        """--style flag loads a JSON style file."""
        import json

        style_file = tmp_path / "style.json"
        style_file.write_text(json.dumps({"dpi": 72}))
        out = tmp_path / "styled.png"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
                "--style",
                str(style_file),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_format_suffix_override(self, runner, tmp_path):
        """--format svg with .png output path should produce .svg suffix."""
        out = tmp_path / "result.png"
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "basic/simple_trio.ped"),
                "-o",
                str(out),
                "--format",
                "svg",
            ],
        )
        assert result.exit_code == 0, result.output
        svg_out = tmp_path / "result.svg"
        assert svg_out.exists()

    def test_malformed_pedigree_error(self, runner, tmp_path):
        """Malformed pedigree file should produce an error exit code."""
        result = runner.invoke(
            main,
            [
                str(TEST_DATA_DIR / "format/malformed_missing_cols.ped"),
                "-o",
                str(tmp_path / "out.png"),
            ],
        )
        assert result.exit_code != 0

    def test_main_module_entrypoint(self):
        """python -m pypedigree --help should work.

        Runs the interpreter of the current test session. Going through
        `uv run` instead re-synced the project's .venv to the lockfile in the
        middle of the session, which broke a run against other dependency
        versions (pandas and pytest disappeared for every later test) and
        made the suite depend on uv being installed.
        """
        import subprocess
        import sys

        result = subprocess.run(
            [sys.executable, "-m", "pypedigree", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "pedigree" in result.stdout.lower()


class TestCLIWarnings:
    """Validation warnings are reported as CLI messages, not raw Python warnings."""

    @pytest.fixture
    def sex_mismatch_ped(self, tmp_path):
        """A father recorded as female, shared by two children."""
        path = tmp_path / "sex_mismatch.ped"
        path.write_text(
            "FAM1\tdad\t0\t0\t2\t1\n"
            "FAM1\tmum\t0\t0\t2\t1\n"
            "FAM1\tkid1\tdad\tmum\t1\t2\n"
            "FAM1\tkid2\tdad\tmum\t2\t1\n"
        )
        return path

    def test_warning_is_reported_cleanly(self, runner, tmp_path, sex_mismatch_ped):
        result = runner.invoke(main, [str(sex_mismatch_ped), "-o", str(tmp_path / "out.png")])
        assert result.exit_code == 0, result.output
        assert "Warning: Individual 'dad' listed as father but has sex=FEMALE" in result.output

    def test_warning_carries_no_internal_paths(self, runner, tmp_path, sex_mismatch_ped):
        result = runner.invoke(main, [str(sex_mismatch_ped), "-o", str(tmp_path / "out.png")])
        assert "UserWarning" not in result.output
        assert "parser.py" not in result.output

    def test_warning_reported_once_per_parent(self, runner, tmp_path, sex_mismatch_ped):
        """Both children share one father: the warning must not repeat."""
        result = runner.invoke(main, [str(sex_mismatch_ped), "-o", str(tmp_path / "out.png")])
        assert result.output.count("listed as father but has sex=FEMALE") == 1
