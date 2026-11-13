"""Tests for pypedigree.config."""

from dataclasses import fields

import pytest

from pypedigree.config import PedigreeStyle
from pypedigree.errors import PedigreeStyleWarning


class TestPedigreeStyleDefaults:
    def test_default_shape_size(self):
        s = PedigreeStyle()
        assert s.shape_size == 0.8

    def test_default_line_color(self):
        s = PedigreeStyle()
        assert s.line_color == "black"

    def test_default_dpi(self):
        s = PedigreeStyle()
        assert s.dpi == 150

    def test_all_defaults_are_set(self):
        s = PedigreeStyle()
        for f in fields(s):
            assert getattr(s, f.name) is not None


class TestFromDict:
    def test_partial_override(self):
        s = PedigreeStyle.from_dict({"shape_size": 1.2, "dpi": 300})
        assert s.shape_size == 1.2
        assert s.dpi == 300
        assert s.line_color == "black"

    def test_unknown_keys_ignored(self):
        with pytest.warns(PedigreeStyleWarning):
            s = PedigreeStyle.from_dict({"nonexistent_field": 42, "shape_size": 1.0})
        assert s.shape_size == 1.0
        assert not hasattr(s, "nonexistent_field")

    def test_unknown_key_warns_with_suggestion(self):
        with pytest.warns(PedigreeStyleWarning, match="did you mean 'phenotype_colors'"):
            PedigreeStyle.from_dict({"phenotype_color": {"IQ": "#111111"}})

    def test_unknown_key_without_close_match_still_warns(self):
        with pytest.warns(PedigreeStyleWarning, match="Unknown style option 'zzz'"):
            PedigreeStyle.from_dict({"zzz": 1})

    def test_warning_can_be_suppressed(self):
        import warnings as _warnings

        with _warnings.catch_warnings():
            _warnings.simplefilter("error", PedigreeStyleWarning)
            # warn_unknown=False must not raise even with the filter escalated
            s = PedigreeStyle.from_dict({"zzz": 1, "dpi": 300}, warn_unknown=False)
        assert s.dpi == 300

    def test_known_keys_do_not_warn(self):
        import warnings as _warnings

        with _warnings.catch_warnings():
            _warnings.simplefilter("error", PedigreeStyleWarning)
            assert PedigreeStyle.from_dict({"dpi": 300}).dpi == 300

    def test_empty_dict_gives_defaults(self):
        s = PedigreeStyle.from_dict({})
        default = PedigreeStyle()
        assert s.shape_size == default.shape_size
        assert s.dpi == default.dpi


class TestToDict:
    def test_round_trip(self):
        original = PedigreeStyle(shape_size=1.5, dpi=300)
        d = original.to_dict()
        restored = PedigreeStyle.from_dict(d)
        assert restored.shape_size == 1.5
        assert restored.dpi == 300

    def test_all_fields_in_dict(self):
        s = PedigreeStyle()
        d = s.to_dict()
        for f in fields(s):
            assert f.name in d
