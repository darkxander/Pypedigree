"""Conftest for visual regression tests."""

import matplotlib

matplotlib.use("Agg")


def pytest_addoption(parser):
    parser.addoption(
        "--regen-golden",
        action="store_true",
        default=False,
        help="Regenerate golden reference images instead of comparing.",
    )
