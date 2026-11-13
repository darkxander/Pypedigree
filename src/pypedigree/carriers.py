"""Carrier-status notation shared by the matplotlib and SVG renderers.

No matplotlib dependency: this module only decides *what* to draw. The
renderers turn each pattern name into a hatch (matplotlib) or an SVG
``<pattern>``.

Two notations are supported, selected by ``PedigreeStyle.carrier_notation``:

- ``"nsgc2022"`` (default): the NSGC focused revision of the standardized
  pedigree nomenclature (Bennett et al., J Genet Couns 2022) no longer uses
  a dot for carriers. The symbol is divided into as many sections as needed,
  each with its own fill pattern explained in the legend, so a carrier of
  two conditions, or an affected carrier, shows every result. Each carrier
  source (the ``carrier`` column, and each ``cr_`` column) gets its own
  pattern; obligate carriers get a sparser version of it.
- ``"legacy"``: the ``carrier`` column is drawn as before 0.4.0, a
  half-filled symbol for carriers and a central dot for obligate carriers.
  ``cr_`` columns are still drawn as patterned sections.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Protocol

from pypedigree.errors import PedigreeError
from pypedigree.models import PhenotypeDisplayMode

CARRIER_NOTATIONS = ("nsgc2022", "legacy")

# Key of the `carrier` column among the carrier sources. Angle brackets keep
# it apart from any `cr_` column name.
LEGACY_CARRIER_KEY = "<carrier>"

# One pattern per carrier source, in order of appearance.
PATTERNS = (
    "diagonal",
    "antidiagonal",
    "vertical",
    "horizontal",
    "crossed",
    "grid",
    "dotted",
    "circles",
)

CARRIER = 1
OBLIGATE_CARRIER = 2


class _HasCarrierData(Protocol):
    """An Individual or a NodeLayout: both carry these attributes."""

    is_carrier: bool
    is_obligate_carrier: bool
    extra_phenotypes_raw: dict[str, float | None]


@dataclass(frozen=True)
class CarrierMark:
    """One patterned section of a symbol.

    Attributes:
        key: Carrier source (``LEGACY_CARRIER_KEY`` or a ``cr_`` column name).
        pattern: Pattern name from ``PATTERNS``.
        obligate: Whether this is an obligate carrier (sparser pattern).
    """

    key: str
    pattern: str
    obligate: bool


@dataclass(frozen=True)
class CarrierScheme:
    """How carrier status is drawn in one figure.

    Attributes:
        notation: ``"nsgc2022"`` or ``"legacy"``.
        keys: Carrier sources drawn as patterned sections, in order.
        obligate_keys: Sources with at least one obligate carrier.
        has_legacy_carriers: Whether the ``carrier`` column marks anyone.
        has_legacy_plain: Whether it marks a (non-obligate) carrier.
        has_legacy_obligate: Whether it marks an obligate carrier.
    """

    notation: str = "nsgc2022"
    keys: tuple[str, ...] = ()
    obligate_keys: frozenset[str] = field(default_factory=frozenset)
    has_legacy_carriers: bool = False
    has_legacy_plain: bool = False
    has_legacy_obligate: bool = False

    @property
    def legacy(self) -> bool:
        """Whether the ``carrier`` column uses the half-fill / dot notation."""
        return self.notation == "legacy"

    def pattern(self, key: str) -> str:
        """Pattern of a carrier source."""
        return PATTERNS[self.keys.index(key) % len(PATTERNS)]

    def marks(self, who: _HasCarrierData) -> list[CarrierMark]:
        """Patterned sections for one individual, in source order."""
        out = []
        for key in self.keys:
            state = carrier_state(who, key)
            if state in (CARRIER, OBLIGATE_CARRIER):
                out.append(CarrierMark(key, self.pattern(key), state == OBLIGATE_CARRIER))
        return out

    def label(self, key: str, obligate: bool) -> str:
        """Legend text of a carrier source."""
        kind = "Obligate carrier" if obligate else "Carrier"
        return kind if key == LEGACY_CARRIER_KEY else f"{kind}: {key}"


def carrier_state(who: _HasCarrierData, key: str) -> int | None:
    """Carrier status of an individual for one source.

    Returns:
        ``OBLIGATE_CARRIER``, ``CARRIER``, 0 when not a carrier, or None
        when the status is unknown (missing value in a ``cr_`` column).
    """
    if key == LEGACY_CARRIER_KEY:
        if who.is_obligate_carrier:
            return OBLIGATE_CARRIER
        return CARRIER if who.is_carrier else 0
    raw = who.extra_phenotypes_raw.get(key)
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return None
    if raw == OBLIGATE_CARRIER:
        return OBLIGATE_CARRIER
    return CARRIER if raw == CARRIER else 0


def resolve_carrier_scheme(
    individuals: Iterable[_HasCarrierData],
    carrier_columns: Iterable[str],
    notation: str,
) -> CarrierScheme:
    """Decide which carrier sources a figure draws, and how.

    Args:
        individuals: Everyone drawn in the figure.
        carrier_columns: ``cr_`` columns (prefix stripped), in file order.
        notation: ``"nsgc2022"`` or ``"legacy"``.

    Returns:
        The figure's carrier scheme.
    """
    people = list(individuals)
    has_legacy_plain = any(p.is_carrier for p in people)
    has_legacy_obligate = any(p.is_obligate_carrier for p in people)
    has_legacy = has_legacy_plain or has_legacy_obligate
    keys: list[str] = []
    if has_legacy and notation != "legacy":
        keys.append(LEGACY_CARRIER_KEY)
    keys.extend(carrier_columns)
    obligate = frozenset(
        k for k in keys if any(carrier_state(p, k) == OBLIGATE_CARRIER for p in people)
    )
    return CarrierScheme(
        notation=notation,
        keys=tuple(keys),
        obligate_keys=obligate,
        has_legacy_carriers=has_legacy,
        has_legacy_plain=has_legacy_plain,
        has_legacy_obligate=has_legacy_obligate,
    )


def carrier_columns(
    phenotype_names: Iterable[str], display_modes: Mapping[str, object]
) -> list[str]:
    """The phenotype columns declared as carrier columns (``cr_`` prefix)."""
    return [n for n in phenotype_names if display_modes.get(n) == PhenotypeDisplayMode.CARRIER]


def check_notation(notation: str) -> None:
    """Raise PedigreeError for an unknown ``carrier_notation`` value."""
    if notation not in CARRIER_NOTATIONS:
        raise PedigreeError(
            f"Unknown carrier_notation {notation!r}: expected one of "
            + ", ".join(repr(n) for n in CARRIER_NOTATIONS)
        )
