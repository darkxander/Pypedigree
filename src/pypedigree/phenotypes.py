"""Phenotype display decisions shared by the matplotlib and SVG renderers.

This module decides *what* each phenotype column looks like: its effective
display mode, the color and text of every band, and the rows of the legend.
The renderers only turn these decisions into matplotlib artists or SVG
elements, so a PNG and an interactive SVG of the same pedigree always agree.

Colors are parsed with ``matplotlib.colors`` (and nothing else from
matplotlib), so a style accepts the same color names in both outputs.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Protocol

import matplotlib.colors as mcolors

from pypedigree.carriers import CarrierScheme
from pypedigree.config import PedigreeStyle, is_pgs_column
from pypedigree.errors import PedigreeError
from pypedigree.models import Individual, NodeLayout, Pedigree, Phenotype, PhenotypeDisplayMode

RGB = tuple[float, ...]
Color = str | RGB

# Distinct colors for categorical phenotype values.
CATEGORICAL_COLORS = [
    "#e74c3c",  # red
    "#3498db",  # blue
    "#2ecc71",  # green
    "#f39c12",  # orange
    "#9b59b6",  # purple
    "#1abc9c",  # teal
    "#e67e22",  # dark orange
    "#34495e",  # dark grey
    "#e91e63",  # pink
    "#00bcd4",  # cyan
]

UNKNOWN_FILL = "#d9d9d9"  # binary 0 / missing value
FALLBACK_COLOR = "#999999"  # column without a phenotype_colors entry
SEPARATOR_COLOR = "#cccccc"  # thin line between bands
NA_TEXT_COLOR = "#999999"  # "NA" written in a band
DARK_TEXT = "#333333"  # band text on a light background
DEFAULT_PGS_COLOR = "#e74c3c"
GRADIENT_STEPS = 20  # swatches making up a legend gradient bar

# Shape of the PGS distribution curve: a gaussian of this width on [0, 1].
PGS_CURVE_SIGMA = 0.18


# --- Colors -----------------------------------------------------------------


def to_rgb(color: str) -> tuple[float, float, float]:
    """Convert a color string to (r, g, b) floats in [0, 1].

    Accepts anything matplotlib accepts: hex (``#e74c3c``), a named color
    (``grey``), or a shorthand (``C0``).

    Raises:
        PedigreeError: If the color is not recognised.
    """
    try:
        return mcolors.to_rgb(color)
    except (ValueError, TypeError) as exc:
        raise PedigreeError(
            f"Invalid color {color!r}: expected a hex value like '#e74c3c' "
            f"or a named color like 'grey'"
        ) from exc


def to_hex(color: Color) -> str:
    """Hex form (``#rrggbb``) of a color string or an RGB tuple."""
    if isinstance(color, str):
        return str(mcolors.to_hex(color))
    r, g, b = color[:3]
    return str(mcolors.to_hex((r, g, b)))


def validate_style_colors(style: PedigreeStyle) -> None:
    """Check every color in a style before any drawing starts.

    Reports all offending values at once, so a style file with several
    typos does not have to be fixed one render at a time.

    Raises:
        PedigreeError: If any color is not recognised, listing each one
            with the style option it came from.
    """
    candidates: list[tuple[str, object]] = [
        ("fill_color", style.fill_color),
        ("outline_color", style.outline_color),
        ("line_color", style.line_color),
    ]
    candidates += [(f"phenotype_colors[{k!r}]", v) for k, v in style.phenotype_colors.items()]
    candidates += [(f"pgs_colors[{k!r}]", v) for k, v in style.pgs_colors.items()]
    for name, conf in style.gradient_diverging.items():
        if not isinstance(conf, dict):
            continue
        for key in ("color_low", "color_high"):
            if key in conf:
                candidates.append((f"gradient_diverging[{name!r}][{key!r}]", conf[key]))

    bad: list[str] = []
    for option, value in candidates:
        try:
            # Deliberately unvalidated: to_rgb raising is how we detect a bad
            # value, so the loose type here is the point, not an oversight.
            mcolors.to_rgb(value)  # type: ignore[arg-type]
        except (ValueError, TypeError):
            bad.append(f"{option} = {value!r}")

    if bad:
        raise PedigreeError(
            "Invalid color(s) in style: "
            + "; ".join(bad)
            + ". Expected hex values like '#e74c3c' or named colors like 'grey'."
        )


def blend_with_white(rgb: Sequence[float], alpha: float) -> RGB:
    """Blend a color toward white by a factor (0=white, 1=full color)."""
    return tuple(1.0 - alpha * (1.0 - c) for c in rgb)


def normalize_to_unit(value: float, low: float, high: float) -> float:
    """Map a value onto [0, 1] against an explicit range, clamping outliers.

    Returns:
        The position in [0, 1]; 0.5 for a degenerate range.
    """
    if high == low:
        return 0.5
    return max(0.0, min(1.0, (value - low) / (high - low)))


def diverging_color_at(conf: Mapping[str, Any], t: float) -> RGB:
    """Color of a diverging (tricolor) gradient at a normalized position.

    Shared by the bands and the legend bar so the two can never drift
    apart. The gradient runs ``color_low`` → white at the midpoint →
    ``color_high``.

    Args:
        conf: A ``gradient_diverging`` entry (``min``, ``mid``, ``max``,
            ``color_low``, ``color_high``; each optional).
        t: Position in [0, 1] between ``min`` and ``max``.
    """
    d_min = conf.get("min", 0)
    d_mid = conf.get("mid", 50)
    d_max = conf.get("max", 100)
    c_low = to_rgb(conf.get("color_low", "#e74c3c"))
    c_high = to_rgb(conf.get("color_high", "#2ecc71"))

    mid_t = normalize_to_unit(d_mid, d_min, d_max)
    if t <= mid_t:
        # Below the midpoint: full color_low at t=0, white at the midpoint.
        frac = t / mid_t if mid_t > 0 else 0.0
        return blend_with_white(c_low, 1.0 - frac)
    # Above the midpoint: white at the midpoint, full color_high at t=1.
    frac = (t - mid_t) / (1.0 - mid_t) if mid_t < 1 else 0.0
    return blend_with_white(c_high, frac)


def diverging_color_for_value(conf: Mapping[str, Any], value: float) -> RGB:
    """Color of a diverging gradient for a raw phenotype value."""
    t = normalize_to_unit(value, conf.get("min", 0), conf.get("max", 100))
    return diverging_color_at(conf, t)


def text_color_on(background: Color) -> str:
    """Dark or white text, whichever reads on a band of this color."""
    if background == "white":
        return DARK_TEXT
    r, g, b = to_rgb(background) if isinstance(background, str) else background[:3]
    brightness = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if brightness < 0.55 else DARK_TEXT


def pgs_curve(n_points: int) -> list[tuple[float, float]]:
    """The PGS distribution curve as (position, density) pairs on [0, 1].

    The highest sample has density 1. Both renderers draw the curve beside
    each node and in the legend from these points.
    """
    ts = [i / (n_points - 1) for i in range(n_points)]
    dens = [math.exp(-0.5 * ((t - 0.5) / PGS_CURVE_SIGMA) ** 2) for t in ts]
    peak = max(dens)
    return [(t, d / peak) for t, d in zip(ts, dens, strict=True)]


# --- Per-figure resolution --------------------------------------------------


@dataclass(frozen=True)
class PhenotypeRender:
    """How each phenotype column is drawn, resolved once per figure.

    Bands and legend read the same resolution, so they can never disagree,
    and nothing is carried over between successive renders in one process.

    Attributes:
        effective_modes: Column name → "binary", "gradient", "categorical",
            "outline", "carrier" or "pgs". AUTO columns are resolved from
            their values.
        categorical_colors: Column name → {category value → hex color}, for
            categorical columns only.
    """

    effective_modes: dict[str, str] = field(default_factory=dict)
    categorical_colors: dict[str, dict[str, str]] = field(default_factory=dict)

    def effective(self, name: str) -> str:
        """Effective display mode of a column ("binary" if unknown)."""
        return self.effective_modes.get(name, "binary")

    def categorical_color(self, name: str, value: str) -> str | None:
        """Color assigned to a category value, or None if unassigned."""
        return self.categorical_colors.get(name, {}).get(value)

    def band_names(self, names: Sequence[str]) -> list[str]:
        """Columns drawn as horizontal bands inside the symbol, in order."""
        return [n for n in names if self.effective(n) not in ("pgs", "outline", "carrier")]

    def outline_names(self, names: Sequence[str]) -> list[str]:
        """Columns drawn as a colored symbol border (``o_`` prefix)."""
        return [n for n in names if self.effective(n) == "outline"]

    def pgs_names(self, names: Sequence[str]) -> list[str]:
        """Columns drawn on the PGS distribution curve."""
        return [n for n in names if self.effective(n) == "pgs"]


def is_na_value(value: float | None) -> bool:
    """Whether a phenotype value carries no data (missing or NaN)."""
    return value is None or (isinstance(value, float) and math.isnan(value))


def detect_auto_mode(name: str, individuals: Sequence[Individual]) -> str:
    """Resolve an AUTO phenotype column from all of its values.

    The decision is taken per column rather than per individual, so a
    column never mixes binary and gradient bands across the diagram.

    Returns:
        "categorical" if any value is a string, "gradient" if any numeric
        value falls outside the 0/1/2 coding, "binary" otherwise.
    """
    for ind in individuals:
        if ind.extra_phenotypes_categorical.get(name) is not None:
            return "categorical"
    for ind in individuals:
        raw = ind.extra_phenotypes_raw.get(name)
        if is_na_value(raw) or raw in (0.0, 1.0, 2.0):
            continue
        return "gradient"
    return "binary"


def assign_categorical_colors(
    name: str,
    individuals: Sequence[Individual],
    style: PedigreeStyle,
) -> dict[str, str]:
    """Assign a distinct color to each category value of one column.

    Values are taken in file order, so the mapping is reproducible for a
    given pedigree. The column's own ``phenotype_colors`` entry is used
    for the first category and then excluded from the palette, otherwise
    two categories could be handed the same color (the palettes overlap).

    Returns:
        Mapping of category value → color.
    """
    values: list[str] = []
    for ind in individuals:
        cat_val = ind.extra_phenotypes_categorical.get(name)
        if cat_val is not None and cat_val not in values:
            values.append(cat_val)

    base = style.phenotype_colors.get(name)
    # Skip colors used by other phenotype columns, so a category swatch is
    # not confused with the band of another column on the same individual.
    taken = {c for n, c in style.phenotype_colors.items() if n != name}
    palette = [c for c in CATEGORICAL_COLORS if c != base and c not in taken]
    if len(palette) < len(values):
        # Too many categories to stay clear of the other columns.
        palette = [c for c in CATEGORICAL_COLORS if c != base]
    if base:
        palette.insert(0, base)

    return {val: palette[i % len(palette)] for i, val in enumerate(values)}


def resolve_phenotype_render(
    pedigree: Pedigree,
    phenotype_names: Sequence[str],
    display_modes: Mapping[str, PhenotypeDisplayMode],
    style: PedigreeStyle,
) -> PhenotypeRender:
    """Resolve how every phenotype column of a pedigree is drawn.

    Args:
        pedigree: The pedigree being rendered.
        phenotype_names: Ordered phenotype column names.
        display_modes: Modes declared by column prefixes.
        style: Resolved style.
    """
    individuals = [ind for fam in pedigree.families.values() for ind in fam.individuals.values()]

    effective_modes: dict[str, str] = {}
    categorical_colors: dict[str, dict[str, str]] = {}

    for name in phenotype_names:
        mode = display_modes.get(name, PhenotypeDisplayMode.AUTO)
        if is_pgs_column(name):
            effective_modes[name] = "pgs"
            continue
        if mode == PhenotypeDisplayMode.OUTLINE:
            effective_modes[name] = "outline"
            continue
        if mode == PhenotypeDisplayMode.CARRIER:
            effective_modes[name] = "carrier"
            continue

        if mode == PhenotypeDisplayMode.CATEGORICAL:
            effective = "categorical"
        elif mode in (PhenotypeDisplayMode.GRADIENT, PhenotypeDisplayMode.GRADIENT_INVERTED):
            effective = "gradient"
        elif mode == PhenotypeDisplayMode.BINARY:
            effective = "binary"
        elif name in style.gradient_diverging or name in style.gradient_range:
            # An explicit gradient configuration implies a continuous column.
            effective = "gradient"
        else:
            effective = detect_auto_mode(name, individuals)
        effective_modes[name] = effective

        if effective == "categorical":
            categorical_colors[name] = assign_categorical_colors(name, individuals, style)

    return PhenotypeRender(effective_modes, categorical_colors)


# --- Bands -------------------------------------------------------------------


class _HasPhenotypes(Protocol):
    """A NodeLayout, as seen by the band decisions."""

    extra_phenotypes: dict[str, float | None]
    extra_phenotypes_raw: dict[str, float | None]
    extra_phenotypes_categorical: dict[str, str | None]


@dataclass(frozen=True)
class BandCell:
    """One phenotype band of one individual.

    Attributes:
        color: Fill color: a color string or an RGB tuple.
        text: Text written in the band, or None.
        is_na: Whether the value is "no data" (-9 / NA); the band is white
            and reads "NA" unless every band of the symbol is NA, in which
            case a single centred "NA" replaces them.
    """

    color: Color
    text: str | None
    is_na: bool

    @property
    def text_color(self) -> str:
        """Color of the band's text."""
        return text_color_on(self.color)


def band_cell(
    node: _HasPhenotypes,
    name: str,
    style: PedigreeStyle,
    render: PhenotypeRender,
    display_modes: Mapping[str, PhenotypeDisplayMode],
) -> BandCell:
    """Color and text of one phenotype band of one individual."""
    color_hex = style.phenotype_colors.get(name, FALLBACK_COLOR)
    norm_val = node.extra_phenotypes.get(name)
    raw_val = node.extra_phenotypes_raw.get(name)
    cat_val = node.extra_phenotypes_categorical.get(name)
    mode = display_modes.get(name, PhenotypeDisplayMode.AUTO)
    inverted = mode == PhenotypeDisplayMode.GRADIENT_INVERTED
    # Effective display type is resolved per column, not per individual,
    # so every band of a column is drawn the same way as its legend row.
    effective = render.effective(name)

    if isinstance(norm_val, float) and math.isnan(norm_val):
        return BandCell("white", "NA", True)

    color: Color
    text: str | None = None
    if effective == "categorical":
        cat_color = render.categorical_color(name, cat_val) if cat_val is not None else None
        color = cat_color if cat_color is not None else UNKNOWN_FILL
        text = cat_val
    elif effective == "binary":
        if norm_val is None:
            color = UNKNOWN_FILL
        elif norm_val == 1.0:
            color = color_hex
        else:
            color = "white"
        if style.show_phenotype_labels and norm_val is not None:
            text = "YES" if norm_val == 1.0 else "NO"
    else:
        rgb = to_rgb(color_hex)
        rv = raw_val if raw_val is not None else 0
        if norm_val is None and raw_val is None:
            color = UNKNOWN_FILL
        elif name in style.gradient_diverging:
            color = diverging_color_for_value(style.gradient_diverging[name], rv)
        elif name in style.gradient_range:
            g_min, g_max = style.gradient_range[name]
            t = normalize_to_unit(rv, g_min, g_max)
            color = blend_with_white(rgb, 1.0 - t if inverted else t)
        else:
            nv = norm_val if norm_val is not None else 0
            color = blend_with_white(rgb, 1.0 - nv if inverted else nv)
        if raw_val is not None:
            # Gradient values are always written.
            text = f"{name}: {raw_val:g}"
    return BandCell(color, text, False)


def pregnancy_loss_fill(
    node: NodeLayout,
    band_names: Sequence[str],
    style: PedigreeStyle,
    render: PhenotypeRender,
    display_modes: Mapping[str, PhenotypeDisplayMode],
) -> str:
    """Fill of a pregnancy-loss triangle, which is too small for bands.

    In a figure with phenotype bands, the color the individual's first band
    would have (white for no data); otherwise the affection status, as for
    any symbol.
    """
    if band_names:
        cell = band_cell(node, band_names[0], style, render, display_modes)
        return "white" if cell.is_na else to_hex(cell.color)
    return style.fill_color if node.phenotype == Phenotype.AFFECTED else "white"


def all_bands_na(node: _HasPhenotypes, band_names: Sequence[str]) -> bool:
    """Whether every phenotype band of a symbol is "no data"."""
    return all(is_na_value(node.extra_phenotypes.get(n)) for n in band_names)


def outline_color(
    node: _HasPhenotypes, outline_names: Sequence[str], style: PedigreeStyle
) -> str | None:
    """Border color from the first ``o_`` column this individual is affected in."""
    for name in outline_names:
        if node.extra_phenotypes.get(name) == 1.0:
            return style.phenotype_colors.get(name, FALLBACK_COLOR)
    return None


def active_pgs(node: _HasPhenotypes, pgs_names: Sequence[str]) -> list[tuple[str, float]]:
    """PGS columns with a value for this individual, with their [0, 1] position."""
    out = []
    for name in pgs_names:
        value = node.extra_phenotypes.get(name)
        if not is_na_value(value):
            out.append((name, value or 0.0))
    return out


# --- Legend ------------------------------------------------------------------

# Legend panel geometry, in the legend's own coordinates (y pointing up). One unit is a
# fixed multiple of the base font size (`style.font_size`), horizontally and
# vertically, so the panel has the same physical size and proportions whatever
# the size of the pedigree above it, and its swatches grow with its text.
# These are internal layout constants, not user-tunable style.
LEG_UNIT_X_EM = 4.2  # one horizontal unit, in font sizes
LEG_UNIT_Y_EM = 2.0  # one vertical unit, in font sizes
LEG_X0 = 0.4  # left margin for swatches and icons
LEG_GRAD_W = 2.6  # width of a gradient bar
LEG_HEADER_H = 1.8  # title + male/female symbol row
LEG_BINARY_HEADER_H = 0.5  # Affected / Unaffected / Unknown captions
LEG_ROW_H = 1.2  # one phenotype or carrier row
LEG_PAD = 0.8  # room below the last row for its captions
LEG_REFERENCE_FONT = 8.0  # font size the marker sizes below are tuned for
LEG_SWATCH_EDGE = "#cccccc"
LEG_CAPTION_COLOR = "#999999"
LEG_CATEGORY_COLOR = "#666666"


@dataclass(frozen=True)
class LegendMetrics:
    """Derived geometry of the legend panel.

    Lengths are in legend units (see `LEG_UNIT_X_EM`); font and marker sizes
    are in points. All of them scale with `style.font_size`, nothing with
    the figure.

    Attributes:
        n_rows: Number of phenotype/carrier rows (at least 1).
        row_h: Height of one row.
        total_h: Height of the whole panel.
        scale: Marker size relative to the reference font size.
        sw_h: Swatch height.
        sw_w: Swatch width.
        sw_gap: Horizontal gap between swatches.
        text_x: Label x for rows built from three swatches.
        icon_text_x: Label x for rows built from a single icon.
        fs: Row font size.
        label_fs: Caption font size (column headers, scale ends).
        marker_size: Size of plot-marker icons.
        unit_x: One horizontal legend unit, in inches.
        unit_y: One vertical legend unit, in inches.
    """

    n_rows: int
    row_h: float
    total_h: float
    scale: float
    sw_h: float
    sw_w: float
    sw_gap: float
    text_x: float
    icon_text_x: float
    fs: float
    label_fs: float
    marker_size: float
    unit_x: float
    unit_y: float

    def row_center(self, index: int, y_top: float) -> float:
        """Vertical centre of row `index`, counting down from `y_top`."""
        return y_top - index * self.row_h - self.row_h / 2

    def swatch_x(self, index: int) -> float:
        """Left edge of the `index`-th swatch in a multi-swatch row."""
        return LEG_X0 + index * (self.sw_w + self.sw_gap)


def legend_metrics(
    n_pheno_rows: int, style: PedigreeStyle, *, binary_headers: bool = False
) -> LegendMetrics:
    """Compute the legend panel geometry for a given number of rows.

    Rows keep a constant height: the panel grows with its row count, and
    the figure grows with the panel (see `_add_legend_axes`).

    Args:
        n_pheno_rows: Number of phenotype and carrier rows to fit.
        style: Style supplying the base font size and family.
        binary_headers: Whether the Affected / Unaffected / Unknown captions
            are drawn above the rows.

    Returns:
        The derived metrics.
    """
    n_rows = max(n_pheno_rows, 1)
    row_h = LEG_ROW_H
    sw_h = row_h * 0.45
    sw_w = 0.8
    sw_gap = 0.35

    # Three-swatch rows put their label just past the last swatch.
    text_x = LEG_X0 + 3 * sw_w + 2 * sw_gap + 0.3
    # Single-icon rows align on the widest icon: PGS curve (ends at 1.3),
    # outline swatch, or carrier marker (~1.2).
    icon_right = max(1.3, LEG_X0 + sw_w, 1.2)
    icon_text_x = icon_right + 0.3

    fs = style.font_size * 0.85
    em = style.font_size / 72  # inches
    header_h = LEG_HEADER_H + (LEG_BINARY_HEADER_H if binary_headers else 0.0)

    return LegendMetrics(
        n_rows=n_rows,
        row_h=row_h,
        total_h=header_h + n_rows * row_h + LEG_PAD,
        scale=style.font_size / LEG_REFERENCE_FONT,
        sw_h=sw_h,
        sw_w=sw_w,
        sw_gap=sw_gap,
        text_x=text_x,
        icon_text_x=icon_text_x,
        fs=fs,
        label_fs=fs * 0.75,
        marker_size=style.font_size,
        unit_x=LEG_UNIT_X_EM * em,
        unit_y=LEG_UNIT_Y_EM * em,
    )


@dataclass(frozen=True)
class LegendRow:
    """One row of the legend panel.

    Attributes:
        kind: "binary", "gradient", "categorical", "outline", "pgs",
            "legacy_carrier" (half-fill / dot) or "carrier" (patterned).
        label: Row text.
        color: Swatch color ("binary", "outline").
        colors: Gradient bar, left to right ("gradient").
        low: Caption under the left end of a gradient bar.
        high: Caption under the right end of a gradient bar.
        categories: Category value → color ("categorical").
        pgs: PGS column → marker color ("pgs").
        pattern: Carrier pattern name ("carrier").
        obligate: Obligate-carrier row ("carrier", "legacy_carrier").
    """

    kind: str
    label: str
    color: str = ""
    colors: tuple[RGB, ...] = ()
    low: str = ""
    high: str = ""
    categories: Mapping[str, str] = field(default_factory=dict)
    pgs: Mapping[str, str] = field(default_factory=dict)
    pattern: str = ""
    obligate: bool = False


def pgs_legend_description(pgs_names: Sequence[str]) -> str:
    """Caption of the PGS row, from its column names."""
    if len(pgs_names) == 1 and pgs_names[0].upper() == "PGS":
        return "PGS distribution (vertical) with individual score"
    traits = [n[4:] if n.upper().startswith("PGS_") else n for n in pgs_names]
    return f"PGS distribution ({', '.join(traits)})"


def legend_rows(
    phenotype_names: Sequence[str],
    style: PedigreeStyle,
    render: PhenotypeRender,
    display_modes: Mapping[str, PhenotypeDisplayMode],
    carriers: CarrierScheme,
) -> list[LegendRow]:
    """The rows of the legend, in drawing order.

    One row per band column, then one per outline column, one shared row
    for all PGS columns, and the carrier rows: the half-fill and dot rows
    with the legacy notation, then one per carrier source and per
    obligate-carrier source.
    """
    rows: list[LegendRow] = []
    names = list(phenotype_names)
    steps = [g / (GRADIENT_STEPS - 1) for g in range(GRADIENT_STEPS)]

    for name in render.band_names(names):
        color = style.phenotype_colors.get(name, FALLBACK_COLOR)
        effective = render.effective(name)
        if name in style.gradient_diverging:
            conf = style.gradient_diverging[name]
            rows.append(
                LegendRow(
                    "gradient",
                    name,
                    colors=tuple(diverging_color_at(conf, t) for t in steps),
                    low=str(conf.get("min", 0)),
                    high=str(conf.get("max", 100)),
                )
            )
        elif effective == "gradient":
            inverted = display_modes.get(name) == PhenotypeDisplayMode.GRADIENT_INVERTED
            rgb = to_rgb(color)
            rows.append(
                LegendRow(
                    "gradient",
                    name,
                    colors=tuple(blend_with_white(rgb, 1.0 - t if inverted else t) for t in steps),
                    low="Low",
                    high="High",
                )
            )
        elif effective == "categorical":
            rows.append(
                LegendRow("categorical", name, categories=render.categorical_colors.get(name, {}))
            )
        else:
            rows.append(LegendRow("binary", name, color=color))

    for name in render.outline_names(names):
        rows.append(
            LegendRow(
                "outline",
                f"Outline: {name} (colored border when affected)",
                color=style.phenotype_colors.get(name, FALLBACK_COLOR),
            )
        )

    pgs_names = render.pgs_names(names)
    if pgs_names:
        rows.append(
            LegendRow(
                "pgs",
                pgs_legend_description(pgs_names),
                pgs={n: style.pgs_colors.get(n, DEFAULT_PGS_COLOR) for n in pgs_names},
            )
        )

    if carriers.legacy:
        if carriers.has_legacy_plain:
            rows.append(LegendRow("legacy_carrier", "Carrier (half-filled = heterozygous carrier)"))
        if carriers.has_legacy_obligate:
            rows.append(
                LegendRow(
                    "legacy_carrier",
                    "Obligate carrier (center dot = inferred from pedigree)",
                    obligate=True,
                )
            )

    for key in carriers.keys:
        for obligate in (False, True):
            if obligate and key not in carriers.obligate_keys:
                continue
            rows.append(
                LegendRow(
                    "carrier",
                    carriers.label(key, obligate),
                    pattern=carriers.pattern(key),
                    obligate=obligate,
                )
            )
    return rows
