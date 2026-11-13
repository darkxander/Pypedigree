"""Configuration and style parameters for pedigree rendering.

All visual parameters are centralized here. No magic numbers
should exist in layout.py or renderer.py.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field, fields, replace
from difflib import get_close_matches
from typing import Any

from pypedigree.errors import PedigreeStyleWarning

DEFAULT_PHENOTYPE_PALETTE = [
    "#e74c3c",  # red
    "#3498db",  # blue
    "#e67e22",  # orange
    "#2ecc71",  # green
    "#9b59b6",  # purple
    "#1abc9c",  # teal
    "#f1c40f",  # yellow
    "#e91e63",  # pink
]


@dataclass
class PedigreeStyle:
    """Visual style configuration for pedigree diagrams.

    All spacing values are in diagram coordinate units.
    """

    # Shape parameters
    shape_size: float = 0.8
    shape_line_width: float = 1.5
    fill_color: str = "black"
    outline_color: str = "black"
    phantom_line_style: str = "--"

    # Line parameters
    line_color: str = "black"
    line_width: float = 1.0
    double_line_gap: float = 0.06
    broken_gap_size: float = 0.15

    # Text parameters
    font_size: float = 8.0
    font_family: str = "sans-serif"
    label_offset: float = 0.15

    # Layout spacing
    generation_gap: float = 2.0
    sibling_gap: float = 1.5
    couple_gap: float = 1.5
    min_node_distance: float = 1.2

    # Carrier notation: "nsgc2022" divides the symbol into patterned sections
    # (Bennett et al. 2022); "legacy" draws the `carrier` column as a
    # half-filled symbol and obligate carriers as a central dot.
    carrier_notation: str = "nsgc2022"
    carrier_dot_radius: float = 0.08

    # Multi-phenotype colors (pheno_name -> hex color)
    phenotype_colors: dict[str, str] = field(default_factory=dict)

    # Show YES/NO text labels inside binary phenotype bands (default: color only)
    show_phenotype_labels: bool = False

    # PGS curve colors (PGS column name -> hex color, e.g. {"PGS_autism": "#3498db"})
    pgs_colors: dict[str, str] = field(default_factory=dict)

    # Custom gradient normalization ranges (pheno_name -> [min, max])
    gradient_range: dict[str, list[float]] = field(default_factory=dict)

    # Diverging (tricolor) gradient config (pheno_name -> {min, mid, max, color_low, color_high})
    gradient_diverging: dict[str, dict[str, Any]] = field(default_factory=dict)

    # Zero-pad individual IDs to this width (0 = no padding)
    zero_pad_ids: int = 0

    # Circular layout parameters
    radial_gap: float = 2.5
    base_radius: float = 3.0

    # Output parameters
    dpi: int = 150
    fig_width: float = 10.0
    fig_height: float = 8.0
    margin: float = 0.5

    @classmethod
    def from_dict(cls, d: dict[str, Any], *, warn_unknown: bool = True) -> PedigreeStyle:
        """Create a style from a dictionary, ignoring unknown keys.

        Unknown keys are skipped rather than raising, so a style file stays
        forward-compatible with newer options. By default each one is
        reported as a `PedigreeStyleWarning` with a suggested spelling,
        since silently ignoring a typo such as `phenotype_color` is
        indistinguishable from the option having no effect.

        Args:
            d: Dictionary of style overrides.
            warn_unknown: Whether to warn about keys that are not style
                fields. Set False to accept a dictionary that legitimately
                carries extra data.

        Returns:
            PedigreeStyle with overrides applied over defaults.
        """
        valid_keys = {f.name for f in fields(cls)}
        unknown = [k for k in d if k not in valid_keys]
        if unknown and warn_unknown:
            for key in unknown:
                close = get_close_matches(str(key), valid_keys, n=1)
                hint = f"; did you mean {close[0]!r}?" if close else ""
                warnings.warn(
                    f"Unknown style option {key!r} ignored{hint}",
                    PedigreeStyleWarning,
                    stacklevel=2,
                )
        filtered = {k: v for k, v in d.items() if k in valid_keys}
        return cls(**filtered)

    def to_dict(self) -> dict[str, Any]:
        """Convert to a plain dictionary."""
        return asdict(self)


def is_pgs_column(name: str) -> bool:
    """Check whether a phenotype column name denotes a polygenic score.

    Args:
        name: Phenotype column name.

    Returns:
        True if the name starts with "PGS" (case-insensitive).
    """
    return name.upper().startswith("PGS")


def resolve_style(
    style: PedigreeStyle | None,
    phenotype_names: Sequence[str] | None = None,
) -> PedigreeStyle:
    """Fill in the pedigree-dependent parts of a style.

    Renderers only draw phenotype bands when ``phenotype_colors`` is
    populated, which depends on the pedigree being rendered rather than on
    the style alone. Every entry point (CLI, PedigreeBuilder,
    render_pedigree) routes through this function to stay consistent.

    The input style is never mutated, and the function is idempotent:
    resolving an already-resolved style returns an equivalent style, so
    repeated renders with the same object are stable.

    Args:
        style: User-supplied style, or None for defaults.
        phenotype_names: Ordered phenotype column names of the pedigree.

    Returns:
        A new PedigreeStyle with auto-assigned phenotype colors.
    """
    resolved = PedigreeStyle() if style is None else replace(style)
    # replace() copies the dataclass but shares its mutable fields.
    resolved.phenotype_colors = dict(resolved.phenotype_colors)
    resolved.pgs_colors = dict(resolved.pgs_colors)
    resolved.gradient_range = {k: list(v) for k, v in resolved.gradient_range.items()}
    resolved.gradient_diverging = {k: dict(v) for k, v in resolved.gradient_diverging.items()}

    # Auto-assign band colors from the palette when the user set none.
    names = list(phenotype_names or [])
    if names and not resolved.phenotype_colors:
        for i, name in enumerate(names):
            if not is_pgs_column(name):
                resolved.phenotype_colors[name] = DEFAULT_PHENOTYPE_PALETTE[
                    i % len(DEFAULT_PHENOTYPE_PALETTE)
                ]

    return resolved


def apply_pgs_spacing(
    style: PedigreeStyle,
    phenotype_names: Sequence[str] | None = None,
) -> PedigreeStyle:
    """Widen node spacing to make room for PGS distribution curves.

    PGS curves are drawn beside each node, so couples and siblings need
    extra horizontal room. Unlike `resolve_style` this is not idempotent —
    each call widens further — so exactly one caller may apply it, which
    is `render_pedigree` just before laying families out.

    Args:
        style: Style to widen.
        phenotype_names: Ordered phenotype column names of the pedigree.

    Returns:
        A new PedigreeStyle, widened if any PGS column is present, else
        `style` unchanged.
    """
    if not any(is_pgs_column(n) for n in (phenotype_names or [])):
        return style

    pgs_extra = style.shape_size * 0.55
    return replace(
        style,
        couple_gap=style.couple_gap + pgs_extra,
        sibling_gap=style.sibling_gap + pgs_extra,
    )
