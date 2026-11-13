"""Matplotlib-based pedigree diagram renderer.

Translates FamilyLayout geometry into matplotlib drawing calls.
Never computes positions — only draws from pre-computed layout data.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
from matplotlib.transforms import Bbox

from pypedigree import symbols
from pypedigree.carriers import (
    CARRIER,
    OBLIGATE_CARRIER,
    CarrierScheme,
    carrier_columns,
    carrier_state,
    check_notation,
    resolve_carrier_scheme,
)
from pypedigree.config import PedigreeStyle, apply_pgs_spacing, resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.layout import compute_layout, compute_layout_circular
from pypedigree.models import (
    DuplicateArc,
    EdgeLayout,
    Family,
    FamilyLayout,
    LineStyle,
    NodeLayout,
    Pedigree,
    Phenotype,
    PhenotypeDisplayMode,
    Pregnancy,
    Sex,
    TwinConnector,
    TwinType,
)
from pypedigree.phenotypes import (
    DEFAULT_PGS_COLOR,
    FALLBACK_COLOR,
    LEG_BINARY_HEADER_H,
    LEG_CAPTION_COLOR,
    LEG_CATEGORY_COLOR,
    LEG_GRAD_W,
    LEG_PAD,
    LEG_SWATCH_EDGE,
    LEG_X0,
    PGS_CURVE_SIGMA,
    SEPARATOR_COLOR,
    UNKNOWN_FILL,
    LegendMetrics,
    LegendRow,
    PhenotypeRender,
    active_pgs,
    all_bands_na,
    band_cell,
    legend_metrics,
    legend_rows,
    outline_color,
    pregnancy_loss_fill,
    resolve_phenotype_render,
    text_color_on,
    validate_style_colors,
)


def render_pedigree(
    pedigree: Pedigree,
    style: PedigreeStyle | None = None,
    output: str | Path | None = None,
    ax: Axes | None = None,
    show_labels: bool = True,
    show_legend: bool = True,
    circular: bool = False,
) -> Figure:
    """Render all families in a pedigree.

    Args:
        pedigree: The parsed pedigree data.
        style: Visual style parameters. Uses defaults if None.
        output: If provided, save the figure to this path.
        ax: If provided, render onto this existing Axes.
        show_labels: Whether to draw sample_id labels under shapes.
        show_legend: Whether to draw the phenotype legend (default: True).
        circular: If True, use circular (radial) layout instead of
            the default hierarchical top-down layout.

    Returns:
        The matplotlib Figure.
    """
    if style is None:
        style = PedigreeStyle()

    validate_style_colors(style)

    if circular:
        style = replace(style, shape_size=style.shape_size * 1.5)

    layout_fn = compute_layout_circular if circular else compute_layout

    families = list(pedigree.families.values())

    if not families:
        fig, ax_obj = plt.subplots(figsize=(style.fig_width, style.fig_height))
        ax_obj.set_xlim(-1, 1)
        ax_obj.set_ylim(-1, 1)
        ax_obj.axis("off")
        return fig

    pheno_names = pedigree.phenotype_names
    display_modes = pedigree.phenotype_display_modes

    # Fill in palette colors (idempotent, so harmless if the caller already
    # resolved), then widen for PGS curves — this is the only place that
    # applies the PGS spacing, so it can never be applied twice.
    style = apply_pgs_spacing(resolve_style(style, pheno_names), pheno_names)
    pheno_render = resolve_phenotype_render(pedigree, pheno_names or [], display_modes, style)
    check_notation(style.carrier_notation)
    carriers = resolve_carrier_scheme(
        (ind for fam in families for ind in fam.individuals.values()),
        carrier_columns(pheno_names or [], display_modes),
        style.carrier_notation,
    )

    if ax is not None:
        fig = _root_figure(ax)
        for fam in families:
            layout = layout_fn(fam, style)
            _render_layout(
                ax,
                layout,
                style,
                show_labels,
                pheno_names,
                display_modes,
                pheno_render,
                carriers=carriers,
            )
        _finalize_axes(ax)
        return fig

    if len(families) == 1:
        layout = layout_fn(families[0], style)
        fig_w, fig_h = _auto_figsize(layout, style)
        if circular:
            # Square figure for circular layout
            size = max(fig_w, fig_h)
            fig_w, fig_h = size, size
        fig, ax_obj = plt.subplots(figsize=(fig_w, fig_h))
        _render_layout(
            ax_obj,
            layout,
            style,
            show_labels,
            pheno_names,
            display_modes,
            pheno_render,
            carriers=carriers,
        )
        _finalize_axes(ax_obj)
        if circular:
            ax_obj.set_aspect("equal")
    else:
        n_fams = len(families)
        fig, axes = plt.subplots(
            1,
            n_fams,
            figsize=(style.fig_width * n_fams, style.fig_height),
        )
        if n_fams == 1:
            axes = [axes]
        for fam, ax_obj in zip(families, axes, strict=True):
            layout = layout_fn(fam, style)
            _render_layout(
                ax_obj,
                layout,
                style,
                show_labels,
                pheno_names,
                display_modes,
                pheno_render,
                carriers=carriers,
            )
            ax_obj.set_title(f"Family: {fam.family_id}", fontsize=style.font_size + 2)
            _finalize_axes(ax_obj)
            if circular:
                ax_obj.set_aspect("equal")

    fig.tight_layout()

    # Draw the legend after tight_layout: it grows the figure downward by its
    # own height instead of taking a share of the pedigree's.
    show_pheno_legend = show_legend and pheno_names and style.phenotype_colors
    show_carrier_legend = show_legend and (carriers.has_legacy_carriers or bool(carriers.keys))
    if show_pheno_legend or show_carrier_legend:
        _draw_legend(
            fig,
            pheno_names or [],
            style,
            carriers=carriers,
            display_modes=display_modes,
            pheno_render=pheno_render,
        )

    if output is not None:
        fig.savefig(str(output), dpi=style.dpi, bbox_inches="tight")

    return fig


def render_family(
    family: Family,
    style: PedigreeStyle | None = None,
    ax: Axes | None = None,
    show_labels: bool = True,
) -> tuple[Figure, FamilyLayout]:
    """Render a single family.

    Args:
        family: The family to render.
        style: Visual style parameters.
        ax: Optional existing Axes.
        show_labels: Whether to draw labels.

    Returns:
        Tuple of (Figure, FamilyLayout).
    """
    if style is None:
        style = PedigreeStyle()

    layout = compute_layout(family, style)

    if ax is None:
        fig_w, fig_h = _auto_figsize(layout, style)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    else:
        fig = _root_figure(ax)

    check_notation(style.carrier_notation)
    carriers = resolve_carrier_scheme(family.individuals.values(), [], style.carrier_notation)
    # No phenotype names/modes: bands need the pedigree-level resolution.
    _render_layout(ax, layout, style, show_labels, carriers=carriers)
    _finalize_axes(ax)

    return fig, layout


def _root_figure(ax: Axes) -> Figure:
    """Return the root Figure an Axes belongs to.

    `Axes.get_figure()` yields a SubFigure when the axes lives inside one.
    A SubFigure's own `figure` attribute walks up to the root, which is what
    `render_pedigree` promises to return. (`get_figure(root=True)` would be
    more direct but only exists from matplotlib 3.10.)

    Args:
        ax: The axes to look up.

    Returns:
        The root Figure.
    """
    figure = ax.get_figure()
    if figure is None:
        raise PedigreeError("The supplied Axes is not attached to a figure")
    if isinstance(figure, SubFigure):
        return figure.figure
    return figure


def _auto_figsize(
    layout: FamilyLayout,
    style: PedigreeStyle,
) -> tuple[float, float]:
    """Compute figure dimensions from layout bounds.

    For small pedigrees, uses the default fig_width/fig_height.
    For large pedigrees, scales the figure so nodes have enough room.

    Returns:
        (width, height) in inches.
    """
    if not layout.nodes:
        return (style.fig_width, style.fig_height)

    xs = [n.x for n in layout.nodes]
    ys = [n.y for n in layout.nodes]
    x_range = max(xs) - min(xs)
    y_range = max(ys) - min(ys)

    # Scale: inches per layout unit, with a minimum of the default size.
    # Horizontal needs more room for labels; vertical for generation spacing.
    h_scale = 0.8
    v_scale = 1.2
    w = max(style.fig_width, x_range * h_scale + 2 * style.margin)
    h = max(style.fig_height, y_range * v_scale + 2 * style.margin)

    # Cap at reasonable maximum to avoid huge files
    w = min(w, 80.0)
    h = min(h, 60.0)

    return (w, h)


def _render_layout(
    ax: Axes,
    layout: FamilyLayout,
    style: PedigreeStyle,
    show_labels: bool,
    phenotype_names: list[str] | None = None,
    display_modes: dict[str, PhenotypeDisplayMode] | None = None,
    pheno_render: PhenotypeRender | None = None,
    *,
    carriers: CarrierScheme | None = None,
) -> None:
    """Render a FamilyLayout onto an Axes.

    Args:
        ax: The matplotlib Axes to draw on.
        layout: Pre-computed layout geometry.
        style: Visual style parameters.
        show_labels: Whether to draw labels.
        phenotype_names: Ordered list of extra phenotype names (for banded rendering).
        display_modes: Map of phenotype name → display mode (from column prefix).
        pheno_render: Per-figure phenotype display resolution. Defaults to an
            empty resolution, which disables banded rendering.
        carriers: Per-figure carrier notation. Defaults to no carrier source
            and the NSGC 2022 notation.
    """
    if display_modes is None:
        display_modes = {}
    if pheno_render is None:
        pheno_render = PhenotypeRender()
    if carriers is None:
        carriers = CarrierScheme()

    # Draw edges first (behind nodes)
    for edge in layout.edges:
        _draw_edge(ax, edge, style)

    # Draw twin connectors
    for tc in layout.twin_connectors:
        _draw_twin_connector(ax, tc, style)

    # Draw duplicate arcs (dotted curves between original and duplicate)
    for arc in layout.duplicate_arcs:
        _draw_duplicate_arc(ax, arc, style)

    # Determine if we use banded rendering
    use_banded = bool(phenotype_names and style.phenotype_colors and pheno_render.effective_modes)
    names = phenotype_names if use_banded and phenotype_names is not None else []
    # Carrier columns are neither: they are drawn from `carriers`.
    band_names = pheno_render.band_names(names)
    outline_names = pheno_render.outline_names(names)
    pgs_names = pheno_render.pgs_names(names)

    # Draw nodes on top
    for node in layout.nodes:
        # Count elements before drawing to apply rotation afterwards
        n_patches_before = len(ax.patches)
        n_lines_before = len(ax.lines)
        n_texts_before = len(ax.texts)

        banded = use_banded and bool(node.extra_phenotypes) and not node.is_phantom
        if symbols.is_pregnancy_loss(node.notation):
            fill = pregnancy_loss_fill(
                node, band_names if banded else [], style, pheno_render, display_modes
            )
            _draw_pregnancy_loss(ax, node, style, fill)
        elif banded:
            _draw_banded_node(
                ax,
                node,
                style,
                band_names,
                display_modes,
                outline_color(node, outline_names, style),
                pheno_render,
                carriers=carriers,
            )
            if pgs_names:
                _draw_pgs_ridge(ax, node, style, pgs_names)
        else:
            _draw_node(ax, node, style, carriers)
        _draw_notation(ax, node, style, carriers, banded)

        # Apply rotation to all patches/lines added for this node
        if node.rotation != 0.0:
            _rotate_node_elements(
                ax,
                node.x,
                node.y,
                node.rotation,
                n_patches_before,
                n_lines_before,
                n_texts_before,
            )

            # Drawn after rotation so it stays upright and points at the
            # rotated outline rather than at the unrotated corner.
            if node.is_proband or node.notation.is_consultand:
                half = style.shape_size / 2
                reach = _shape_boundary_radius(node.sex, half, node.rotation, _PROBAND_DIRECTION)
                _draw_proband_marker(
                    ax,
                    node.x + reach * math.cos(_PROBAND_DIRECTION),
                    node.y + reach * math.sin(_PROBAND_DIRECTION),
                    half,
                    style,
                    letter=node.is_proband,
                )

        if show_labels and node.label:
            _draw_label(ax, node, style)


# The proband arrow points at the lower-left of the shape (Bennett et al. 2008).
_PROBAND_DIRECTION = math.radians(225)


def _shape_boundary_radius(sex: Sex, half: float, rotation: float, direction: float) -> float:
    """Distance from a shape's centre to its outline in a given direction.

    Args:
        sex: Decides the shape: square (male), circle (female) or diamond.
        half: Half the shape size.
        rotation: Rotation applied to the shape, in radians.
        direction: Screen direction to measure along, in radians.

    Returns:
        The distance to the outline along `direction`.
    """
    return symbols.shape_reach(sex, half, direction, rotation)


def _shape_lower_extent(sex: Sex, half: float, rotation: float) -> float:
    """How far a rotated shape reaches below its centre.

    Args:
        sex: Decides the shape: square (male), circle (female) or diamond.
        half: Half the shape size.
        rotation: Rotation applied to the shape, in radians.

    Returns:
        The vertical distance from the centre to the shape's lowest point.
    """
    c, s = abs(math.cos(rotation)), abs(math.sin(rotation))
    if sex == Sex.MALE:
        # A square's lowest point is a corner once rotated.
        return half * (c + s)
    if sex == Sex.FEMALE:
        return half
    return half * _DIAMOND_SCALE * max(c, s)


def _rotate_node_elements(
    ax: Axes,
    cx: float,
    cy: float,
    angle: float,
    n_patches_before: int,
    n_lines_before: int,
    n_texts_before: int,
) -> None:
    """Rotate all patches, lines, and texts added after the given counts around (cx, cy)."""
    import math

    from matplotlib.transforms import Affine2D

    rot = Affine2D().rotate_around(cx, cy, angle)
    rot_transform = rot + ax.transData

    # First pass: rotate all patches and track clip source patches
    new_patches = ax.patches[n_patches_before:]
    for patch in new_patches:
        patch.set_transform(rot_transform)

    # Second pass: re-apply clip paths so they pick up the rotated transform
    for patch in new_patches:
        cp = patch.get_clip_path()
        if cp is not None and hasattr(cp, "_patch"):
            # Re-set the clip path to force matplotlib to use updated transform
            patch.set_clip_path(cp._patch)

    # Rotate lines (deceased marker, proband arrow, etc.)
    for line in ax.lines[n_lines_before:]:
        xs = np.asarray(line.get_xdata(), dtype=float)
        ys = np.asarray(line.get_ydata(), dtype=float)
        if xs.size > 0:
            pts = rot.transform(np.column_stack((xs, ys)))
            line.set_data(pts[:, 0], pts[:, 1])

    # Rotate text positions (band labels like "IQ: 105"). Text follows the
    # shape but is never drawn upside down: past a quarter turn it is
    # flipped by 180°, which keeps it centred in the same band.
    text_angle = (math.degrees(angle) + 180.0) % 360.0 - 180.0
    if text_angle > 90.0:
        text_angle -= 180.0
    elif text_angle < -90.0:
        text_angle += 180.0
    for txt in ax.texts[n_texts_before:]:
        tx, ty = txt.get_position()
        new_x, new_y = rot.transform_point((tx, ty))
        txt.set_position((float(new_x), float(new_y)))
        txt.set_rotation(text_angle)


def _draw_node(
    ax: Axes, node: NodeLayout, style: PedigreeStyle, carriers: CarrierScheme | None = None
) -> None:
    """Draw an individual's shape at the node position.

    With the NSGC 2022 carrier notation (the default), an individual with
    something to show besides affection status gets a symbol divided into
    vertical sections, left to right: affected (filled) first, then one
    patterned section per carrier result (see ``pypedigree.carriers``).

    Args:
        ax: The Axes to draw on.
        node: The node to draw.
        style: Visual style.
        carriers: Per-figure carrier notation.
    """
    if carriers is None:
        carriers = CarrierScheme()
    legacy_carrier = carriers.legacy and (node.is_carrier or node.is_obligate_carrier)
    affected = node.phenotype == Phenotype.AFFECTED and not legacy_carrier

    edgecolor = style.outline_color
    linewidth = style.shape_line_width
    linestyle = style.phantom_line_style if node.is_phantom else "-"
    half = style.shape_size / 2

    marks = carriers.marks(node)
    sections: list[_Section] = []
    if marks:
        if affected:
            sections.append(_Section(fill=style.fill_color))
        sections.extend(_Section(hatch=_hatch(mk.pattern, mk.obligate)) for mk in marks)
        facecolor = "white"
    else:
        facecolor = style.fill_color if affected else "white"

    _draw_shape(ax, node.sex, node.x, node.y, half, facecolor, edgecolor, linewidth, linestyle)
    if sections:
        _draw_sections(ax, node, half, sections, style)
        outline = _draw_shape(
            ax, node.sex, node.x, node.y, half, "none", edgecolor, linewidth, linestyle
        )
        outline.set_zorder(6)

    if carriers.legacy:
        if node.is_carrier:
            _draw_carrier_half_fill(ax, node.x, node.y, half, node.sex, style)
        if node.is_obligate_carrier:
            _draw_obligate_carrier_dot(ax, node.x, node.y, style)

    # Deceased marker
    if node.is_deceased:
        _draw_deceased_marker(ax, node.x, node.y, half, style)

    # Proband marker. A rotated node gets its marker after rotation instead,
    # so the "P" stays upright (see _render_layout).
    if (node.is_proband or node.notation.is_consultand) and node.rotation == 0.0:
        _draw_proband_marker(ax, node.x - half, node.y - half, half, style, letter=node.is_proband)


@dataclass(frozen=True)
class _Section:
    """One vertical section of a divided symbol: a fill color or a hatch."""

    fill: str = "white"
    hatch: str | None = None


# Matplotlib hatch for each carrier pattern. Obligate carriers get the same
# pattern at half the density, so a reader sees which source they belong to.
_MPL_HATCHES = {
    "diagonal": "////",
    "antidiagonal": "\\\\\\\\",
    "vertical": "||||",
    "horizontal": "----",
    "crossed": "xxxx",
    "grid": "++++",
    "dotted": "....",
    "circles": "oooo",
}


def _hatch(pattern: str, obligate: bool) -> str:
    """Matplotlib hatch string for a carrier pattern."""
    hatch = _MPL_HATCHES[pattern]
    return hatch[: len(hatch) // 2] if obligate else hatch


def _shape_half_width(sex: Sex, half: float) -> float:
    """Horizontal half-extent of a shape (the diamond reaches further)."""
    return symbols.shape_half_width(sex, half)


def _draw_shape(
    ax: Axes,
    sex: Sex,
    cx: float,
    cy: float,
    half: float,
    facecolor: str,
    edgecolor: str,
    linewidth: float,
    linestyle: str,
) -> mpatches.Patch:
    """Draw the square, circle or diamond of an individual; return the patch."""
    if sex == Sex.MALE:
        return _draw_square(ax, cx, cy, half, facecolor, edgecolor, linewidth, linestyle)
    if sex == Sex.FEMALE:
        return _draw_circle(ax, cx, cy, half, facecolor, edgecolor, linewidth, linestyle)
    return _draw_diamond(ax, cx, cy, half, facecolor, edgecolor, linewidth, linestyle)


def _draw_sections(
    ax: Axes, node: NodeLayout, half: float, sections: list[_Section], style: PedigreeStyle
) -> None:
    """Fill a symbol with equal vertical sections, left to right."""
    clip_patch = _get_shape_clip_path(node.x, node.y, half, node.sex)
    ax.add_patch(clip_patch)
    reach = _shape_half_width(node.sex, half)
    width = 2 * reach / len(sections)
    for i, section in enumerate(sections):
        left = node.x - reach + i * width
        rect = mpatches.Rectangle(
            (left, node.y - reach * 1.5),
            width,
            reach * 3,
            facecolor=section.fill,
            edgecolor=style.fill_color if section.hatch else "none",
            linewidth=0,
            hatch=section.hatch,
            zorder=4,
        )
        rect.set_clip_path(clip_patch)
        ax.add_patch(rect)
        if i > 0:
            ax.plot(
                [left, left],
                [node.y - reach * 1.5, node.y + reach * 1.5],
                color=style.outline_color,
                linewidth=style.shape_line_width * 0.5,
                zorder=5,
                clip_path=clip_patch,
                clip_on=True,
            )


def _draw_square(
    ax: Axes,
    cx: float,
    cy: float,
    half: float,
    facecolor: str,
    edgecolor: str,
    linewidth: float,
    linestyle: str,
) -> mpatches.Patch:
    """Draw a square (male) centered at (cx, cy)."""
    rect = mpatches.FancyBboxPatch(
        (cx - half, cy - half),
        half * 2,
        half * 2,
        boxstyle="square,pad=0",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        linestyle=linestyle,
    )
    ax.add_patch(rect)
    return rect


def _draw_circle(
    ax: Axes,
    cx: float,
    cy: float,
    half: float,
    facecolor: str,
    edgecolor: str,
    linewidth: float,
    linestyle: str,
) -> mpatches.Patch:
    """Draw a circle (female) centered at (cx, cy)."""
    circle = mpatches.Circle(
        (cx, cy),
        radius=half,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        linestyle=linestyle,
    )
    ax.add_patch(circle)
    return circle


# The unknown-sex diamond's vertices sit this far from its centre, relative
# to half the shape size, so it reads as the same size as a square.
_DIAMOND_SCALE = symbols.DIAMOND_SCALE


def _draw_diamond(
    ax: Axes,
    cx: float,
    cy: float,
    half: float,
    facecolor: str,
    edgecolor: str,
    linewidth: float,
    linestyle: str,
) -> mpatches.Patch:
    """Draw a diamond (unknown sex) centered at (cx, cy)."""
    diamond = mpatches.RegularPolygon(
        (cx, cy),
        numVertices=4,
        radius=half * _DIAMOND_SCALE,
        orientation=0,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        linestyle=linestyle,
    )
    ax.add_patch(diamond)
    return diamond


def _draw_deceased_marker(
    ax: Axes,
    cx: float,
    cy: float,
    half: float,
    style: PedigreeStyle,
) -> None:
    """Draw a diagonal line from bottom-left to top-right, extending slightly past shape."""
    overshoot = half * 0.2
    ax.plot(
        [cx - half - overshoot, cx + half + overshoot],
        [cy - half - overshoot, cy + half + overshoot],
        color=style.outline_color,
        linewidth=style.line_width,
        zorder=10,
    )


def _draw_proband_marker(
    ax: Axes,
    tip_x: float,
    tip_y: float,
    half: float,
    style: PedigreeStyle,
    *,
    letter: bool = True,
) -> None:
    """Draw the 'P' label and a short arrow pointing up-right at (tip_x, tip_y).

    The tip is the bottom-left corner of an unrotated shape, or the point of
    a rotated shape's outline in that direction. The marker itself is always
    drawn upright.

    Args:
        ax: The Axes to draw on.
        tip_x: X coordinate the arrow points at.
        tip_y: Y coordinate the arrow points at.
        half: Half the shape size, which scales the marker.
        style: Visual style.
        letter: Write the "P" of the proband; the consultand's arrow has none.
    """
    arrow_start_x = tip_x - half * 0.5
    arrow_start_y = tip_y - half * 0.5

    # Draw the "P"
    if letter:
        ax.text(
            arrow_start_x - half * 0.3,
            arrow_start_y - half * 0.1,
            "P",
            fontsize=style.font_size,
            fontfamily=style.font_family,
            fontweight="bold",
            ha="center",
            va="center",
            color=style.line_color,
        )

    # Small diagonal arrow pointing up-right toward the shape
    ax.annotate(
        "",
        xy=(tip_x, tip_y),
        xytext=(arrow_start_x, arrow_start_y),
        arrowprops=dict(
            arrowstyle="->",
            color=style.line_color,
            lw=style.line_width,
        ),
    )


def _draw_pregnancy_loss(ax: Axes, node: NodeLayout, style: PedigreeStyle, fill: str) -> None:
    """A triangle for a miscarriage, an ectopic pregnancy or a termination.

    Filled as decided by ``phenotypes.pregnancy_loss_fill``; a termination
    is struck through. Bands and carrier sections are not drawn on a
    triangle.
    """
    half = style.shape_size / 2
    points = symbols.triangle_points(node.x, node.y, half)
    ax.add_patch(
        mpatches.Polygon(
            points,
            closed=True,
            facecolor=fill,
            edgecolor=style.outline_color,
            linewidth=style.shape_line_width,
            joinstyle="miter",
        )
    )
    if node.notation.pregnancy == Pregnancy.TERMINATION:
        (x1, y1), (x2, y2) = symbols.termination_stroke(node.x, node.y, half)
        ax.plot(
            [x1, x2], [y1, y2], color=style.outline_color, linewidth=style.line_width, zorder=10
        )
    if (node.is_proband or node.notation.is_consultand) and node.rotation == 0.0:
        tip_x, tip_y = points[1]
        _draw_proband_marker(ax, tip_x, tip_y, half, style, letter=node.is_proband)


def _draw_notation(
    ax: Axes, node: NodeLayout, style: PedigreeStyle, carriers: CarrierScheme, banded: bool
) -> None:
    """Draw adoption brackets and the letter inside a symbol (P, D or S).

    The letter is white on a filled symbol and gets a white backing on
    carrier sections. A symbol with phenotype bands has text of its own, so
    the letter goes outside, beside its upper-right corner.
    """
    half = style.shape_size / 2
    triangle = symbols.is_pregnancy_loss(node.notation)
    if node.notation.adoption is not None:
        for line in symbols.bracket_lines(node.x, node.y, half, node.sex, triangle):
            ax.plot(
                [p[0] for p in line],
                [p[1] for p in line],
                color=style.line_color,
                linewidth=style.line_width,
                solid_joinstyle="miter",
                zorder=10,
            )
    letter = symbols.inner_letter(node.notation)
    if letter is None or triangle:
        return
    divided = bool(carriers.marks(node)) or (
        carriers.legacy and (node.is_carrier or node.is_obligate_carrier)
    )
    x, y = symbols.letter_position(node.x, node.y, half, node.sex, outside=banded)
    kwargs: dict[str, Any] = {"color": style.line_color}
    if banded:
        kwargs.update(ha="left", va="baseline")
    else:
        kwargs.update(ha="center", va="center_baseline")
        if divided:
            kwargs["bbox"] = _LETTER_BACKING
        elif node.phenotype == Phenotype.AFFECTED:
            kwargs["color"] = text_color_on(style.fill_color)
    ax.text(
        x,
        y,
        letter,
        fontsize=style.font_size,
        fontfamily=style.font_family,
        zorder=11,
        **kwargs,
    )


# White backing behind a letter written over bands or carrier sections.
_LETTER_BACKING = {
    "boxstyle": "round,pad=0.15",
    "facecolor": "white",
    "edgecolor": "none",
    "alpha": 0.9,
}


def _draw_carrier_half_fill(
    ax: Axes,
    cx: float,
    cy: float,
    half: float,
    sex: Sex,
    style: PedigreeStyle,
) -> None:
    """Draw a half-filled shape (left half filled) for carrier individuals.

    Uses a rectangle clipped to the shape outline to fill the left half.
    """
    clip_patch = _get_shape_clip_path(cx, cy, half, sex)
    ax.add_patch(clip_patch)

    # Left-half rectangle, extending past shape bounds vertically to ensure full coverage
    left_rect = mpatches.FancyBboxPatch(
        (cx - half * 1.5, cy - half * 1.5),
        half * 1.5,  # width: from left edge to center
        half * 3,
        boxstyle="square,pad=0",
        facecolor=style.fill_color,
        edgecolor="none",
        zorder=5,
    )
    left_rect.set_clip_path(clip_patch)
    ax.add_patch(left_rect)


def _draw_obligate_carrier_dot(
    ax: Axes,
    cx: float,
    cy: float,
    style: PedigreeStyle,
) -> None:
    """Draw a small filled dot at center for obligate carrier individuals."""
    dot = mpatches.Circle(
        (cx, cy),
        radius=style.carrier_dot_radius,
        facecolor=style.fill_color,
        edgecolor="none",
        zorder=5,
    )
    ax.add_patch(dot)


def _draw_edge(ax: Axes, edge: EdgeLayout, style: PedigreeStyle) -> None:
    """Draw a connection edge.

    Args:
        ax: The Axes to draw on.
        edge: The edge to draw.
        style: Visual style.
    """
    x1, y1 = edge.start
    x2, y2 = edge.end
    color = style.line_color
    lw = style.line_width

    # If waypoints are present, draw a polyline through them
    if edge.waypoints and edge.line_style == LineStyle.SOLID:
        all_pts = [(x1, y1)] + list(edge.waypoints) + [(x2, y2)]
        xs = [p[0] for p in all_pts]
        ys = [p[1] for p in all_pts]
        ax.plot(xs, ys, color=color, linewidth=lw, zorder=1)
        return

    if edge.line_style == LineStyle.SOLID:
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, zorder=1)

    elif edge.line_style == LineStyle.DOUBLE:
        # Two parallel lines
        gap = style.double_line_gap
        if y1 == y2:  # Horizontal
            ax.plot([x1, x2], [y1 + gap, y2 + gap], color=color, linewidth=lw, zorder=1)
            ax.plot([x1, x2], [y1 - gap, y2 - gap], color=color, linewidth=lw, zorder=1)
        else:  # Vertical
            ax.plot([x1 + gap, x2 + gap], [y1, y2], color=color, linewidth=lw, zorder=1)
            ax.plot([x1 - gap, x2 - gap], [y1, y2], color=color, linewidth=lw, zorder=1)

    elif edge.line_style == LineStyle.BROKEN:
        # Line with a gap in the middle and "//" marks
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        gap = style.broken_gap_size

        # Left segment
        ax.plot([x1, mid_x - gap], [y1, mid_y], color=color, linewidth=lw, zorder=1)
        # Right segment
        ax.plot([mid_x + gap, x2], [mid_y, y2], color=color, linewidth=lw, zorder=1)
        # "//" marks
        slash_h = gap * 0.6
        ax.plot(
            [mid_x - gap * 0.3, mid_x + gap * 0.3],
            [mid_y - slash_h, mid_y + slash_h],
            color=color,
            linewidth=lw * 0.7,
            zorder=2,
        )
        ax.plot(
            [mid_x - gap * 0.6, mid_x],
            [mid_y - slash_h, mid_y + slash_h],
            color=color,
            linewidth=lw * 0.7,
            zorder=2,
        )

    elif edge.line_style == LineStyle.DASHED:
        ax.plot(
            [x1, x2],
            [y1, y2],
            color=color,
            linewidth=lw,
            linestyle="--",
            zorder=1,
        )


def _draw_twin_connector(ax: Axes, tc: TwinConnector, style: PedigreeStyle) -> None:
    """Draw twin connection lines from convergence point to each twin.

    Visual elements per twin type:
    - All types: V-lines from convergence point down to each twin's shape
    - MZ (monozygotic): horizontal bar across the V-lines
    - DZ (dizygotic): plain V-lines only
    - Unknown: "?" mark above the convergence point

    Args:
        ax: The Axes to draw on.
        tc: The twin connector geometry.
        style: Visual style parameters.
    """
    cx, cy = tc.convergence_point
    color = style.line_color
    lw = style.line_width

    # Draw V-lines from convergence point to each child
    for child_x, child_y in tc.child_positions:
        ax.plot(
            [cx, child_x],
            [cy, child_y],
            color=color,
            linewidth=lw,
            zorder=1,
        )

    if tc.twin_type == TwinType.MONOZYGOTIC:
        # Horizontal bar across the V-lines at midpoint
        if len(tc.child_positions) >= 2:
            # Bar at 50% of the way from convergence to children
            bar_frac = 0.5
            bar_points = []
            for child_x, child_y in tc.child_positions:
                bx = cx + bar_frac * (child_x - cx)
                by = cy + bar_frac * (child_y - cy)
                bar_points.append((bx, by))
            # Draw bar connecting the midpoints
            bar_xs = [p[0] for p in bar_points]
            bar_ys = [p[1] for p in bar_points]
            ax.plot(bar_xs, bar_ys, color=color, linewidth=lw, zorder=2)

    elif tc.twin_type == TwinType.UNKNOWN and len(tc.child_positions) >= 2:
        # "?" mark in the center of the V, slightly below convergence
        bar_frac = 0.35
        avg_child_x = sum(p[0] for p in tc.child_positions) / len(tc.child_positions)
        avg_child_y = sum(p[1] for p in tc.child_positions) / len(tc.child_positions)
        q_x = cx + bar_frac * (avg_child_x - cx)
        q_y = cy + bar_frac * (avg_child_y - cy)
        ax.text(
            q_x,
            q_y,
            "?",
            fontsize=style.font_size,
            fontfamily=style.font_family,
            fontweight="bold",
            ha="center",
            va="center",
            color=color,
            zorder=3,
        )

    # DZ: just V-lines, nothing extra


def _draw_duplicate_arc(ax: Axes, arc: DuplicateArc, style: PedigreeStyle) -> None:
    """Draw a dotted arc connecting an original node to its duplicate.

    The arc curves above (or below) the straight line between the two
    nodes, creating a visible connection between instances of the same
    individual.

    Args:
        ax: The Axes to draw on.
        arc: The duplicate arc geometry.
        style: Visual style parameters.
    """
    x1, y1 = arc.from_pos
    x2, y2 = arc.to_pos
    color = style.line_color

    # Compute arc as a curved path using a FancyArrowPatch with
    # connectionstyle="arc3" for a smooth curve.
    # The curvature sign is pre-computed so the arc bows away from
    # the tree centre (see layout._compute_duplicate_arcs).
    rad = arc.curvature
    arc_patch = mpatches.FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        connectionstyle=f"arc3,rad={rad}",
        linestyle=":",
        linewidth=style.line_width,
        edgecolor=color,
        facecolor="none",
        arrowstyle="-",
        zorder=0,
    )
    ax.add_patch(arc_patch)


def _draw_label(ax: Axes, node: NodeLayout, style: PedigreeStyle) -> None:
    """Draw the sample_id label and any annotations below a node.

    A rotated node (circular layout) needs two adjustments. Its shape can
    reach further down than an upright one — a square at 45° by √2 × half —
    so the label moves down only as far as needed to stay clear of it. And
    lines in a radial layout run in every direction, so its labels get a
    white backing that keeps a line passing behind from striking through.
    """
    half = style.shape_size / 2
    if node.label_angle is not None:
        _draw_radial_label(ax, node, style)
        return
    text_kwargs: dict[str, Any] = {}
    if node.rotation == 0.0:
        y_cursor = node.y - half - style.label_offset - node.marker_depth
    else:
        extent = _shape_lower_extent(node.sex, half, node.rotation)
        y_cursor = node.y - max(
            half + style.label_offset,
            extent + style.label_offset * _ROTATED_LABEL_MIN_GAP,
        )
        text_kwargs["bbox"] = _ROTATED_LABEL_BACKING
    ax.text(
        node.x,
        y_cursor,
        node.label,
        fontsize=style.font_size,
        fontfamily=style.font_family,
        ha="center",
        va="top",
        color=style.line_color,
        **text_kwargs,
    )
    # Notation notes ("SB", "AMAB"...), then annotations, below the ID label
    for text, kwargs in _lines_below_label(node, style):
        y_cursor -= symbols.TEXT_LINE_STEP
        ax.text(
            node.x,
            y_cursor,
            text,
            fontsize=style.font_size * 0.7,
            fontfamily=style.font_family,
            ha="center",
            va="top",
            **kwargs,
            **text_kwargs,
        )


def _lines_below_label(node: NodeLayout, style: PedigreeStyle) -> list[tuple[str, dict[str, Any]]]:
    """Text below the ID label: notation notes, then italic grey annotations."""
    lines: list[tuple[str, dict[str, Any]]] = [
        (note, {"color": style.line_color}) for note in node.notes
    ]
    lines += [(a, {"color": "#666666", "fontstyle": "italic"}) for a in node.annotations]
    return lines


def _draw_radial_label(ax: Axes, node: NodeLayout, style: PedigreeStyle) -> None:
    """Draw the label and annotations of a circular-layout node, outward.

    Lines reach a node of the circular layout from the centre, so the label
    is set on the far side of the shape along the node's radial direction,
    aligned so the text grows away from it. Annotations follow further out.
    """
    assert node.label_angle is not None
    half = style.shape_size / 2
    c, s = math.cos(node.label_angle), math.sin(node.label_angle)
    ha = "left" if c > 0.35 else "right" if c < -0.35 else "center"
    va = "bottom" if s > 0.35 else "top" if s < -0.35 else "center"
    reach = _shape_boundary_radius(node.sex, half, node.rotation, node.label_angle)
    d = reach + style.label_offset + node.marker_depth
    lines: list[tuple[str, float, dict[str, Any]]] = [
        (node.label, style.font_size, {"color": style.line_color})
    ]
    lines += [(t, style.font_size * 0.7, kw) for t, kw in _lines_below_label(node, style)]
    for text, size, kwargs in lines:
        ax.text(
            node.x + d * c,
            node.y + d * s,
            text,
            fontsize=size,
            fontfamily=style.font_family,
            ha=ha,
            va=va,
            bbox=_ROTATED_LABEL_BACKING,
            **kwargs,
        )
        d += 0.25 if va == "center" and ha == "center" else 0.3


# A rotated node's label keeps at least this fraction of label_offset between
# the lowest point of its shape and the top of the text.
_ROTATED_LABEL_MIN_GAP = 0.4

# White backing behind rotated nodes' labels, so lines do not strike through.
_ROTATED_LABEL_BACKING = {
    "boxstyle": "round,pad=0.12",
    "facecolor": "white",
    "edgecolor": "none",
    "alpha": 0.85,
}


def _get_shape_clip_path(cx: float, cy: float, half: float, sex: Sex) -> mpatches.Patch:
    """Create a clip patch for the given shape type."""
    if sex == Sex.MALE:
        return mpatches.FancyBboxPatch(
            (cx - half, cy - half),
            half * 2,
            half * 2,
            boxstyle="square,pad=0",
            facecolor="none",
            edgecolor="none",
        )
    elif sex == Sex.FEMALE:
        return mpatches.Circle((cx, cy), radius=half, facecolor="none", edgecolor="none")
    else:
        return mpatches.RegularPolygon(
            (cx, cy),
            numVertices=4,
            radius=half * _DIAMOND_SCALE,
            orientation=0,
            facecolor="none",
            edgecolor="none",
        )


def _draw_banded_node(
    ax: Axes,
    node: NodeLayout,
    style: PedigreeStyle,
    band_names: list[str],
    display_modes: dict[str, PhenotypeDisplayMode] | None = None,
    outline_color: str | None = None,
    pheno_render: PhenotypeRender | None = None,
    *,
    carriers: CarrierScheme | None = None,
) -> None:
    """Draw a node with horizontal color bands for multiple phenotypes.

    Each band occupies an equal vertical slice of the shape. Binary phenotypes
    show full color (1) or white (0). Continuous phenotypes show a gradient
    from white (low/0) to full color (high/1).

    Args:
        ax: The matplotlib Axes to draw on.
        node: Node geometry and phenotype values.
        style: Visual style parameters.
        band_names: Phenotype columns to draw as bands, top to bottom.
        display_modes: Map of phenotype name → declared display mode.
        outline_color: Shape border color from an `o_` phenotype, if any.
        pheno_render: Per-figure phenotype display resolution.
        carriers: Per-figure carrier notation. Each carrier source adds a
            band below the phenotype bands, patterned for carriers.
    """
    if pheno_render is None:
        pheno_render = PhenotypeRender()
    if carriers is None:
        carriers = CarrierScheme()
    half = style.shape_size / 2
    cx, cy = node.x, node.y
    edgecolor = style.outline_color
    linewidth = style.shape_line_width
    linestyle = style.phantom_line_style if node.is_phantom else "-"
    carrier_keys = list(carriers.keys)
    n_bands = len(band_names) + len(carrier_keys)

    if not band_names:
        # No phenotype band: affection status and carrier sections only.
        _draw_node(ax, node, style, carriers)
        return

    # Create clip patch and add to axes
    clip_patch = _get_shape_clip_path(cx, cy, half, node.sex)
    ax.add_patch(clip_patch)

    # Draw white background first (covers the whole shape)
    bg = _get_shape_clip_path(cx, cy, half, node.sex)
    bg.set_facecolor("white")
    bg.set_edgecolor("none")
    bg.set_zorder(3)
    ax.add_patch(bg)

    # Band geometry: divide the shape height into n_bands equal parts
    band_height = (half * 2) / n_bands
    shape_top = cy + half

    if display_modes is None:
        display_modes = {}

    # Check if ALL bands are NA — if so, show a single centered "NA"
    all_na = not carrier_keys and all_bands_na(node, band_names)

    for i, name in enumerate(band_names):
        cell = band_cell(node, name, style, pheno_render, display_modes)

        # Draw band rectangle (will be clipped to shape)
        band_top = shape_top - i * band_height
        band_rect = mpatches.FancyBboxPatch(
            (cx - half * 1.5, band_top - band_height),
            half * 3,
            band_height,
            boxstyle="square,pad=0",
            facecolor=cell.color,
            edgecolor="none",
            zorder=4,
        )
        band_rect.set_clip_path(clip_patch)
        ax.add_patch(band_rect)

        # Draw band separator line (thin grey)
        if i > 0:
            ax.plot(
                [cx - half * 1.5, cx + half * 1.5],
                [band_top, band_top],
                color=SEPARATOR_COLOR,
                linewidth=0.5,
                zorder=5,
                clip_path=clip_patch,
                clip_on=True,
            )

        # Text inside the band; a symbol with no data at all gets a single
        # centred "NA" instead (drawn below).
        if cell.text is not None and not (cell.is_na and all_na):
            ax.text(
                cx,
                band_top - band_height / 2,
                cell.text,
                fontsize=style.font_size * 0.75,
                fontfamily=style.font_family,
                fontweight="bold",
                ha="center",
                va="center",
                color=cell.text_color,
                zorder=6,
            )

    for j, key in enumerate(carrier_keys):
        band_top = shape_top - (len(band_names) + j) * band_height
        state = carrier_state(node, key)
        hatch = None
        if state in (CARRIER, OBLIGATE_CARRIER):
            hatch = _hatch(carriers.pattern(key), state == OBLIGATE_CARRIER)
        carrier_rect = mpatches.Rectangle(
            (cx - half * 1.5, band_top - band_height),
            half * 3,
            band_height,
            facecolor="white",
            edgecolor=style.fill_color if hatch else "none",
            linewidth=0,
            hatch=hatch,
            zorder=4,
        )
        carrier_rect.set_clip_path(clip_patch)
        ax.add_patch(carrier_rect)
        ax.plot(
            [cx - half * 1.5, cx + half * 1.5],
            [band_top, band_top],
            color="#cccccc",
            linewidth=0.5,
            zorder=5,
            clip_path=clip_patch,
            clip_on=True,
        )
        if state is None:
            ax.text(
                cx,
                band_top - band_height / 2,
                "NA",
                fontsize=style.font_size * 0.75,
                fontfamily=style.font_family,
                fontweight="bold",
                ha="center",
                va="center",
                color="#999999",
                zorder=6,
            )

    # If all bands are NA, draw a single centered "NA" label
    if all_na and band_names:
        ax.text(
            cx,
            cy,
            "NA",
            fontsize=style.font_size * 0.85,
            fontfamily=style.font_family,
            fontweight="bold",
            ha="center",
            va="center",
            color="#999999",
            zorder=6,
        )

    # Draw shape outline on top (use outline_color from o_ phenotype if set)
    oc = outline_color if outline_color else edgecolor
    olw = linewidth * 2.5 if outline_color else linewidth
    if node.sex == Sex.MALE:
        _draw_square(ax, cx, cy, half, "none", oc, olw, linestyle)
    elif node.sex == Sex.FEMALE:
        _draw_circle(ax, cx, cy, half, "none", oc, olw, linestyle)
    else:
        _draw_diamond(ax, cx, cy, half, "none", oc, olw, linestyle)

    # Legacy carrier notation: half-fill and dot, drawn on top of bands
    if carriers.legacy:
        if node.is_carrier:
            _draw_carrier_half_fill(ax, cx, cy, half, node.sex, style)
        if node.is_obligate_carrier:
            _draw_obligate_carrier_dot(ax, cx, cy, style)

    # Deceased marker
    if node.is_deceased:
        _draw_deceased_marker(ax, cx, cy, half, style)

    # Proband marker. A rotated node gets its marker after rotation instead,
    # so the "P" stays upright (see _render_layout).
    if (node.is_proband or node.notation.is_consultand) and node.rotation == 0.0:
        _draw_proband_marker(ax, cx - half, cy - half, half, style, letter=node.is_proband)


def _draw_pgs_ridge(
    ax: Axes,
    node: NodeLayout,
    style: PedigreeStyle,
    pgs_names: list[str],
) -> None:
    """Draw a single gaussian distribution with multiple colored markers.

    One shared curve, each PGS trait's individual score shown as a colored
    dashed line + dot at the corresponding y-position.
    """
    half = style.shape_size / 2
    cx, cy = node.x, node.y

    # Filter to PGS columns that have data for this node
    active = active_pgs(node, pgs_names)
    if not active:
        return

    box_width = half * 0.8
    box_height = half * 2
    gap = half * 0.15
    box_left = cx - half - box_width - gap
    box_bottom = cy - half

    # Shared box outline
    box_rect = mpatches.FancyBboxPatch(
        (box_left, box_bottom),
        box_width,
        box_height,
        boxstyle="round,pad=0.02",
        facecolor=(1, 1, 1, 0.5),
        edgecolor="#cccccc",
        linewidth=0.5,
        zorder=3,
    )
    ax.add_patch(box_rect)

    # Single gaussian curve (grey/neutral)
    n_points = 50
    y_vals = np.linspace(0, 1, n_points)
    density = np.exp(-0.5 * ((y_vals - 0.5) / 0.18) ** 2)
    density = density / density.max()

    plot_x = box_left + density * box_width * 0.8 + box_width * 0.1
    plot_y = box_bottom + y_vals * box_height

    ax.plot(plot_x, plot_y, color="#aaaaaa", linewidth=0.6, zorder=4)
    ax.fill_betweenx(
        plot_y,
        box_left + box_width * 0.1,
        plot_x,
        alpha=0.08,
        color="#999999",
        zorder=4,
    )

    # One colored marker per PGS trait
    for pgs_name, norm_val in active:
        marker_color = style.pgs_colors.get(pgs_name, DEFAULT_PGS_COLOR)

        marker_y = box_bottom + norm_val * box_height
        ax.plot(
            [box_left, box_left + box_width],
            [marker_y, marker_y],
            color=marker_color,
            linewidth=1.0,
            linestyle="--",
            zorder=5,
        )
        ax.plot(
            box_left + box_width * 0.5,
            marker_y,
            "o",
            color=marker_color,
            markersize=3,
            zorder=6,
        )

    # Single "PGS" label
    ax.text(
        box_left + box_width / 2,
        box_bottom - half * 0.15,
        "PGS",
        fontsize=style.font_size * 0.6,
        fontfamily=style.font_family,
        ha="center",
        va="top",
        color="#666666",
        zorder=4,
    )


def _legend_swatch(
    ax: Axes,
    x: float,
    y_center: float,
    width: float,
    height: float,
    facecolor: str | tuple[float, ...],
    *,
    edgecolor: str = LEG_SWATCH_EDGE,
    linewidth: float | None = 0.5,
) -> None:
    """Add one rectangular swatch centred vertically on `y_center`.

    Args:
        ax: Legend axes.
        x: Left edge.
        y_center: Vertical centre.
        width: Swatch width.
        height: Swatch height.
        facecolor: Fill color, or "none".
        edgecolor: Border color, or "none" for no border.
        linewidth: Border width; None leaves matplotlib's default.
    """
    kwargs = {} if linewidth is None else {"linewidth": linewidth}
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (x, y_center - height / 2),
            width,
            height,
            boxstyle="square,pad=0",
            facecolor=facecolor,
            edgecolor=edgecolor,
            **kwargs,
        )
    )


def _legend_caption(
    ax: Axes,
    x: float,
    y: float,
    text: str,
    m: LegendMetrics,
    style: PedigreeStyle,
    *,
    ha: str = "left",
    color: str = LEG_CAPTION_COLOR,
    size_factor: float = 1.0,
) -> None:
    """Add a small grey caption (column header, scale end, category name)."""
    ax.text(
        x,
        y,
        text,
        fontsize=m.label_fs * size_factor,
        fontfamily=style.font_family,
        va="top",
        ha=ha,
        color=color,
    )


def _legend_row_label(
    ax: Axes, y_center: float, text: str, x: float, m: LegendMetrics, style: PedigreeStyle
) -> None:
    """Add the description text of a legend row."""
    ax.text(x, y_center, text, fontsize=m.fs, fontfamily=style.font_family, va="center")


def _add_legend_axes(fig: Figure, m: LegendMetrics) -> Axes:
    """Grow the figure downward by the panel height and add the legend axes.

    The axes already on the figure keep their size in inches; the legend
    gets a strip of exactly `m.total_h` vertical units underneath, and one
    legend unit is `m.unit_x` by `m.unit_y` inches.

    Args:
        fig: Figure to attach the legend axes to.
        m: Panel metrics.

    Returns:
        The legend axes, with data coordinates in legend units.
    """
    width, height = fig.get_size_inches()
    panel_h = m.total_h * m.unit_y
    new_height = height + panel_h
    for other in fig.axes:
        x0, y0, w, h = other.get_position().bounds
        other.set_position((x0, (y0 * height + panel_h) / new_height, w, h * height / new_height))
    fig.set_size_inches(width, new_height)

    left, ax_w = 0.03, 0.94
    ax = fig.add_axes((left, 0.0, ax_w, panel_h / new_height))
    ax.set_xlim(0, ax_w * width / m.unit_x)
    ax.set_ylim(0, m.total_h)
    ax.axis("off")
    return ax


def _draw_legend_header(ax: Axes, m: LegendMetrics, style: PedigreeStyle) -> float:
    """Draw the panel title and the sex symbols.

    Args:
        ax: Legend axes.
        m: Panel metrics.
        style: Visual style.

    Returns:
        The y cursor just below the header.
    """
    y = m.total_h - LEG_PAD
    ax.text(
        LEG_X0,
        y,
        "Legend",
        fontsize=max(m.fs + 1, style.font_size * 0.7),
        fontfamily=style.font_family,
        fontweight="bold",
        va="top",
    )
    y -= 1.0

    # Male / female reference symbols, drawn as plot markers so they stay
    # circular/square regardless of the axes aspect ratio.
    sym_gap = 0.25
    for x, marker, label in [(0.6, "s", "= Male"), (2.5, "o", "= Female")]:
        ax.plot(
            x,
            y,
            marker=marker,
            markersize=m.marker_size,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=1.0,
            linestyle="none",
        )
        ax.text(
            x + sym_gap,
            y,
            label,
            fontsize=m.fs,
            fontfamily=style.font_family,
            va="center",
        )

    return y - 0.8


def _draw_legend_border(ax: Axes, m: LegendMetrics) -> None:
    """Frame the panel content, which is drawn first so it can be measured.

    The border hugs the content horizontally rather than spanning the whole
    figure, so a wide pedigree does not get a wide, mostly empty panel.
    """
    artists = [*ax.texts, *ax.patches, *ax.lines, *ax.collections]
    # Before the first draw, matplotlib 3.7 reports a collection's extent as
    # the null box (infinite bounds), which would make the union infinite.
    extents = [
        e for e in (a.get_window_extent() for a in artists) if np.isfinite([e.x0, e.x1]).all()
    ]
    right = LEG_X0 + 4.0
    if extents:
        content = Bbox.union(extents)
        right = max(right, ax.transData.inverted().transform((content.x1, 0))[0])
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (0.1, 0.1),
            right + LEG_X0 - 0.1,
            m.total_h - 0.2,
            boxstyle="round,pad=0.1",
            facecolor="white",
            edgecolor=LEG_SWATCH_EDGE,
            linewidth=0.8,
            zorder=-1,
            clip_on=False,
        )
    )


def _draw_legend_binary_headers(
    ax: Axes, y: float, m: LegendMetrics, style: PedigreeStyle
) -> float:
    """Label the three binary swatch columns once, above all rows.

    Returns:
        The y cursor below the headers.
    """
    for i, text in enumerate(("Affected", "Unaffected", "Unknown")):
        _legend_caption(ax, m.swatch_x(i) + m.sw_w / 2, y, text, m, style, ha="center")
    return y - LEG_BINARY_HEADER_H


def _draw_legend_gradient_bar(
    ax: Axes,
    y_center: float,
    colors: list[tuple[float, ...]],
    low_label: str,
    high_label: str,
    m: LegendMetrics,
    style: PedigreeStyle,
) -> None:
    """Draw a gradient bar from a list of colors, with end captions."""
    step_w = LEG_GRAD_W / len(colors)
    for i, color in enumerate(colors):
        # Steps carry no edge of their own, so the bar reads as continuous.
        _legend_swatch(
            ax,
            LEG_X0 + i * step_w,
            y_center,
            step_w,
            m.sw_h,
            color,
            edgecolor="none",
            linewidth=None,
        )
    # Single border drawn over the seamless fill.
    _legend_swatch(ax, LEG_X0, y_center, LEG_GRAD_W, m.sw_h, "none")

    caption_y = y_center - m.sw_h / 2 - 0.05
    _legend_caption(ax, LEG_X0, caption_y, low_label, m, style)
    _legend_caption(ax, LEG_X0 + LEG_GRAD_W, caption_y, high_label, m, style, ha="right")


def _draw_legend_categorical_row(
    ax: Axes,
    y_center: float,
    categories: dict[str, str],
    m: LegendMetrics,
    style: PedigreeStyle,
) -> None:
    """Draw one small swatch per category value, each captioned."""
    x = LEG_X0
    for value, color in categories.items():
        _legend_swatch(ax, x, y_center, m.sw_w * 0.8, m.sw_h, color)
        _legend_caption(
            ax,
            x + m.sw_w * 0.4,
            y_center - m.sw_h / 2 - 0.05,
            value,
            m,
            style,
            ha="center",
            color=LEG_CATEGORY_COLOR,
            size_factor=0.8,
        )
        x += m.sw_w + m.sw_gap * 0.5


def _draw_legend_binary_row(ax: Axes, y_center: float, color_hex: str, m: LegendMetrics) -> None:
    """Draw the affected / unaffected / unknown swatches of a binary column."""
    for i, fill in enumerate((color_hex, "white", UNKNOWN_FILL)):
        _legend_swatch(ax, m.swatch_x(i), y_center, m.sw_w, m.sw_h, fill)


def _draw_legend_outline_row(
    ax: Axes, y_center: float, row: LegendRow, m: LegendMetrics, style: PedigreeStyle
) -> None:
    """Draw the thick-bordered swatch standing for an `o_` column."""
    _legend_swatch(
        ax, LEG_X0, y_center, m.sw_w, m.sw_h, "white", edgecolor=row.color, linewidth=2.5
    )
    _legend_row_label(ax, y_center, row.label, m.icon_text_x, m, style)


def _draw_legend_pgs_row(
    ax: Axes,
    y_center: float,
    row: LegendRow,
    m: LegendMetrics,
    style: PedigreeStyle,
) -> None:
    """Draw the shared distribution-curve icon and its marker lines."""
    xs = np.linspace(0, 1, 30)
    dens = np.exp(-0.5 * ((xs - 0.5) / PGS_CURVE_SIGMA) ** 2)
    dens = dens / dens.max() * m.sw_h * 0.9
    plot_x = LEG_X0 + xs * 0.9
    baseline = y_center - m.sw_h / 2
    plot_y = baseline + dens

    ax.plot(plot_x, plot_y, color="#aaaaaa", linewidth=0.6)
    ax.fill_between(plot_x, baseline, plot_y, alpha=0.08, color=FALLBACK_COLOR)

    for k, color in enumerate(row.pgs.values()):
        marker_x = LEG_X0 + (0.3 + k * 0.2) * 0.9
        ax.plot(
            [marker_x, marker_x],
            [baseline, y_center + m.sw_h / 2],
            color=color,
            linewidth=0.8,
            linestyle="--",
        )

    _legend_row_label(ax, y_center, row.label, m.icon_text_x, m, style)


def _draw_legend_carrier_row(
    ax: Axes,
    y_center: float,
    row: LegendRow,
    m: LegendMetrics,
    style: PedigreeStyle,
) -> None:
    """Draw the carrier or obligate-carrier icon and its caption."""
    icon_cx = 0.9
    ax.plot(
        icon_cx,
        y_center,
        marker="s",
        markersize=9 * m.scale,
        markerfacecolor="white",
        markeredgecolor="black",
        markeredgewidth=1.0,
        linestyle="none",
    )
    if row.obligate:
        # Centre dot.
        ax.plot(
            icon_cx,
            y_center,
            marker="o",
            markersize=3 * m.scale,
            markerfacecolor="black",
            markeredgecolor="none",
            linestyle="none",
        )
    else:
        # Left half filled, approximated with a smaller solid marker.
        ax.plot(
            icon_cx - m.sw_h * 0.45 * 0.25,
            y_center,
            marker="s",
            markersize=5 * m.scale,
            markerfacecolor="black",
            markeredgecolor="none",
            linestyle="none",
            clip_on=True,
        )

    _legend_row_label(ax, y_center, row.label, m.icon_text_x, m, style)


def _draw_legend_hatch_row(
    ax: Axes,
    y_center: float,
    hatch: str,
    text: str,
    m: LegendMetrics,
    style: PedigreeStyle,
) -> None:
    """Draw a patterned swatch for one carrier source (NSGC 2022 notation)."""
    ax.add_patch(
        mpatches.Rectangle(
            (LEG_X0, y_center - m.sw_h / 2),
            m.sw_w,
            m.sw_h,
            facecolor="white",
            edgecolor=style.fill_color,
            linewidth=0.5,
            hatch=hatch,
        )
    )
    _legend_row_label(ax, y_center, text, m.icon_text_x, m, style)


def _draw_legend(
    fig: Figure,
    phenotype_names: list[str],
    style: PedigreeStyle,
    *,
    carriers: CarrierScheme | None = None,
    display_modes: dict[str, PhenotypeDisplayMode] | None = None,
    pheno_render: PhenotypeRender | None = None,
) -> None:
    """Draw the legend explaining the phenotype encoding and carrier notation.

    The rows come from `legend_rows`, which the SVG legend reads too, and
    each row's appearance from `pheno_render`, the same resolution the
    bands were drawn from.

    Args:
        fig: The figure to attach the legend to.
        phenotype_names: Ordered phenotype column names.
        style: Visual style.
        carriers: Per-figure carrier notation.
        display_modes: Modes declared by column prefixes.
        pheno_render: Per-figure phenotype display resolution.
    """
    rows = legend_rows(
        phenotype_names,
        style,
        pheno_render or PhenotypeRender(),
        display_modes or {},
        carriers or CarrierScheme(),
    )
    if not rows:
        return

    binary_headers = any(r.kind == "binary" for r in rows)
    m = legend_metrics(len(rows), style, binary_headers=binary_headers)
    ax = _add_legend_axes(fig, m)
    y_top = _draw_legend_header(ax, m, style)

    if binary_headers:
        y_top = _draw_legend_binary_headers(ax, y_top, m, style)

    for i, row in enumerate(rows):
        y_center = m.row_center(i, y_top)
        if row.kind == "gradient":
            _draw_legend_gradient_bar(ax, y_center, list(row.colors), row.low, row.high, m, style)
        elif row.kind == "categorical":
            _draw_legend_categorical_row(ax, y_center, dict(row.categories), m, style)
        elif row.kind == "binary":
            _draw_legend_binary_row(ax, y_center, row.color, m)
        elif row.kind == "outline":
            _draw_legend_outline_row(ax, y_center, row, m, style)
        elif row.kind == "pgs":
            _draw_legend_pgs_row(ax, y_center, row, m, style)
        elif row.kind == "legacy_carrier":
            _draw_legend_carrier_row(ax, y_center, row, m, style)
        else:
            _draw_legend_hatch_row(
                ax, y_center, _hatch(row.pattern, row.obligate), row.label, m, style
            )
        if row.kind in ("gradient", "categorical", "binary"):
            # Colors are self-explanatory, so the row label is just the name.
            _legend_row_label(ax, y_center, row.label, m.text_x, m, style)

    _draw_legend_border(ax, m)


def _finalize_axes(ax: Axes) -> None:
    """Clean up axes: equal aspect, no ticks, auto-scale."""
    ax.set_aspect("equal")
    ax.autoscale()
    ax.margins(0.15)
    ax.axis("off")
