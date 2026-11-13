"""Interactive SVG renderer for pedigree diagrams.

Generates standalone SVG files with embedded JavaScript for:
- Tooltips on hover (sample_id, sex, phenotype, DOB, extra phenotypes)
- Click to highlight connected relatives (parents, children, siblings)
- Click on background to reset
- Mouse wheel zoom and drag to pan

Uses the same FamilyLayout pipeline as the matplotlib renderer, and the
same phenotype and carrier decisions (`pypedigree.phenotypes`,
`pypedigree.carriers`), but writes SVG directly: from matplotlib it only
uses `matplotlib.colors`, through `pypedigree.phenotypes`, to read colors.
"""

from __future__ import annotations

import html
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from pypedigree import symbols
from pypedigree.carriers import (
    CARRIER,
    LEGACY_CARRIER_KEY,
    OBLIGATE_CARRIER,
    CarrierScheme,
    carrier_columns,
    carrier_state,
    check_notation,
    resolve_carrier_scheme,
)
from pypedigree.config import PedigreeStyle, apply_pgs_spacing, resolve_style
from pypedigree.layout import compute_layout
from pypedigree.models import (
    DuplicateArc,
    EdgeLayout,
    Family,
    FamilyLayout,
    Individual,
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
    LEG_UNIT_X_EM,
    LEG_UNIT_Y_EM,
    LEG_X0,
    NA_TEXT_COLOR,
    SEPARATOR_COLOR,
    UNKNOWN_FILL,
    LegendRow,
    PhenotypeRender,
    active_pgs,
    all_bands_na,
    band_cell,
    legend_metrics,
    legend_rows,
    outline_color,
    pgs_curve,
    pregnancy_loss_fill,
    resolve_phenotype_render,
    text_color_on,
    to_hex,
    validate_style_colors,
)


def render_interactive_svg(
    pedigree: Pedigree,
    style: PedigreeStyle | None = None,
    output: str | Path | None = None,
    show_labels: bool = True,
    show_legend: bool = True,
) -> str:
    """Render pedigree as interactive SVG with tooltips and click highlighting.

    The drawing matches `render_pedigree`: phenotype bands, outline colors,
    PGS curves, annotations, carrier notation and the legend all come from
    the same per-figure decisions (`pypedigree.phenotypes`,
    `pypedigree.carriers`). Several families are drawn side by side.

    Args:
        pedigree: The parsed pedigree data.
        style: Visual style parameters. Uses defaults if None.
        output: If provided, write the SVG to this path.
        show_labels: Whether to draw sample_id labels under shapes.
        show_legend: Whether to draw the phenotype and carrier legend.

    Returns:
        The SVG content as a string.
    """
    if style is None:
        style = PedigreeStyle()
    validate_style_colors(style)
    names = list(pedigree.phenotype_names or [])
    modes = pedigree.phenotype_display_modes
    # Same style resolution as render_pedigree, so both outputs share colors
    # and spacing. Each renderer applies the PGS spacing exactly once.
    style = apply_pgs_spacing(resolve_style(style, names), names)
    check_notation(style.carrier_notation)
    families = list(pedigree.families.values())
    carriers = resolve_carrier_scheme(
        (ind for fam in families for ind in fam.individuals.values()),
        carrier_columns(names, modes),
        style.carrier_notation,
    )
    render = resolve_phenotype_render(pedigree, names, modes, style)
    banded = bool(names and style.phenotype_colors and render.effective_modes)
    ctx = _SvgContext(
        style=style,
        carriers=carriers,
        render=render,
        display_modes=modes,
        band_names=render.band_names(names) if banded else [],
        outline_names=render.outline_names(names) if banded else [],
        pgs_names=render.pgs_names(names) if banded else [],
        banded=banded,
        show_labels=show_labels,
    )

    if not families:
        return _write(_wrap_svg("", (0.0, 0.0, 100.0, 100.0), style), output)

    parts: list[str] = []
    if carriers.keys:
        parts.append(_svg_carrier_patterns(carriers, style))

    extents: list[_Extent] = []
    cursor = 0.0
    titled = len(families) > 1
    for index, family in enumerate(families):
        layout = _flip_vertical(compute_layout(family, style))
        x0, y0, x1, y1 = _family_extent(layout, ctx)
        if index > 0:
            # Side by side, left to right, like render_pedigree's subplots.
            dx = cursor - x0

            def shift(x: float, y: float, dx: float = dx) -> tuple[float, float]:
                return (x + dx, y)

            layout = _map_layout(layout, shift)
            x0, x1 = x0 + dx, x1 + dx
        if titled:
            title_fs = (style.font_size + 2) * _PT
            parts.append(
                _svg_text(
                    (x0 + x1) / 2,
                    y0 - title_fs * 0.5,
                    f"Family: {family.family_id}",
                    title_fs,
                    style,
                    anchor="middle",
                    baseline="auto",
                    cls="family-title",
                )
            )
            y0 -= title_fs * 2
        parts.append(_render_family_svg(layout, family, _build_adjacency(family), ctx, f"f{index}"))
        extents.append((x0, y0, x1, y1))
        cursor = x1 + style.shape_size * 3

    bounds = _union(extents)
    show_pheno_legend = show_legend and names and style.phenotype_colors
    show_carrier_legend = show_legend and (carriers.has_legacy_carriers or bool(carriers.keys))
    if show_pheno_legend or show_carrier_legend:
        rows = legend_rows(names, style, render, modes, carriers)
        if rows:
            legend, legend_extent = _svg_legend(
                rows, style, carriers, bounds[0], bounds[3] + style.shape_size / 2
            )
            parts.append(legend)
            bounds = _union([bounds, legend_extent])

    return _write(_wrap_svg("\n".join(parts), bounds, style), output)


# Layout units per typographic point: font sizes and line widths, given in
# points in PedigreeStyle, are scaled by this in the SVG.
_PT = 0.02

_Extent = tuple[float, float, float, float]  # x0, y0, x1, y1 (y down)


@dataclass(frozen=True)
class _SvgContext:
    """Per-figure decisions every node of the SVG is drawn from."""

    style: PedigreeStyle
    carriers: CarrierScheme
    render: PhenotypeRender
    display_modes: Mapping[str, PhenotypeDisplayMode]
    band_names: list[str]
    outline_names: list[str]
    pgs_names: list[str]
    banded: bool
    show_labels: bool


def _write(svg: str, output: str | Path | None) -> str:
    if output is not None:
        Path(output).write_text(svg, encoding="utf-8")
    return svg


def _union(extents: Sequence[_Extent]) -> _Extent:
    return (
        min(e[0] for e in extents),
        min(e[1] for e in extents),
        max(e[2] for e in extents),
        max(e[3] for e in extents),
    )


def _map_layout(
    layout: FamilyLayout, f: Callable[[float, float], tuple[float, float]]
) -> FamilyLayout:
    """Apply a point transformation to every coordinate of a layout."""

    def p(pt: tuple[float, float]) -> tuple[float, float]:
        return f(pt[0], pt[1])

    nodes = []
    for n in layout.nodes:
        x, y = f(n.x, n.y)
        nodes.append(replace(n, x=x, y=y))
    return replace(
        layout,
        nodes=nodes,
        edges=[
            replace(e, start=p(e.start), end=p(e.end), waypoints=[p(w) for w in e.waypoints])
            for e in layout.edges
        ],
        twin_connectors=[
            replace(
                tc,
                convergence_point=p(tc.convergence_point),
                child_positions=[p(c) for c in tc.child_positions],
            )
            for tc in layout.twin_connectors
        ],
        duplicate_arcs=[
            replace(a, from_pos=p(a.from_pos), to_pos=p(a.to_pos)) for a in layout.duplicate_arcs
        ],
    )


def _flip_vertical(layout: FamilyLayout) -> FamilyLayout:
    """Mirror a layout top to bottom for SVG coordinates.

    The layout puts generation 0 at y = 0 and descendants at negative y
    (matplotlib's axis points up); SVG's y axis points down, so without
    this the founders were drawn at the bottom. `_svg_duplicate_arc` bows
    to the opposite side of matplotlib's arc3 for the same curvature, which
    the mirroring compensates exactly.
    """
    return _map_layout(layout, lambda x, y: (x, -y))


def _family_extent(layout: FamilyLayout, ctx: _SvgContext) -> _Extent:
    """Bounding box of everything drawn for one family."""
    style = ctx.style
    half = style.shape_size / 2
    reach = half * 1.2  # the diamond is the widest shape
    xs = [n.x for n in layout.nodes]
    ys = [n.y for n in layout.nodes]
    for e in layout.edges:
        for x, y in (e.start, e.end, *e.waypoints):
            xs.append(x)
            ys.append(y)
    if not xs:
        return (0.0, 0.0, 1.0, 1.0)
    x0, x1 = min(xs) - reach, max(xs) + reach
    y0, y1 = min(ys) - reach, max(ys) + reach
    if ctx.pgs_names:
        x0 -= half * 0.95  # PGS box on the left of each symbol
    if ctx.show_labels and layout.nodes:
        label_bottom = max(
            n.y
            + half
            + style.label_offset
            + n.marker_depth
            + style.font_size * _PT * 1.2
            + (len(n.notes) + len(n.annotations)) * _ANNOTATION_STEP
            for n in layout.nodes
        )
        y1 = max(y1, label_bottom)
    return (x0, y0, x1, y1)


def _build_adjacency(family: Family) -> dict[str, dict[str, list[str]]]:
    """Build adjacency data for each individual (parents, children, siblings).

    Returns:
        Dict mapping sample_id -> {"parents": [...], "children": [...], "siblings": [...]}.
    """
    adj: dict[str, dict[str, list[str]]] = {}

    for sid in family.individuals:
        adj[sid] = {"parents": [], "children": [], "siblings": []}

    for ind in family.individuals.values():
        if ind.father_id and ind.father_id in family.individuals:
            adj[ind.sample_id]["parents"].append(ind.father_id)
        if ind.mother_id and ind.mother_id in family.individuals:
            adj[ind.sample_id]["parents"].append(ind.mother_id)

    for union in family.unions:
        p1, p2 = union.partner_ids
        for cid in union.children_ids:
            if cid in adj:
                if p1 in family.individuals:
                    adj[p1]["children"].append(cid)
                if p2 in family.individuals:
                    adj[p2]["children"].append(cid)
        # Siblings: all children in same union
        for i, c1 in enumerate(union.children_ids):
            for c2 in union.children_ids[i + 1 :]:
                if c1 in adj and c2 in adj:
                    adj[c1]["siblings"].append(c2)
                    adj[c2]["siblings"].append(c1)

    return adj


def _render_family_svg(
    layout: FamilyLayout,
    family: Family,
    adjacency: dict[str, dict[str, list[str]]],
    ctx: _SvgContext,
    prefix: str,
) -> str:
    """Render one family's layout as SVG elements.

    Args:
        prefix: Makes the element ids of this family unique in the document
            (sample ids may repeat across families, and need not be valid
            XML ids).
    """
    parts: list[str] = []
    style = ctx.style

    for edge in layout.edges:
        parts.append(_svg_edge(edge, style))
    for tc in layout.twin_connectors:
        parts.append(_svg_twin_connector(tc, style))
    for arc in layout.duplicate_arcs:
        parts.append(_svg_duplicate_arc(arc, style))

    empty: dict[str, list[str]] = {"parents": [], "children": [], "siblings": []}
    for i, node in enumerate(layout.nodes):
        ind = family.individuals.get(node.individual_id)
        adj = adjacency.get(node.individual_id, empty)
        parts.append(_svg_node(node, ind, adj, ctx, f"{prefix}n{i}"))

    return "\n".join(parts)


# SVG drawing of each carrier pattern: (rotation in degrees, lines, dots).
# Lines are "v" (vertical) and/or "h" (horizontal) before rotation.
_SVG_PATTERNS: dict[str, tuple[int, str, str]] = {
    "diagonal": (45, "v", ""),
    "antidiagonal": (-45, "v", ""),
    "vertical": (0, "v", ""),
    "horizontal": (0, "h", ""),
    "crossed": (45, "vh", ""),
    "grid": (0, "vh", ""),
    "dotted": (0, "", "dot"),
    "circles": (0, "", "ring"),
}


def _svg_pattern_id(pattern: str, obligate: bool) -> str:
    return f"carrier-{pattern}-{'sparse' if obligate else 'dense'}"


def _svg_carrier_patterns(carriers: CarrierScheme, style: PedigreeStyle) -> str:
    """``<pattern>`` definitions for every carrier source of the figure."""
    out = ["<defs>"]
    stroke = html.escape(style.fill_color)
    for key in carriers.keys:
        pattern = carriers.pattern(key)
        for obligate in (False, True):
            if obligate and key not in carriers.obligate_keys:
                continue
            size = style.shape_size * (0.25 if obligate else 0.125)
            angle, lines, dots = _SVG_PATTERNS[pattern]
            body = ""
            paint = f'stroke="{stroke}" stroke-width="0.02"'
            centre = f'cx="{size / 2:.3f}" cy="{size / 2:.3f}"'
            if "v" in lines:
                body += f'<line x1="0" y1="0" x2="0" y2="{size:.3f}" {paint}/>'
            if "h" in lines:
                body += f'<line x1="0" y1="0" x2="{size:.3f}" y2="0" {paint}/>'
            if dots == "dot":
                body += f'<circle {centre} r="{size / 8:.3f}" fill="{stroke}"/>'
            elif dots == "ring":
                body += (
                    f'<circle cx="{size / 2:.3f}" cy="{size / 2:.3f}" r="{size / 4:.3f}" '
                    f'fill="none" stroke="{stroke}" stroke-width="0.015"/>'
                )
            out.append(
                f'<pattern id="{_svg_pattern_id(pattern, obligate)}" patternUnits="userSpaceOnUse" '
                f'width="{size:.3f}" height="{size:.3f}" patternTransform="rotate({angle})">'
                f"{body}</pattern>"
            )
    out.append("</defs>")
    return "".join(out)


_PX_PER_UNIT = 60.0
_MAX_PX_W = 1600
_MAX_PX_H = 1200
_MIN_PX_W = 600


def _wrap_svg(content: str, bounds: _Extent, style: PedigreeStyle) -> str:
    """Wrap SVG content in a complete SVG document with interactivity."""
    min_x, min_y, max_x, max_y = bounds
    pad = style.margin
    vb_x = min_x - pad
    vb_y = min_y - pad
    vb_w = max_x - min_x + 2 * pad
    vb_h = max_y - min_y + 2 * pad

    # Default display size: 60 px per layout unit, scaled down to fit in
    # 1600 x 1200 and up to be at least 600 px wide (unless that would be
    # taller than 1200). One scale for both axes, so the size has the
    # drawing's proportions; clamping each axis on its own letterboxed wide
    # pedigrees with bands of white above and below.
    fit = min(_MAX_PX_W / vb_w, _MAX_PX_H / vb_h)
    scale = min(max(_PX_PER_UNIT, _MIN_PX_W / vb_w), fit)
    px_w = round(vb_w * scale)
    px_h = round(vb_h * scale)

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     viewBox="{vb_x:.2f} {vb_y:.2f} {vb_w:.2f} {vb_h:.2f}"
     width="{px_w}" height="{px_h}"
     style="background: white; font-family: {style.font_family};">
<defs>
<marker id="arrowhead" viewBox="0 0 10 10" refX="9" refY="5"
        markerWidth="5" markerHeight="5" orient="auto-start-reverse">
<path d="M0,0 L10,5 L0,10 z" fill="{style.line_color}"/></marker>
<style>
.node {{ cursor: pointer; transition: opacity 0.3s; }}
.node.dimmed {{ opacity: 0.15; }}
.node:hover .shape {{ filter: brightness(0.9); }}
.edge {{ transition: opacity 0.3s; }}
.edge.dimmed {{ opacity: 0.1; }}
.tooltip-box {{ pointer-events: none; }}
text {{ user-select: none; }}
</style>
</defs>
<g id="main-group">
{content}
</g>
<!-- Tooltip -->
<g id="tooltip" class="tooltip-box" visibility="hidden">
  <rect id="tooltip-bg" rx="0.06" ry="0.06" fill="rgba(0,0,0,0.85)" />
  <text id="tooltip-text" fill="white" font-size="0.28" />
</g>
{_INTERACTIVE_JS}
</svg>"""


# Shift from the requested anchor to the alphabetic baseline, in font sizes
# (cap height of a typical sans-serif font is about 0.7 em).
_BASELINE_SHIFT = {"hanging": 0.75, "central": 0.35, "auto": 0.0}


def _svg_text(
    x: float,
    y: float,
    text: str,
    size: float,
    style: PedigreeStyle,
    *,
    anchor: str = "middle",
    baseline: str = "central",
    color: str = "black",
    weight: str = "",
    italic: bool = False,
    cls: str = "",
) -> str:
    """One `<text>` element; `size` is in layout units.

    `baseline` says where `y` falls on the text: "hanging" (top),
    "central" (middle) or "auto" (baseline). It is turned into an explicit
    offset rather than `dominant-baseline`, which several SVG renderers
    (librsvg, older Inkscape) do not honour.
    """
    y += size * _BASELINE_SHIFT[baseline]
    extra = ""
    if weight:
        extra += f' font-weight="{weight}"'
    if italic:
        extra += ' font-style="italic"'
    if cls:
        extra += f' class="{cls}"'
    return (
        f'<text x="{x:.3f}" y="{y:.3f}" font-size="{size:.3f}" '
        f'text-anchor="{anchor}" fill="{html.escape(color)}"{extra}>{html.escape(text)}</text>'
    )


_GENDER_LABELS = {"man": "Man", "woman": "Woman", "nonbinary": "Non-binary"}
_PREGNANCY_LABELS = {
    "P": "Pregnancy",
    "SAB": "Miscarriage (SAB)",
    "TOP": "Termination of pregnancy (TOP)",
    "ECT": "Ectopic pregnancy (ECT)",
    "SB": "Stillbirth (SB)",
}
_CHILDLESS_LABELS = {"choice": "No children by choice", "infertility": "Infertility"}
_ROLE_LABELS = {"D": "Gamete donor", "S": "Surrogate"}


def _tooltip(node: NodeLayout, ind: Individual | None, carriers: CarrierScheme) -> str:
    """Text of a node's tooltip, lines separated by " | "."""
    lines = [f"ID: {node.individual_id}"]
    sex = ind.sex if ind is not None else node.sex
    lines.append(f"Sex: {({Sex.MALE: 'Male', Sex.FEMALE: 'Female'}).get(sex, 'Unknown')}")
    notation = node.notation
    if notation.gender is not None:
        lines.append(f"Gender: {_GENDER_LABELS[notation.gender.value]}")
    pheno_label = {
        Phenotype.AFFECTED: "Affected",
        Phenotype.UNAFFECTED: "Unaffected",
    }.get(node.phenotype, "Unknown")
    lines.append(f"Phenotype: {pheno_label}")

    if node.is_carrier:
        lines.append("Status: Carrier")
    if node.is_obligate_carrier:
        lines.append("Status: Obligate carrier")
    for key in carriers.keys:
        state = carrier_state(node, key)
        if key != LEGACY_CARRIER_KEY and state in (CARRIER, OBLIGATE_CARRIER):
            lines.append(carriers.label(key, state == OBLIGATE_CARRIER))

    if notation.pregnancy is not None:
        lines.append(_PREGNANCY_LABELS[notation.pregnancy.value])
    if notation.adoption is not None:
        lines.append(f"Adopted {notation.adoption.value}")
    if notation.childlessness is not None:
        lines.append(_CHILDLESS_LABELS[notation.childlessness.value])
    if notation.reproductive_role is not None:
        lines.append(_ROLE_LABELS[notation.reproductive_role.value])
    if notation.is_consultand:
        lines.append("Consultand")

    if ind is not None and ind.date_of_birth:
        lines.append(f"DOB: {ind.date_of_birth}")

    raw = node.extra_phenotypes_raw
    categorical = node.extra_phenotypes_categorical
    for name in dict.fromkeys([*raw, *categorical]):
        if name in carriers.keys:
            continue
        category = categorical.get(name)
        value = raw.get(name)
        if category is not None:
            lines.append(f"{name}: {category}")
        elif value is not None:
            lines.append(f"{name}: {'NA' if math.isnan(value) else f'{value:g}'}")
    return " | ".join(lines)


def _svg_node(
    node: NodeLayout,
    ind: Individual | None,
    adj: dict[str, list[str]],
    ctx: _SvgContext,
    uid: str,
) -> str:
    """Generate SVG for one node with data attributes for interactivity.

    Mirrors `renderer._render_layout`: a node with phenotype data is drawn
    with bands (and its PGS curve) when the figure has any, otherwise with
    its affection status and carrier sections.
    """
    style = ctx.style
    half = style.shape_size / 2
    parts = [
        f'<g class="node" data-id="{html.escape(node.individual_id)}" '
        f'data-parents="{html.escape(",".join(adj.get("parents", [])))}" '
        f'data-children="{html.escape(",".join(adj.get("children", [])))}" '
        f'data-siblings="{html.escape(",".join(adj.get("siblings", [])))}">'
    ]

    banded = ctx.banded and bool(node.extra_phenotypes) and not node.is_phantom
    triangle = symbols.is_pregnancy_loss(node.notation)
    if triangle:
        names = ctx.band_names if banded else []
        fill = pregnancy_loss_fill(node, names, style, ctx.render, ctx.display_modes)
        parts.append(_svg_pregnancy_loss(node, style, fill))
    elif banded and ctx.band_names:
        parts.append(_svg_banded_body(node, ctx, uid))
    else:
        parts.append(_svg_plain_body(node, ctx, uid))

    if not triangle:
        parts.append(_svg_legacy_carrier(node, ctx, uid))
    parts.append(_svg_notation(node, ctx, banded and bool(ctx.band_names)))

    sw = style.shape_line_width * _PT
    if node.is_deceased and not triangle:
        o = half * 0.2
        parts.append(
            f'<line x1="{node.x - half - o:.3f}" y1="{node.y + half + o:.3f}" '
            f'x2="{node.x + half + o:.3f}" y2="{node.y - half - o:.3f}" '
            f'stroke="{style.outline_color}" stroke-width="{sw:.3f}"/>'
        )

    if node.is_proband or node.notation.is_consultand:
        if triangle:
            tip_x, tip_y = symbols.triangle_points(node.x, node.y, half, up=-1.0)[1]
        else:
            tip_x, tip_y = node.x - half, node.y + half
        ax = tip_x - half * 0.5
        ay = tip_y + half * 0.5
        if node.is_proband:
            parts.append(
                _svg_text(
                    ax - half * 0.3,
                    ay + half * 0.1,
                    "P",
                    style.font_size * _PT,
                    style,
                    color=style.line_color,
                    weight="bold",
                )
            )
        parts.append(
            f'<line x1="{ax:.3f}" y1="{ay:.3f}" '
            f'x2="{tip_x:.3f}" y2="{tip_y:.3f}" '
            f'stroke="{style.line_color}" stroke-width="{sw:.3f}" '
            f'marker-end="url(#arrowhead)"/>'
        )

    if banded and ctx.pgs_names:
        parts.append(_svg_pgs(node, ctx))

    if ctx.show_labels and node.label:
        parts.append(_svg_labels(node, style))

    parts.append(f"<title>{html.escape(_tooltip(node, ind, ctx.carriers))}</title>")
    parts.append("</g>")
    return "\n".join(p for p in parts if p)


def _svg_pregnancy_loss(node: NodeLayout, style: PedigreeStyle, fill: str) -> str:
    """Triangle of a miscarriage, ectopic pregnancy or termination.

    Mirrors `renderer._draw_pregnancy_loss`: filled as decided by
    ``phenotypes.pregnancy_loss_fill``, struck through for a termination.
    """
    half = style.shape_size / 2
    points = symbols.triangle_points(node.x, node.y, half, up=-1.0)
    sw = style.shape_line_width * _PT
    out = [
        f'<polygon class="shape" points="{" ".join(f"{x:.3f},{y:.3f}" for x, y in points)}" '
        f'fill="{html.escape(fill)}" stroke="{html.escape(style.outline_color)}" '
        f'stroke-width="{sw:.3f}"/>'
    ]
    if node.notation.pregnancy == Pregnancy.TERMINATION:
        (x1, y1), (x2, y2) = symbols.termination_stroke(node.x, node.y, half, up=-1.0)
        out.append(
            f'<line class="termination" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" '
            f'y2="{y2:.3f}" stroke="{html.escape(style.outline_color)}" '
            f'stroke-width="{style.line_width * _PT:.3f}"/>'
        )
    return "\n".join(out)


def _svg_notation(node: NodeLayout, ctx: _SvgContext, banded: bool) -> str:
    """Adoption brackets and the letter inside a symbol (P, D or S).

    Mirrors `renderer._draw_notation`.
    """
    style, carriers = ctx.style, ctx.carriers
    half = style.shape_size / 2
    triangle = symbols.is_pregnancy_loss(node.notation)
    out = []
    if node.notation.adoption is not None:
        for line in symbols.bracket_lines(node.x, node.y, half, node.sex, triangle):
            pts = " ".join(f"{x:.3f},{y:.3f}" for x, y in line)
            out.append(
                f'<polyline class="adoption-bracket" points="{pts}" fill="none" '
                f'stroke="{html.escape(style.line_color)}" '
                f'stroke-width="{style.line_width * _PT:.3f}"/>'
            )
    letter = symbols.inner_letter(node.notation)
    if letter is not None and not triangle:
        divided = bool(carriers.marks(node)) or (
            carriers.legacy and (node.is_carrier or node.is_obligate_carrier)
        )
        color = style.line_color
        size = style.font_size * _PT
        if banded:
            x, y = symbols.letter_position(node.x, node.y, half, node.sex, True, up=-1.0)
            out.append(
                _svg_text(
                    x,
                    y,
                    letter,
                    size,
                    style,
                    anchor="start",
                    baseline="auto",
                    color=color,
                    cls="symbol-letter",
                )
            )
            return "\n".join(out)
        if divided:
            out.append(
                f'<rect class="symbol-letter-backing" x="{node.x - size * 0.55:.3f}" '
                f'y="{node.y - size * 0.6:.3f}" width="{size * 1.1:.3f}" '
                f'height="{size * 1.2:.3f}" rx="{size * 0.25:.3f}" fill="white" '
                f'fill-opacity="0.9"/>'
            )
        elif node.phenotype == Phenotype.AFFECTED:
            color = text_color_on(style.fill_color)
        out.append(_svg_text(node.x, node.y, letter, size, style, color=color, cls="symbol-letter"))
    return "\n".join(out)


def _svg_plain_body(node: NodeLayout, ctx: _SvgContext, uid: str) -> str:
    """Shape filled by affection status, divided into carrier sections.

    With the NSGC 2022 notation, someone with carrier results gets a symbol
    divided into vertical sections (see pypedigree.carriers).
    """
    style, carriers = ctx.style, ctx.carriers
    half = style.shape_size / 2
    legacy_carrier = carriers.legacy and (node.is_carrier or node.is_obligate_carrier)
    affected = node.phenotype == Phenotype.AFFECTED and not legacy_carrier
    marks = carriers.marks(node)
    sections: list[str] = []
    if marks:
        if affected:
            sections.append(html.escape(style.fill_color))
        sections.extend(f"url(#{_svg_pattern_id(mk.pattern, mk.obligate)})" for mk in marks)
    fill = style.fill_color if affected and not sections else "white"

    stroke = style.outline_color
    sw = style.shape_line_width * _PT
    dash = 'stroke-dasharray="0.1 0.05"' if node.is_phantom else ""
    out = [
        f'{_svg_shape_outline(node, half, cls="shape")} fill="{html.escape(fill)}" '
        f'stroke="{html.escape(stroke)}" stroke-width="{sw:.3f}" {dash}/>'
    ]
    if sections:
        out.append(_svg_sections(node, uid, half, sections, stroke, sw, dash))
    return "\n".join(out)


def _svg_banded_body(node: NodeLayout, ctx: _SvgContext, uid: str) -> str:
    """Horizontal phenotype bands, then one band per carrier source.

    Mirrors `renderer._draw_banded_node`, from the same `band_cell`
    decisions.
    """
    style, carriers = ctx.style, ctx.carriers
    half = style.shape_size / 2
    cx, cy = node.x, node.y
    keys = list(carriers.keys)
    n_bands = len(ctx.band_names) + len(keys)
    band_h = 2 * half / n_bands
    top = cy - half
    left, width = cx - half * 1.5, half * 3
    clip = f"bands-{uid}"
    clip_attr = f'clip-path="url(#{clip})"'
    sep = f'stroke="{SEPARATOR_COLOR}" stroke-width="{0.5 * _PT:.3f}" {clip_attr}'
    band_fs = style.font_size * 0.75 * _PT

    def separator(y: float) -> str:
        return f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{left + width:.3f}" y2="{y:.3f}" {sep}/>'

    out = [
        f'<defs><clipPath id="{clip}">{_svg_shape_outline(node, half)}/></clipPath></defs>',
        f'{_svg_shape_outline(node, half)} fill="white" stroke="none"/>',
    ]
    all_na = not keys and all_bands_na(node, ctx.band_names)

    for i, name in enumerate(ctx.band_names):
        cell = band_cell(node, name, style, ctx.render, ctx.display_modes)
        y = top + i * band_h
        out.append(
            f'<rect class="band" x="{left:.3f}" y="{y:.3f}" width="{width:.3f}" '
            f'height="{band_h:.3f}" fill="{to_hex(cell.color)}" {clip_attr}/>'
        )
        if i > 0:
            out.append(
                f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{left + width:.3f}" y2="{y:.3f}" {sep}/>'
            )
        if cell.text is not None and not (cell.is_na and all_na):
            out.append(
                _svg_text(
                    cx,
                    y + band_h / 2,
                    cell.text,
                    band_fs,
                    style,
                    color=cell.text_color,
                    weight="bold",
                )
            )

    for j, key in enumerate(keys):
        y = top + (len(ctx.band_names) + j) * band_h
        state = carrier_state(node, key)
        geometry = f'x="{left:.3f}" y="{y:.3f}" width="{width:.3f}" height="{band_h:.3f}"'
        out.append(f'<rect class="band" {geometry} fill="white" {clip_attr}/>')
        if state in (CARRIER, OBLIGATE_CARRIER):
            pid = _svg_pattern_id(carriers.pattern(key), state == OBLIGATE_CARRIER)
            out.append(f'<rect {geometry} fill="url(#{pid})" {clip_attr}/>')
        out.append(
            f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{left + width:.3f}" y2="{y:.3f}" {sep}/>'
        )
        if state is None:
            out.append(
                _svg_text(
                    cx, y + band_h / 2, "NA", band_fs, style, color=NA_TEXT_COLOR, weight="bold"
                )
            )

    if all_na:
        out.append(
            _svg_text(
                cx,
                cy,
                "NA",
                style.font_size * 0.85 * _PT,
                style,
                color=NA_TEXT_COLOR,
                weight="bold",
            )
        )

    border = outline_color(node, ctx.outline_names, style)
    sw = style.shape_line_width * _PT * (2.5 if border else 1.0)
    dash = 'stroke-dasharray="0.1 0.05"' if node.is_phantom else ""
    out.append(
        f'{_svg_shape_outline(node, half, cls="shape")} fill="none" '
        f'stroke="{html.escape(border or style.outline_color)}" stroke-width="{sw:.3f}" {dash}/>'
    )
    return "\n".join(out)


def _svg_legacy_carrier(node: NodeLayout, ctx: _SvgContext, uid: str) -> str:
    """Legacy notation: left half filled for a carrier, central dot for an obligate one."""
    style = ctx.style
    if not ctx.carriers.legacy:
        return ""
    half = style.shape_size / 2
    out = []
    if node.is_carrier:
        clip = f"half-{uid}"
        out.append(
            f'<defs><clipPath id="{clip}">{_svg_shape_outline(node, half)}/></clipPath></defs>'
        )
        out.append(
            f'<rect class="carrier-half-fill" x="{node.x - half * 1.5:.3f}" '
            f'y="{node.y - half * 1.5:.3f}" width="{half * 1.5:.3f}" height="{half * 3:.3f}" '
            f'fill="{html.escape(style.fill_color)}" clip-path="url(#{clip})"/>'
        )
    if node.is_obligate_carrier:
        out.append(
            f'<circle class="obligate-carrier-dot" cx="{node.x:.3f}" cy="{node.y:.3f}" '
            f'r="{style.carrier_dot_radius:.3f}" fill="{html.escape(style.fill_color)}"/>'
        )
    return "\n".join(out)


def _svg_pgs(node: NodeLayout, ctx: _SvgContext) -> str:
    """PGS distribution curve beside a symbol, one marker per score.

    Mirrors `renderer._draw_pgs_ridge`: the curve runs bottom to top, with
    each score's position on [0, 1] marked by a dashed line and a dot.
    """
    style = ctx.style
    active = active_pgs(node, ctx.pgs_names)
    if not active:
        return ""
    half = style.shape_size / 2
    box_w, box_h = half * 0.8, half * 2
    left = node.x - half - box_w - half * 0.15
    bottom = node.y + half
    thin = 0.5 * _PT

    curve = [(left + d * box_w * 0.8 + box_w * 0.1, bottom - t * box_h) for t, d in pgs_curve(50)]
    base_x = left + box_w * 0.1
    pts = " ".join(f"{x:.3f},{y:.3f}" for x, y in curve)
    area = f"{base_x:.3f},{bottom:.3f} {pts} {base_x:.3f},{bottom - box_h:.3f}"
    out = [
        '<g class="pgs">',
        f'<rect x="{left:.3f}" y="{bottom - box_h:.3f}" width="{box_w:.3f}" height="{box_h:.3f}" '
        f'rx="0.02" fill="white" fill-opacity="0.5" stroke="{SEPARATOR_COLOR}" '
        f'stroke-width="{thin:.3f}"/>',
        f'<polygon points="{area}" fill="{FALLBACK_COLOR}" fill-opacity="0.08"/>',
        f'<polyline points="{pts}" fill="none" stroke="#aaaaaa" stroke-width="{0.6 * _PT:.3f}"/>',
    ]
    for name, value in active:
        color = html.escape(style.pgs_colors.get(name, DEFAULT_PGS_COLOR))
        y = bottom - value * box_h
        out.append(
            f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{left + box_w:.3f}" y2="{y:.3f}" '
            f'stroke="{color}" stroke-width="{_PT:.3f}" stroke-dasharray="0.06 0.03"/>'
        )
        out.append(
            f'<circle cx="{left + box_w / 2:.3f}" cy="{y:.3f}" r="{1.5 * _PT:.3f}" fill="{color}"/>'
        )
    out.append(
        _svg_text(
            left + box_w / 2,
            bottom + half * 0.15,
            "PGS",
            style.font_size * 0.6 * _PT,
            style,
            baseline="hanging",
            color="#666666",
        )
    )
    out.append("</g>")
    return "\n".join(out)


# Vertical step between annotation lines under a label, in layout units.
_ANNOTATION_STEP = symbols.TEXT_LINE_STEP


def _svg_labels(node: NodeLayout, style: PedigreeStyle) -> str:
    """The ID label under a symbol, its notes ("SB", "AMAB"...), then its annotations."""
    y = node.y + style.shape_size / 2 + style.label_offset + node.marker_depth
    out = [
        _svg_text(
            node.x,
            y,
            node.label,
            style.font_size * _PT,
            style,
            baseline="hanging",
            color=style.line_color,
        )
    ]
    for text in node.notes:
        y += _ANNOTATION_STEP
        out.append(
            _svg_text(
                node.x,
                y,
                text,
                style.font_size * 0.7 * _PT,
                style,
                baseline="hanging",
                color=style.line_color,
                cls="note",
            )
        )
    for text in node.annotations:
        y += _ANNOTATION_STEP
        out.append(
            _svg_text(
                node.x,
                y,
                text,
                style.font_size * 0.7 * _PT,
                style,
                baseline="hanging",
                color="#666666",
                italic=True,
                cls="annotation",
            )
        )
    return "\n".join(out)


def _svg_shape_outline(node: NodeLayout, half: float, cls: str = "") -> str:
    """Geometry of an individual's shape, without paint attributes.

    The element is left open, so callers append attributes and ``/>``.
    """
    attr = f' class="{cls}"' if cls else ""
    if node.sex == Sex.MALE:
        return (
            f'<rect{attr} x="{node.x - half:.3f}" y="{node.y - half:.3f}" '
            f'width="{half * 2:.3f}" height="{half * 2:.3f}"'
        )
    if node.sex == Sex.FEMALE:
        return f'<circle{attr} cx="{node.x:.3f}" cy="{node.y:.3f}" r="{half:.3f}"'
    r = half * 1.2
    points = (
        f"{node.x:.3f},{node.y - r:.3f} {node.x + r:.3f},{node.y:.3f} "
        f"{node.x:.3f},{node.y + r:.3f} {node.x - r:.3f},{node.y:.3f}"
    )
    return f'<polygon{attr} points="{points}"'


def _svg_sections(
    node: NodeLayout,
    uid: str,
    half: float,
    fills: list[str],
    stroke: str,
    sw: float,
    dash: str,
) -> str:
    """Equal vertical sections clipped to the shape, then the outline on top."""
    clip_id = f"sections-{uid}"
    reach = half * 1.2 if node.sex not in (Sex.MALE, Sex.FEMALE) else half
    width = 2 * reach / len(fills)
    out = [f'<defs><clipPath id="{clip_id}">{_svg_shape_outline(node, half)}/></clipPath></defs>']
    for i, fill in enumerate(fills):
        left = node.x - reach + i * width
        out.append(
            f'<rect x="{left:.3f}" y="{node.y - reach:.3f}" width="{width:.3f}" '
            f'height="{reach * 2:.3f}" fill="{fill}" clip-path="url(#{clip_id})"/>'
        )
        if i > 0:
            out.append(
                f'<line x1="{left:.3f}" y1="{node.y - reach:.3f}" x2="{left:.3f}" '
                f'y2="{node.y + reach:.3f}" stroke="{stroke}" stroke-width="{sw / 2:.3f}" '
                f'clip-path="url(#{clip_id})"/>'
            )
    out.append(
        f'{_svg_shape_outline(node, half)} fill="none" stroke="{stroke}" '
        f'stroke-width="{sw:.3f}" {dash}/>'
    )
    return "".join(out)


# Average advance width of a character, in font sizes. Python cannot measure
# text, so the legend border is sized from this estimate.
_CHAR_WIDTH_EM = 0.55


def _svg_legend(
    rows: Sequence[LegendRow],
    style: PedigreeStyle,
    carriers: CarrierScheme,
    x0: float,
    y0: float,
) -> tuple[str, _Extent]:
    """The legend panel, its top-left corner at (x0, y0).

    Same rows (`legend_rows`) and same geometry (`legend_metrics`, in units
    of the font size) as the matplotlib legend.

    Returns:
        The SVG group and its bounding box.
    """
    binary_headers = any(r.kind == "binary" for r in rows)
    m = legend_metrics(len(rows), style, binary_headers=binary_headers)
    ux = LEG_UNIT_X_EM * style.font_size * _PT
    uy = LEG_UNIT_Y_EM * style.font_size * _PT
    fs = m.fs * _PT
    label_fs = m.label_fs * _PT
    body: list[str] = []
    right = LEG_X0 + 4.0  # rightmost content, in legend units

    def X(u: float) -> float:
        return x0 + u * ux

    def Y(v: float) -> float:
        # Legend units have y pointing up from the bottom of the panel.
        return y0 + (m.total_h - v) * uy

    def text(u: float, v: float, s: str, size: float, anchor: str = "start", **kw: Any) -> None:
        nonlocal right
        width = len(s) * _CHAR_WIDTH_EM * size / ux
        end = {"start": u + width, "middle": u + width / 2}.get(anchor, u)
        right = max(right, end)
        body.append(_svg_text(X(u), Y(v), s, size, style, anchor=anchor, **kw))

    def swatch(
        u: float,
        vc: float,
        w: float,
        h: float,
        fill: str,
        stroke: str = LEG_SWATCH_EDGE,
        width_pt: float = 0.5,
    ) -> None:
        stroke_attr = (
            f'stroke="{html.escape(stroke)}" stroke-width="{width_pt * _PT:.3f}"'
            if stroke != "none"
            else 'stroke="none"'
        )
        body.append(
            f'<rect x="{X(u):.3f}" y="{Y(vc + h / 2):.3f}" width="{w * ux:.3f}" '
            f'height="{h * uy:.3f}" fill="{html.escape(fill)}" {stroke_attr}/>'
        )

    def caption(u: float, v: float, s: str, anchor: str = "start", **kw: Any) -> None:
        kw.setdefault("color", LEG_CAPTION_COLOR)
        size = label_fs * kw.pop("size_factor", 1.0)
        text(u, v, s, size, anchor, baseline="hanging", **kw)

    # Header: title and sex symbols.
    v = m.total_h - LEG_PAD
    text(
        LEG_X0,
        v,
        "Legend",
        max(m.fs + 1, style.font_size * 0.7) * _PT,
        baseline="hanging",
        weight="bold",
    )
    v -= 1.0
    side = m.marker_size * _PT
    body.append(
        f'<rect x="{X(0.6) - side / 2:.3f}" y="{Y(v) - side / 2:.3f}" width="{side:.3f}" '
        f'height="{side:.3f}" fill="white" stroke="black" stroke-width="{_PT:.3f}"/>'
    )
    body.append(
        f'<circle cx="{X(2.5):.3f}" cy="{Y(v):.3f}" r="{side / 2:.3f}" fill="white" '
        f'stroke="black" stroke-width="{_PT:.3f}"/>'
    )
    text(0.6 + 0.25, v, "= Male", fs)
    text(2.5 + 0.25, v, "= Female", fs)
    y_top = v - 0.8

    if binary_headers:
        for i, s in enumerate(("Affected", "Unaffected", "Unknown")):
            caption(m.swatch_x(i) + m.sw_w / 2, y_top, s, "middle")
        y_top -= LEG_BINARY_HEADER_H

    for i, row in enumerate(rows):
        vc = m.row_center(i, y_top)
        if row.kind == "gradient":
            step = LEG_GRAD_W / len(row.colors)
            for k, rgb in enumerate(row.colors):
                swatch(LEG_X0 + k * step, vc, step * 1.02, m.sw_h, to_hex(rgb), stroke="none")
            swatch(LEG_X0, vc, LEG_GRAD_W, m.sw_h, "none")
            cap_v = vc - m.sw_h / 2 - 0.05
            caption(LEG_X0, cap_v, row.low)
            caption(LEG_X0 + LEG_GRAD_W, cap_v, row.high, "end")
        elif row.kind == "categorical":
            u = LEG_X0
            for value, color in row.categories.items():
                swatch(u, vc, m.sw_w * 0.8, m.sw_h, color)
                caption(
                    u + m.sw_w * 0.4,
                    vc - m.sw_h / 2 - 0.05,
                    value,
                    "middle",
                    color=LEG_CATEGORY_COLOR,
                    size_factor=0.8,
                )
                u += m.sw_w + m.sw_gap * 0.5
        elif row.kind == "binary":
            for k, fill in enumerate((row.color, "white", UNKNOWN_FILL)):
                swatch(m.swatch_x(k), vc, m.sw_w, m.sw_h, fill)
        elif row.kind == "outline":
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, "white", stroke=row.color, width_pt=2.5)
        elif row.kind == "pgs":
            baseline = vc - m.sw_h / 2
            curve = [(LEG_X0 + t * 0.9, baseline + d * m.sw_h * 0.9) for t, d in pgs_curve(30)]
            pts = " ".join(f"{X(u):.3f},{Y(w):.3f}" for u, w in curve)
            body.append(
                f'<polygon points="{X(LEG_X0):.3f},{Y(baseline):.3f} {pts} '
                f'{X(LEG_X0 + 0.9):.3f},{Y(baseline):.3f}" fill="{FALLBACK_COLOR}" '
                f'fill-opacity="0.08"/>'
            )
            body.append(
                f'<polyline points="{pts}" fill="none" stroke="#aaaaaa" '
                f'stroke-width="{0.6 * _PT:.3f}"/>'
            )
            for k, color in enumerate(row.pgs.values()):
                u = LEG_X0 + (0.3 + k * 0.2) * 0.9
                body.append(
                    f'<line x1="{X(u):.3f}" y1="{Y(baseline):.3f}" x2="{X(u):.3f}" '
                    f'y2="{Y(vc + m.sw_h / 2):.3f}" stroke="{html.escape(color)}" '
                    f'stroke-width="{0.8 * _PT:.3f}" stroke-dasharray="0.04 0.02"/>'
                )
        elif row.kind == "legacy_carrier":
            big = 9 * m.scale * _PT
            body.append(
                f'<rect x="{X(0.9) - big / 2:.3f}" y="{Y(vc) - big / 2:.3f}" width="{big:.3f}" '
                f'height="{big:.3f}" fill="white" stroke="black" stroke-width="{_PT:.3f}"/>'
            )
            if row.obligate:
                body.append(
                    f'<circle cx="{X(0.9):.3f}" cy="{Y(vc):.3f}" r="{1.5 * m.scale * _PT:.3f}" '
                    f'fill="black"/>'
                )
            else:
                body.append(
                    f'<rect x="{X(0.9) - big / 2:.3f}" y="{Y(vc) - big / 2:.3f}" '
                    f'width="{big / 2:.3f}" height="{big:.3f}" fill="black"/>'
                )
        else:  # carrier: patterned swatch
            pid = _svg_pattern_id(row.pattern, row.obligate)
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, "white", stroke=style.fill_color)
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, f"url(#{pid})", stroke="none")

        label_u = m.text_x if row.kind in ("gradient", "categorical", "binary") else m.icon_text_x
        text(label_u, vc, row.label, fs, baseline="central")

    width = (right + LEG_X0 - 0.1) * ux
    border = (
        f'<rect x="{X(0.1):.3f}" y="{Y(m.total_h - 0.1):.3f}" width="{width:.3f}" '
        f'height="{(m.total_h - 0.2) * uy:.3f}" rx="{0.1 * uy:.3f}" fill="white" '
        f'stroke="{LEG_SWATCH_EDGE}" stroke-width="{0.8 * _PT:.3f}"/>'
    )
    svg = '<g class="legend">\n' + "\n".join([border, *body]) + "\n</g>"
    return svg, (x0, y0, X(0.1) + width, y0 + m.total_h * uy)


def _svg_edge(edge: EdgeLayout, style: PedigreeStyle) -> str:
    """Generate SVG line for an edge."""
    x1, y1 = edge.start
    x2, y2 = edge.end
    sw = style.line_width * 0.02

    if edge.line_style == LineStyle.DOUBLE:
        gap = style.double_line_gap
        if y1 == y2:  # Horizontal
            return (
                f'<g class="edge">'
                f'<line x1="{x1:.3f}" y1="{y1 + gap:.3f}" '
                f'x2="{x2:.3f}" y2="{y2 + gap:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f'<line x1="{x1:.3f}" y1="{y1 - gap:.3f}" '
                f'x2="{x2:.3f}" y2="{y2 - gap:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f"</g>"
            )
        else:
            return (
                f'<g class="edge">'
                f'<line x1="{x1 + gap:.3f}" y1="{y1:.3f}" '
                f'x2="{x2 + gap:.3f}" y2="{y2:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f'<line x1="{x1 - gap:.3f}" y1="{y1:.3f}" '
                f'x2="{x2 - gap:.3f}" y2="{y2:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f"</g>"
            )

    dash = ""
    if edge.line_style == LineStyle.DASHED:
        dash = 'stroke-dasharray="0.1 0.05"'
    elif edge.line_style == LineStyle.BROKEN:
        dash = 'stroke-dasharray="0.15 0.08"'

    return (
        f'<line class="edge" x1="{x1:.3f}" y1="{y1:.3f}" '
        f'x2="{x2:.3f}" y2="{y2:.3f}" '
        f'stroke="{style.line_color}" stroke-width="{sw:.3f}" {dash}/>'
    )


def _svg_twin_connector(tc: TwinConnector, style: PedigreeStyle) -> str:
    """Generate SVG for twin connector lines."""
    cx, cy = tc.convergence_point
    sw = style.line_width * 0.02
    parts = ['<g class="edge">']

    for child_x, child_y in tc.child_positions:
        parts.append(
            f'<line x1="{cx:.3f}" y1="{cy:.3f}" '
            f'x2="{child_x:.3f}" y2="{child_y:.3f}" '
            f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
        )

    if tc.twin_type == TwinType.MONOZYGOTIC and len(tc.child_positions) >= 2:
        frac = 0.5
        bar_points = []
        for child_x, child_y in tc.child_positions:
            bx = cx + frac * (child_x - cx)
            by = cy + frac * (child_y - cy)
            bar_points.append((bx, by))
        parts.append(
            f'<line x1="{bar_points[0][0]:.3f}" y1="{bar_points[0][1]:.3f}" '
            f'x2="{bar_points[-1][0]:.3f}" y2="{bar_points[-1][1]:.3f}" '
            f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
        )

    if tc.twin_type == TwinType.UNKNOWN and len(tc.child_positions) >= 2:
        frac = 0.35
        avg_y = sum(p[1] for p in tc.child_positions) / len(tc.child_positions)
        qy = cy + frac * (avg_y - cy)
        parts.append(_svg_text(cx, qy, "?", 0.3, style, color=style.line_color, weight="bold"))

    parts.append("</g>")
    return "\n".join(parts)


def _svg_duplicate_arc(arc: DuplicateArc, style: PedigreeStyle) -> str:
    """Generate SVG for duplicate arc (dotted curve)."""
    x1, y1 = arc.from_pos
    x2, y2 = arc.to_pos
    sw = style.line_width * 0.02

    # Compute control point for quadratic bezier
    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    dx = x2 - x1
    dy = y2 - y1
    # Perpendicular offset for curvature
    ctrl_x = mid_x - dy * arc.curvature
    ctrl_y = mid_y + dx * arc.curvature

    return (
        f'<path class="edge" d="M {x1:.3f} {y1:.3f} Q {ctrl_x:.3f} {ctrl_y:.3f} '
        f'{x2:.3f} {y2:.3f}" '
        f'fill="none" stroke="{style.line_color}" stroke-width="{sw:.3f}" '
        f'stroke-dasharray="0.05 0.05"/>'
    )


_INTERACTIVE_JS = """<script>
<![CDATA[
(function() {
  const svg = document.querySelector('svg');
  const mainGroup = document.getElementById('main-group');
  const tooltip = document.getElementById('tooltip');
  const tooltipBg = document.getElementById('tooltip-bg');
  const tooltipText = document.getElementById('tooltip-text');
  const nodes = document.querySelectorAll('.node');
  const edges = document.querySelectorAll('.edge');
  let activeNode = null;

  // Click on node: highlight relatives
  nodes.forEach(node => {
    node.addEventListener('click', function(e) {
      e.stopPropagation();
      const id = this.dataset.id;
      const parents = (this.dataset.parents || '').split(',').filter(Boolean);
      const children = (this.dataset.children || '').split(',').filter(Boolean);
      const siblings = (this.dataset.siblings || '').split(',').filter(Boolean);
      const connected = new Set([id, ...parents, ...children, ...siblings]);

      if (activeNode === id) {
        // Toggle off
        nodes.forEach(n => n.classList.remove('dimmed'));
        edges.forEach(e => e.classList.remove('dimmed'));
        activeNode = null;
      } else {
        nodes.forEach(n => {
          if (connected.has(n.dataset.id)) {
            n.classList.remove('dimmed');
          } else {
            n.classList.add('dimmed');
          }
        });
        edges.forEach(e => e.classList.add('dimmed'));
        activeNode = id;
      }
    });

    // Hover: show tooltip
    node.addEventListener('mouseenter', function(e) {
      const title = this.querySelector('title');
      if (!title) return;
      const text = title.textContent;
      const lines = text.split(' | ');

      // Anchor above the node, in the SVG's own coordinates, so the tooltip
      // follows the node through zoom and pan.
      const bbox = this.getBBox();
      const toSvg = svg.getScreenCTM().inverse().multiply(this.getScreenCTM());
      const anchor = new DOMPoint(bbox.x + bbox.width / 2, bbox.y).matrixTransform(toSvg);

      tooltip.removeAttribute('transform');
      tooltipText.textContent = '';
      lines.forEach((line, i) => {
        const tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
        tspan.textContent = line;
        tspan.setAttribute('x', 0);
        tspan.setAttribute('dy', i === 0 ? '0' : '0.32');
        tooltipText.appendChild(tspan);
      });

      const tbox = tooltipText.getBBox();
      const pad = 0.1;
      tooltipBg.setAttribute('x', tbox.x - pad);
      tooltipBg.setAttribute('y', tbox.y - pad);
      tooltipBg.setAttribute('width', tbox.width + 2 * pad);
      tooltipBg.setAttribute('height', tbox.height + 2 * pad);

      // Centre over the node, kept inside the visible area.
      const vb = svg.viewBox.baseVal;
      let left = anchor.x - tbox.width / 2;
      left = Math.min(left, vb.x + vb.width - tbox.width - 2 * pad);
      left = Math.max(left, vb.x + pad);
      let top = anchor.y - 0.15 - tbox.height - pad;
      if (top < vb.y + pad) top = vb.y + pad;
      tooltip.setAttribute('transform',
        'translate(' + (left - tbox.x) + ',' + (top - tbox.y) + ')');
      tooltip.setAttribute('visibility', 'visible');
    });

    node.addEventListener('mouseleave', function() {
      tooltip.setAttribute('visibility', 'hidden');
    });
  });

  // Click on background: reset
  svg.addEventListener('click', function() {
    nodes.forEach(n => n.classList.remove('dimmed'));
    edges.forEach(e => e.classList.remove('dimmed'));
    activeNode = null;
  });

  // Zoom with mouse wheel
  let scale = 1;
  let panX = 0, panY = 0;
  let isPanning = false;
  let startX, startY;

  function updateTransform() {
    mainGroup.setAttribute('transform',
      'translate(' + panX + ',' + panY + ') scale(' + scale + ')');
  }

  svg.addEventListener('wheel', function(e) {
    e.preventDefault();
    const factor = e.deltaY > 0 ? 0.9 : 1.1;
    scale *= factor;
    scale = Math.max(0.1, Math.min(10, scale));
    updateTransform();
  }, {passive: false});

  // Pan with middle mouse button or ctrl+drag
  svg.addEventListener('mousedown', function(e) {
    if (e.button === 1 || (e.button === 0 && e.ctrlKey)) {
      isPanning = true;
      startX = e.clientX;
      startY = e.clientY;
      e.preventDefault();
    }
  });

  window.addEventListener('mousemove', function(e) {
    if (!isPanning) return;
    const vb = svg.viewBox.baseVal;
    const rect = svg.getBoundingClientRect();
    const sx = vb.width / rect.width;
    const sy = vb.height / rect.height;
    panX += (e.clientX - startX) * sx;
    panY += (e.clientY - startY) * sy;
    startX = e.clientX;
    startY = e.clientY;
    updateTransform();
  });

  window.addEventListener('mouseup', function() {
    isPanning = false;
  });
})();
]]>
</script>"""
