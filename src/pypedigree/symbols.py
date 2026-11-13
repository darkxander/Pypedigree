"""How standard pedigree notation is drawn, shared by both renderers.

Decisions and geometry for the notation of ``models.Notation``: pregnancies
and pregnancy losses, adoption, childlessness, assisted reproduction, gender
identity and the consultand. Both the PNG renderer (``renderer.py``) and the
interactive SVG (``svg_renderer.py``) draw from these functions, so the two
outputs cannot drift apart. No matplotlib here.

Geometry is in layout units with the y axis pointing up (as in
``FamilyLayout``); functions whose result is not symmetric top to bottom
take an ``up`` sign, so the SVG renderer, whose y axis points down, passes
``up=-1``.

Reference: Bennett RL, French KS, Resta RG, Austin J. Standardized pedigree
nomenclature update centered on sex and gender inclusivity. J Genet Couns
2022;31(6):1238-48. Figures 1-5. The 2022 revision kept the 2008 symbols for
pregnancies and assisted reproduction (Bennett et al., J Genet Couns
2008;17:424-33).
"""

from __future__ import annotations

import math

from pypedigree.models import (
    Childlessness,
    Gender,
    Individual,
    Notation,
    Pregnancy,
    Sex,
)

# The unknown-sex diamond's vertices sit this far from its centre, relative
# to half the shape size, so it reads as the same size as a square.
DIAMOND_SCALE = 1.2

# A pregnancy loss is a triangle with its apex where the line from the
# parents arrives (the top of an ordinary symbol). Its base lies this far
# below the centre and it is this wide, both relative to half the shape
# size: an equilateral triangle a little shorter than a square.
TRIANGLE_BASE = 0.7
TRIANGLE_HALF_WIDTH = 0.98

# Childlessness marker: a stem below the symbol ending in a bar, two bars
# for infertility. Lengths in layout units, except the bar's half-width,
# relative to half the shape size.
CHILDLESS_STEM = 0.25
CHILDLESS_BAR_GAP = 0.07
CHILDLESS_BAR_HALF_WIDTH = 0.75
# Gap between the marker and the label below it.
CHILDLESS_LABEL_GAP = 0.05

# Adoption brackets: gap to the symbol, length of the ticks pointing
# inward, and how far they reach above and below the symbol.
BRACKET_GAP = 0.08
BRACKET_TICK = 0.1
BRACKET_OVERHANG = 0.08

# Vertical step between lines of text below the label (notes, annotations).
TEXT_LINE_STEP = 0.25

_GENDER_SHAPE = {Gender.MAN: Sex.MALE, Gender.WOMAN: Sex.FEMALE, Gender.NONBINARY: Sex.UNKNOWN}
_ASSIGNED_AT_BIRTH = {Sex.MALE: "AMAB", Sex.FEMALE: "AFAB"}
_PREGNANCY_LOSSES = frozenset({Pregnancy.MISCARRIAGE, Pregnancy.TERMINATION, Pregnancy.ECTOPIC})


def symbol_sex(sex: Sex, notation: Notation) -> Sex:
    """The sex whose shape is drawn: square, circle or diamond.

    The 2022 revision draws gender identity, not the sex assigned at birth:
    a square for a man, a circle for a woman, a diamond for a non-binary or
    gender-diverse person.
    """
    if notation.gender is None:
        return sex
    return _GENDER_SHAPE[notation.gender]


def sex_at_birth_note(sex: Sex, notation: Notation) -> str | None:
    """ "AMAB" or "AFAB" when the symbol does not show the sex at birth.

    Nothing is written for a cisgender person, nor when the sex assigned at
    birth is unknown.
    """
    if notation.gender is None or sex not in _ASSIGNED_AT_BIRTH:
        return None
    if _GENDER_SHAPE[notation.gender] == sex:
        return None
    return _ASSIGNED_AT_BIRTH[sex]


def notes(ind: Individual) -> list[str]:
    """Notation text written below the label, top to bottom.

    "SB" for a stillbirth and "ECT" for an ectopic pregnancy, then the sex
    assigned at birth when the symbol does not show it.
    """
    out = []
    if ind.notation.pregnancy in (Pregnancy.STILLBIRTH, Pregnancy.ECTOPIC):
        out.append(ind.notation.pregnancy.value)
    note = sex_at_birth_note(ind.sex, ind.notation)
    if note:
        out.append(note)
    return out


def is_pregnancy_loss(notation: Notation) -> bool:
    """True for a miscarriage, an ectopic pregnancy or a termination: a triangle."""
    return notation.pregnancy in _PREGNANCY_LOSSES


def is_struck_through(ind: Individual) -> bool:
    """Whether the symbol gets the deceased line; a stillbirth always does."""
    return ind.is_deceased or ind.notation.pregnancy == Pregnancy.STILLBIRTH


def inner_letter(notation: Notation) -> str | None:
    """Letter written inside the symbol: "P" (pregnancy), "D" or "S"."""
    if notation.pregnancy == Pregnancy.ONGOING:
        return "P"
    if notation.reproductive_role is not None:
        return notation.reproductive_role.value
    return None


def letter_position(
    x: float, y: float, half: float, sex: Sex, outside: bool, up: float = 1.0
) -> tuple[float, float]:
    """Where the inner letter goes: the centre, or beside the upper right.

    A symbol divided into phenotype bands has text of its own in the
    middle, so the letter moves outside, next to its upper-right corner;
    the renderer anchors it there by its lower-left corner.
    """
    if not outside:
        return (x, y)
    return (x + shape_half_width(sex, half) + LETTER_OFFSET, y + up * half * 0.55)


# Gap between a symbol and a letter written beside it.
LETTER_OFFSET = 0.04


def triangle_points(x: float, y: float, half: float, up: float = 1.0) -> list[tuple[float, float]]:
    """Vertices of a pregnancy-loss triangle: apex, base left, base right."""
    base_y = y - up * half * TRIANGLE_BASE
    w = half * TRIANGLE_HALF_WIDTH
    return [(x, y + up * half), (x - w, base_y), (x + w, base_y)]


def termination_stroke(
    x: float, y: float, half: float, up: float = 1.0
) -> tuple[tuple[float, float], tuple[float, float]]:
    """The line struck through a termination's triangle, lower left to upper right."""
    w = half * TRIANGLE_HALF_WIDTH
    return (x - w * 1.1, y - up * half * 0.95), (x + w * 0.8, y + up * half * 0.75)


def shape_half_width(sex: Sex, half: float) -> float:
    """Horizontal half-extent of a shape (the diamond reaches further)."""
    return half * DIAMOND_SCALE if sex not in (Sex.MALE, Sex.FEMALE) else half


def shape_reach(sex: Sex, half: float, direction: float, rotation: float = 0.0) -> float:
    """Distance from a shape's centre to its outline in a given direction.

    Args:
        sex: Decides the shape: square (male), circle (female) or diamond.
        half: Half the shape size.
        direction: Direction to measure along, in radians.
        rotation: Rotation applied to the shape, in radians.
    """
    theta = direction - rotation  # the direction in the shape's own frame
    c, s = abs(math.cos(theta)), abs(math.sin(theta))
    if sex == Sex.MALE:
        return half / max(c, s)
    if sex == Sex.FEMALE:
        return half
    # The diamond is a square standing on a vertex: its outline is |x|+|y| = R.
    return half * DIAMOND_SCALE / (c + s)


def bracket_lines(
    x: float, y: float, half: float, sex: Sex, triangle: bool = False
) -> list[list[tuple[float, float]]]:
    """The two adoption brackets around a symbol, as polylines.

    Symmetric top to bottom, so the same points serve both y directions.
    """
    if triangle:
        reach_x, reach_y = half * TRIANGLE_HALF_WIDTH, half
    else:
        # Square and circle reach `half` both ways; the diamond further.
        reach_x = reach_y = shape_half_width(sex, half)
    left = x - reach_x - BRACKET_GAP
    right = x + reach_x + BRACKET_GAP
    top = y + reach_y + BRACKET_OVERHANG
    bottom = y - reach_y - BRACKET_OVERHANG
    lt, rt = left + BRACKET_TICK, right - BRACKET_TICK
    return [
        [(lt, top), (left, top), (left, bottom), (lt, bottom)],
        [(rt, top), (right, top), (right, bottom), (rt, bottom)],
    ]


def childless_depth(notation: Notation) -> float:
    """Room the childlessness marker takes below the symbol, label gap included."""
    if notation.childlessness is None:
        return 0.0
    depth = CHILDLESS_STEM + CHILDLESS_LABEL_GAP
    if notation.childlessness == Childlessness.INFERTILITY:
        depth += CHILDLESS_BAR_GAP
    return depth


def childless_segments(
    x: float,
    y: float,
    sex: Sex,
    half: float,
    childlessness: Childlessness,
    direction: float = -math.pi / 2,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Segments of a childlessness marker: the stem, then one or two bars.

    The stem leaves the symbol in ``direction`` (straight down by default;
    outward in the circular layout) and the bars cross its end.
    """
    c, s = math.cos(direction), math.sin(direction)
    start = shape_reach(sex, half, direction)
    end = start + CHILDLESS_STEM
    w = half * CHILDLESS_BAR_HALF_WIDTH

    def at(d: float, across: float) -> tuple[float, float]:
        # `across` runs perpendicular to the stem.
        return (x + d * c - across * s, y + d * s + across * c)

    segments = [(at(start, 0.0), at(end, 0.0)), (at(end, -w), at(end, w))]
    if childlessness == Childlessness.INFERTILITY:
        d = end + CHILDLESS_BAR_GAP
        segments.append((at(d, -w), at(d, w)))
    return segments


def text_lines_below(ind: Individual) -> int:
    """Lines of text below the label: notes, then annotations."""
    return len(notes(ind)) + sum(1 for v in ind.annotations.values() if v)
