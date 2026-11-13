"""Data model for pedigree structures.

Pure data structures with no I/O or rendering dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class PhenotypeDisplayMode(Enum):
    """How an extra phenotype column should be rendered.

    Controlled by column name prefix in the .ped header:
    - ``b_`` → BINARY (color/white/grey bands)
    - ``g_`` → GRADIENT (white→color continuous scale)
    - ``gi_`` → GRADIENT_INVERTED (color→white: high value = white, low = color)
    - ``o_`` → OUTLINE (colored shape border)
    - ``c_`` → CATEGORICAL (distinct color per value)
    - ``cr_`` → CARRIER (carrier status of one gene or condition: 0 = not a
      carrier, 1 = carrier, 2 = obligate carrier; drawn as a patterned
      section, see ``pypedigree.carriers``)
    - no prefix → AUTO (auto-detect from data)
    """

    AUTO = "auto"
    BINARY = "binary"
    GRADIENT = "gradient"
    GRADIENT_INVERTED = "gradient_inverted"
    OUTLINE = "outline"
    CARRIER = "carrier"
    CATEGORICAL = "categorical"


class Sex(Enum):
    """Biological sex of an individual."""

    MALE = 1
    FEMALE = 2
    UNKNOWN = 0

    @classmethod
    def from_raw(cls, value: int | str) -> Sex:
        """Parse sex from raw file value.

        Args:
            value: 1=male, 2=female, 0 or -9=unknown.
        """
        try:
            v = int(value)
        except (ValueError, TypeError):
            return cls.UNKNOWN
        if v == 1:
            return cls.MALE
        if v == 2:
            return cls.FEMALE
        return cls.UNKNOWN


class Phenotype(Enum):
    """Phenotypic status of an individual."""

    UNAFFECTED = 1
    AFFECTED = 2
    UNKNOWN = 0

    @classmethod
    def from_raw(cls, value: int | str) -> Phenotype:
        """Parse phenotype from raw file value.

        Args:
            value: 1=unaffected, 2=affected, 0 or -9=unknown.
        """
        try:
            v = int(value)
        except (ValueError, TypeError):
            return cls.UNKNOWN
        if v == 1:
            return cls.UNAFFECTED
        if v == 2:
            return cls.AFFECTED
        return cls.UNKNOWN


class TwinType(Enum):
    """Zygosity type for twins."""

    MONOZYGOTIC = "MZ"
    DIZYGOTIC = "DZ"
    UNKNOWN = "TW"


class Pregnancy(Enum):
    """A pregnancy, or how one ended (Bennett et al. 2022, figure 4)."""

    ONGOING = "P"
    MISCARRIAGE = "SAB"
    TERMINATION = "TOP"
    ECTOPIC = "ECT"
    STILLBIRTH = "SB"


class Adoption(Enum):
    """Adoption status.

    ``IN``: the parents given are the adoptive parents. ``OUT``: the
    parents given are the biological parents.
    """

    IN = "in"
    OUT = "out"


class Childlessness(Enum):
    """No children by choice (or for an unknown reason), or infertility."""

    BY_CHOICE = "choice"
    INFERTILITY = "infertility"


class ReproductiveRole(Enum):
    """Role in assisted reproduction: gamete donor or surrogate."""

    DONOR = "D"
    SURROGATE = "S"


class Gender(Enum):
    """Gender identity, which decides the symbol (NSGC 2022 revision)."""

    MAN = "man"
    WOMAN = "woman"
    NONBINARY = "nonbinary"


@dataclass(frozen=True)
class Notation:
    """Standard pedigree notation beyond sex, affection status and death.

    Each attribute follows the NSGC nomenclature (Bennett et al., J Genet
    Couns 2008 and 2022); ``pypedigree.symbols`` decides how it is drawn.

    Attributes:
        pregnancy: An ongoing pregnancy ("P" inside the symbol), or a
            pregnancy loss: a triangle for a miscarriage, an ectopic
            pregnancy ("ECT" below) or a termination (struck through), and
            a struck-through symbol with "SB" below for a stillbirth.
        adoption: Brackets around the symbol; the line to adoptive parents
            is dashed.
        childlessness: A line below the symbol ending in one bar (no
            children by choice) or two (infertility).
        reproductive_role: "D" (gamete donor) or "S" (surrogate) inside the
            symbol.
        gender: Gender identity. The symbol follows it, and the sex
            assigned at birth (``Individual.sex``) is written below as AMAB
            or AFAB when the two differ.
        is_consultand: An arrow without "P" points at the symbol.
    """

    pregnancy: Pregnancy | None = None
    adoption: Adoption | None = None
    childlessness: Childlessness | None = None
    reproductive_role: ReproductiveRole | None = None
    gender: Gender | None = None
    is_consultand: bool = False


@dataclass
class TwinInfo:
    """Twin relationship information for an individual."""

    twin_type: TwinType
    sibling_id: str


@dataclass
class Individual:
    """A single person in a pedigree.

    Attributes:
        sample_id: Unique identifier within the family.
        family_id: Family this individual belongs to.
        sex: Biological sex.
        phenotype: Affected/unaffected status.
        father_id: Sample ID of the father, or None if unknown.
        mother_id: Sample ID of the mother, or None if unknown.
        date_of_birth: Date of birth string (YYYYMMDD, YYYYMM, or YYYY), or None.
        twin_info: List of twin relationships, or empty.
        is_deceased: Whether the individual is deceased.
        is_proband: Whether this is the proband.
        is_phantom: Whether this is an inferred (missing) parent.
        is_duplicate: Whether this is a duplicate node for layout purposes.
        duplicate_of: Sample ID of the original if this is a duplicate.
        notation: Pregnancy, adoption, childlessness, assisted reproduction,
            gender identity and consultand status (see ``Notation``).
    """

    sample_id: str
    family_id: str
    sex: Sex = Sex.UNKNOWN
    phenotype: Phenotype = Phenotype.UNKNOWN
    father_id: str | None = None
    mother_id: str | None = None
    date_of_birth: str | None = None
    twin_info: list[TwinInfo] = field(default_factory=list)
    is_deceased: bool = False
    is_proband: bool = False
    is_phantom: bool = False
    is_duplicate: bool = False
    is_carrier: bool = False
    is_obligate_carrier: bool = False
    duplicate_of: str | None = None
    extra_phenotypes: dict[str, float | None] = field(default_factory=dict)
    extra_phenotypes_raw: dict[str, float | None] = field(default_factory=dict)
    extra_phenotypes_categorical: dict[str, str | None] = field(default_factory=dict)
    annotations: dict[str, str] = field(default_factory=dict)
    notation: Notation = field(default_factory=Notation)

    @property
    def is_founder(self) -> bool:
        """True if this individual has no known parents."""
        return self.father_id is None and self.mother_id is None


@dataclass
class Union:
    """A mating relationship between two individuals.

    Attributes:
        partner_ids: Tuple of (partner1_id, partner2_id). By convention,
            male is first when sexes are known.
        children_ids: List of sample IDs of children from this union.
        is_consanguinous: Whether the partners share ancestors.
        is_broken: Whether this is a past/dissolved relationship.
    """

    partner_ids: tuple[str, str]
    children_ids: list[str] = field(default_factory=list)
    is_consanguinous: bool = False
    is_broken: bool = False


@dataclass
class Family:
    """A family unit containing individuals and their relationships.

    Attributes:
        family_id: Unique family identifier.
        individuals: Map of sample_id to Individual.
        unions: List of Union relationships.
    """

    family_id: str
    individuals: dict[str, Individual] = field(default_factory=dict)
    unions: list[Union] = field(default_factory=list)

    def get_founders(self) -> list[Individual]:
        """Return individuals with no known parents."""
        return [ind for ind in self.individuals.values() if ind.is_founder]

    def get_children(self, parent_id: str) -> list[Individual]:
        """Return all children of a given parent."""
        return [
            ind
            for ind in self.individuals.values()
            if ind.father_id == parent_id or ind.mother_id == parent_id
        ]

    def get_union_for_parents(self, father_id: str, mother_id: str) -> Union | None:
        """Find the union for a specific pair of parents."""
        pair = {father_id, mother_id}
        for union in self.unions:
            if set(union.partner_ids) == pair:
                return union
        return None

    def get_unions_for_individual(self, sample_id: str) -> list[Union]:
        """Return all unions involving a given individual."""
        return [u for u in self.unions if sample_id in u.partner_ids]


@dataclass
class Pedigree:
    """A collection of families parsed from a pedigree file.

    Attributes:
        families: Map of family_id to Family.
    """

    families: dict[str, Family] = field(default_factory=dict)
    phenotype_names: list[str] = field(default_factory=list)
    phenotype_display_modes: dict[str, PhenotypeDisplayMode] = field(default_factory=dict)
    annotation_names: list[str] = field(default_factory=list)

    def get_family(self, family_id: str) -> Family | None:
        """Get a family by its ID."""
        return self.families.get(family_id)

    @property
    def family_ids(self) -> list[str]:
        """List of all family IDs."""
        return list(self.families.keys())


# --- Layout result types (decoupled from rendering) ---


class EdgeType(Enum):
    """Types of edges in the pedigree diagram."""

    RELATIONSHIP = "relationship"
    DESCENT = "descent"
    SIBSHIP = "sibship"
    INDIVIDUAL = "individual"
    CHILDLESS = "childless"  # no children by choice, or infertility


class LineStyle(Enum):
    """Visual styles for edges."""

    SOLID = "solid"
    DOUBLE = "double"
    BROKEN = "broken"
    DASHED = "dashed"


@dataclass
class NodeLayout:
    """Computed position and visual attributes for one individual.

    Attributes:
        individual_id: The sample_id this node represents.
        x: Horizontal position in diagram coordinates.
        y: Vertical position in diagram coordinates.
        sex: Decides the shape: the gender identity when one is given,
            otherwise the sex.
        phenotype: Phenotype for fill style.
        is_deceased: Whether to draw deceased marker.
        is_proband: Whether to draw proband marker.
        is_phantom: Whether to draw with dashed border.
        is_duplicate: Whether this is a duplicate node.
        label: Text to display under the shape.
        notation: Notation drawn on or around the symbol (see ``Notation``).
        notes: Notation text written below the label, such as "SB" or "AMAB".
        marker_depth: Room taken below the symbol by a childlessness marker;
            the label is moved down by it.
    """

    individual_id: str
    x: float
    y: float
    sex: Sex
    phenotype: Phenotype
    is_deceased: bool = False
    is_proband: bool = False
    is_phantom: bool = False
    is_duplicate: bool = False
    is_carrier: bool = False
    is_obligate_carrier: bool = False
    label: str = ""
    extra_phenotypes: dict[str, float | None] = field(default_factory=dict)
    extra_phenotypes_raw: dict[str, float | None] = field(default_factory=dict)
    extra_phenotypes_categorical: dict[str, str | None] = field(default_factory=dict)
    rotation: float = 0.0  # Rotation angle in radians
    # Direction of the label from the node, in radians (circular layout: radially
    # outward, away from the lines arriving from the centre). None = below.
    label_angle: float | None = None
    pgs_side: str = "left"  # "left" or "right" — which side to draw PGS curve
    annotations: list[str] = field(default_factory=list)
    notation: Notation = field(default_factory=Notation)
    notes: list[str] = field(default_factory=list)
    marker_depth: float = 0.0


@dataclass
class EdgeLayout:
    """Computed position for one edge/line segment.

    Attributes:
        start: (x, y) start point.
        end: (x, y) end point.
        edge_type: What kind of connection this represents.
        line_style: How to draw the line.
    """

    start: tuple[float, float]
    end: tuple[float, float]
    edge_type: EdgeType
    line_style: LineStyle = LineStyle.SOLID
    waypoints: list[tuple[float, float]] = field(default_factory=list)


@dataclass
class TwinConnector:
    """Layout for twin connection lines.

    Twins connect from a convergence point (on the sibship line) down
    to each twin's shape with V-shaped lines.

    Attributes:
        convergence_point: (x, y) where the V-lines meet on the sibship line.
        child_positions: List of (x, y) top-center of each twin's shape.
        twin_type: Zygosity type (affects visual: MZ gets horizontal bar).
    """

    convergence_point: tuple[float, float]
    child_positions: list[tuple[float, float]]
    twin_type: TwinType


@dataclass
class DuplicateArc:
    """Layout for a dotted arc connecting duplicate instances of an individual.

    Attributes:
        from_pos: (x, y) position of the original node.
        to_pos: (x, y) position of the duplicate node.
        curvature: Signed curvature for the arc. Positive curves left
            of the from→to direction, negative curves right. The sign
            is chosen so the arc bows away from the tree centre.
    """

    from_pos: tuple[float, float]
    to_pos: tuple[float, float]
    curvature: float = 0.3


@dataclass
class FamilyLayout:
    """Complete layout result for one family.

    Contains all geometric information needed by the renderer.
    No matplotlib dependency — pure numeric data.

    Attributes:
        family_id: The family this layout represents.
        nodes: List of positioned nodes.
        edges: List of positioned edges.
        twin_connectors: List of twin connection geometries.
        duplicate_arcs: List of duplicate connection arcs.
    """

    family_id: str
    nodes: list[NodeLayout] = field(default_factory=list)
    edges: list[EdgeLayout] = field(default_factory=list)
    twin_connectors: list[TwinConnector] = field(default_factory=list)
    duplicate_arcs: list[DuplicateArc] = field(default_factory=list)
