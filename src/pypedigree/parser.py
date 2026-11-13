"""Pedigree file parser.

Reads pedigree files in various formats (TSV, CSV, space-separated),
validates the data, infers unions from parent references, and builds
the data model.
"""

from __future__ import annotations

import warnings
from collections import deque
from dataclasses import dataclass, field, replace
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    import pandas as pd

from pypedigree.errors import (
    PedigreeDataWarning,
    PedigreeParseError,
    PedigreeValidationError,
)
from pypedigree.models import (
    Adoption,
    Childlessness,
    Family,
    Gender,
    Individual,
    Notation,
    Pedigree,
    Phenotype,
    PhenotypeDisplayMode,
    Pregnancy,
    ReproductiveRole,
    Sex,
    TwinInfo,
    TwinType,
    Union,
)

# Known header tokens (case-insensitive) for auto-detection
_HEADER_TOKENS = frozenset(
    {
        "fid",
        "iid",
        "family_id",
        "sample_id",
        "paternal_id",
        "maternal_id",
        "father",
        "mother",
        "sex",
        "phenotype",
        "pat",
        "mat",
        "pid",
        "mid",
    }
)

# Column name aliases → canonical role (case-insensitive)
_COLUMN_ROLES: dict[str, str] = {
    "fid": "family_id",
    "family_id": "family_id",
    "family": "family_id",
    "iid": "sample_id",
    "id": "sample_id",
    "sample_id": "sample_id",
    "individual_id": "sample_id",
    "ind": "sample_id",
    "paternal_id": "father",
    "father": "father",
    "pat": "father",
    "pid": "father",
    "father_id": "father",
    "dad": "father",
    "maternal_id": "mother",
    "mother": "mother",
    "mat": "mother",
    "mid": "mother",
    "mother_id": "mother",
    "mom": "mother",
    "sex": "sex",
    "gender": "sex",
    "phenotype": "phenotype",
    "pheno": "phenotype",
    "affected": "phenotype",
    "status": "phenotype",
    "date_of_birth": "dob",
    "dob": "dob",
    "birth": "dob",
    "twins": "twins",
    "twin": "twins",
    "deceased": "deceased",
    "dead": "deceased",
    "proband": "proband",
    "carrier": "carrier",
    # Standard notation (see models.Notation). "gender" means the sex
    # column, unless the header also has one (see _parse_header_columns).
    "pregnancy": "pregnancy",
    "pregnancy_outcome": "pregnancy",
    "outcome": "pregnancy",
    "adopted": "adopted",
    "adoption": "adopted",
    "childless": "childless",
    "childlessness": "childless",
    "infertility": "childless",
    "art": "art",
    "assisted_reproduction": "art",
    "gender_identity": "gender_identity",
    "identity": "gender_identity",
    "consultand": "consultand",
}

# Columns read into models.Notation.
_NOTATION_ROLES = ("pregnancy", "adopted", "childless", "art", "gender_identity", "consultand")

# Standard structural columns (order matters for the fixed-position fallback)
_STRUCTURAL_ROLES = {
    "family_id",
    "sample_id",
    "father",
    "mother",
    "sex",
    "phenotype",
    "dob",
    "twins",
    "deceased",
    "proband",
    "carrier",
    *_NOTATION_ROLES,
}


@dataclass
class ColumnMapping:
    """Maps column indices to their roles based on header analysis."""

    family_id: int = 0
    sample_id: int = 1
    father: int = 2
    mother: int = 3
    sex: int = 4
    phenotype: int | None = 5
    dob: int | None = None
    twins: int | None = None
    deceased: int | None = None
    proband: int | None = None
    carrier: int | None = None
    extra_phenotype_cols: list[tuple[int, str]] = field(default_factory=list)
    """List of (column_index, phenotype_name) for extra phenotype columns."""
    annotation_cols: list[tuple[int, str]] = field(default_factory=list)
    """List of (column_index, annotation_name) for a_ prefixed annotation columns."""
    notation_cols: dict[str, int] = field(default_factory=dict)
    """Map of notation role (pregnancy, adopted, ...) → column index."""
    phenotype_display_modes: dict[str, PhenotypeDisplayMode] = field(default_factory=dict)
    """Map of phenotype_name → display mode detected from prefix."""

    @property
    def is_header_mapped(self) -> bool:
        """True if this mapping was built from an actual header row."""
        return len(self.extra_phenotype_cols) > 0 or self.phenotype is None

    @property
    def is_positional(self) -> bool:
        """True if every structural column is where the positional parser reads it.

        The positional parser reads columns 1-10 in the standard order and a
        legacy ``carrier`` column anywhere from column 11. A header that
        names a structural column elsewhere (``Twins`` as column 7, a
        ``Deceased`` column with no ``Twins`` before it) must be read by name.
        """
        for role, position in _POSITIONAL_COLUMNS.items():
            index = getattr(self, role)
            if index is not None and index != position:
                return False
        return self.carrier is None or self.carrier >= len(_POSITIONAL_COLUMNS)


# Column index each structural role has in the headerless (positional) format.
_POSITIONAL_COLUMNS = {
    "family_id": 0,
    "sample_id": 1,
    "father": 2,
    "mother": 3,
    "sex": 4,
    "phenotype": 5,
    "dob": 6,
    "twins": 7,
    "deceased": 8,
    "proband": 9,
}


_DISPLAY_MODE_PREFIXES: dict[str, PhenotypeDisplayMode] = {
    "cr_": PhenotypeDisplayMode.CARRIER,
    "b_": PhenotypeDisplayMode.BINARY,
    "gi_": PhenotypeDisplayMode.GRADIENT_INVERTED,
    "g_": PhenotypeDisplayMode.GRADIENT,
    "o_": PhenotypeDisplayMode.OUTLINE,
    "c_": PhenotypeDisplayMode.CATEGORICAL,
}


def _strip_display_prefix(
    name: str,
) -> tuple[str, PhenotypeDisplayMode]:
    """Strip display mode prefix from a phenotype column name.

    Args:
        name: Raw column name (e.g. "b_ASD", "g_IQ", "ASD").

    Returns:
        Tuple of (stripped_name, display_mode).
    """
    name_lower = name.lower()
    for prefix, mode in _DISPLAY_MODE_PREFIXES.items():
        if name_lower.startswith(prefix) and len(name) > len(prefix):
            return name[len(prefix) :], mode
    return name, PhenotypeDisplayMode.AUTO


def _parse_header_columns(header_fields: list[str]) -> ColumnMapping:
    """Parse header column names and build a ColumnMapping.

    Recognises known structural column names (FID, ID, Father, Mother, Sex,
    Phenotype, DOB, Twins, Deceased, Proband, Carrier) by alias lookup.
    Any unrecognised column name after the core 5 (FID..Sex) is treated as
    an extra phenotype column.

    Args:
        header_fields: List of column name strings from the header row.

    Returns:
        A ColumnMapping describing which column index maps to which role.
    """
    mapping = ColumnMapping()
    role_assigned: dict[str, int] = {}
    extra_phenos: list[tuple[int, str]] = []
    annotation_cols: list[tuple[int, str]] = []
    display_modes: dict[str, PhenotypeDisplayMode] = {}
    # With both a "sex" and a "gender" column, "gender" is the gender
    # identity; alone, it is the sex column, as it always was.
    lowered = {f.strip().lower() for f in header_fields}
    roles = dict(_COLUMN_ROLES)
    if "gender" in lowered and lowered & {"sex"}:
        roles["gender"] = "gender_identity"

    for i, raw_name in enumerate(header_fields):
        name = raw_name.strip()
        name_lower = name.lower()
        # Check for annotation prefix (a_) before anything else
        if name_lower.startswith("a_") and len(name) > 2 and i >= 5:
            annotation_cols.append((i, name[2:]))  # Strip a_ prefix
            continue
        # Check if this is a known structural column (using stripped name too)
        stripped_name, mode = _strip_display_prefix(name)
        role = roles.get(name_lower)

        if role and role not in role_assigned:
            role_assigned[role] = i
        elif role and role in role_assigned:
            # Duplicate structural column → treat as phenotype
            extra_phenos.append((i, stripped_name))
            if mode != PhenotypeDisplayMode.AUTO:
                display_modes[stripped_name] = mode
        else:
            # Not a known structural column → extra phenotype
            if i >= 5:  # Only cols after the core 5 can be phenotypes
                extra_phenos.append((i, stripped_name))
                if mode != PhenotypeDisplayMode.AUTO:
                    display_modes[stripped_name] = mode

    # Assign structural columns from role_assigned
    if "family_id" in role_assigned:
        mapping.family_id = role_assigned["family_id"]
    if "sample_id" in role_assigned:
        mapping.sample_id = role_assigned["sample_id"]
    if "father" in role_assigned:
        mapping.father = role_assigned["father"]
    if "mother" in role_assigned:
        mapping.mother = role_assigned["mother"]
    if "sex" in role_assigned:
        mapping.sex = role_assigned["sex"]
    mapping.phenotype = role_assigned.get("phenotype")
    if "dob" in role_assigned:
        mapping.dob = role_assigned["dob"]
    if "twins" in role_assigned:
        mapping.twins = role_assigned["twins"]
    if "deceased" in role_assigned:
        mapping.deceased = role_assigned["deceased"]
    if "proband" in role_assigned:
        mapping.proband = role_assigned["proband"]
    if "carrier" in role_assigned:
        mapping.carrier = role_assigned["carrier"]
    mapping.notation_cols = {r: role_assigned[r] for r in _NOTATION_ROLES if r in role_assigned}

    mapping.extra_phenotype_cols = extra_phenos
    mapping.annotation_cols = annotation_cols
    mapping.phenotype_display_modes = display_modes
    return mapping


_MISSING_VALUES = frozenset({"0", "-9", "nan", "na", ".", ""})

# For extra phenotype columns (11+), "0" is a valid value (binary = unaffected)
_PHENO_MISSING_VALUES = frozenset({"-9", "nan", "na", ".", ""})


def _is_missing(value: str) -> bool:
    """Check if a value represents missing data."""
    return value.strip().lower() in _MISSING_VALUES


def _is_pheno_missing(value: str) -> bool:
    """Check if a phenotype value represents missing data (0 is valid)."""
    return value.strip().lower() in _PHENO_MISSING_VALUES


def _detect_separator(line: str) -> str | None:
    """Auto-detect the field separator from a line of text.

    Tries tab first, then comma, then whitespace.

    Args:
        line: A single line from the pedigree file.

    Returns:
        The detected separator string, or None to signal splitting on any
        whitespace (what `str.split()` with no argument does).

    Raises:
        PedigreeParseError: If no separator can produce >= 6 fields.
    """
    for sep in ["\t", ","]:
        parts = line.split(sep)
        if len(parts) >= 6:
            return sep

    # Try whitespace (split with no args splits on any whitespace)
    parts = line.split()
    if len(parts) >= 6:
        return None  # None signals whitespace splitting

    raise PedigreeParseError(f"Cannot detect separator: line has fewer than 6 fields: {line!r}")


def _detect_header(first_line: str, separator: str | None) -> bool:
    """Detect whether the first line is a header row.

    A line is considered a header if its first field (lowercased) matches
    a known header token.

    Args:
        first_line: The first line of the file.
        separator: The detected separator (None for whitespace).

    Returns:
        True if the first line appears to be a header.
    """
    parts = first_line.split(separator) if separator else first_line.split()
    if not parts:
        return False
    return parts[0].strip().lower() in _HEADER_TOKENS


def _split_line(line: str, separator: str | None) -> list[str]:
    """Split a line using the detected separator."""
    if separator is None:
        return line.split()
    return line.split(separator)


def _parse_individual(
    fields: list[str],
    line_num: int,
    phenotype_names: list[str] | None = None,
    carrier_col_index: int | None = None,
) -> Individual:
    """Parse one row of a pedigree file into an Individual.

    Handles the standard 6-column format and the extended format (up to 10 columns):
    - Columns 1-6: family_id, sample_id, paternal_id, maternal_id, sex, phenotype
    - Column 7 (optional): date_of_birth
    - Column 8 (optional): twins (raw string, parsed in later phases)
    - Column 9 (optional): deceased (1=deceased, 0=alive)
    - Column 10 (optional): proband (1=proband, 0=not)

    Args:
        fields: List of field values (at least 6 elements).
        line_num: Line number in the file (for error messages).

    Returns:
        An Individual instance.

    Raises:
        PedigreeParseError: If the row has fewer than 6 fields.
    """
    if len(fields) < 6:
        raise PedigreeParseError(f"Line {line_num}: expected at least 6 fields, got {len(fields)}")

    family_id = fields[0].strip()
    sample_id = fields[1].strip()
    father_raw = fields[2].strip()
    mother_raw = fields[3].strip()
    sex_raw = fields[4].strip()
    phenotype_raw = fields[5].strip()

    father_id = None if _is_missing(father_raw) else father_raw
    mother_id = None if _is_missing(mother_raw) else mother_raw

    # Extended columns (all optional)
    date_of_birth = None
    is_deceased = False
    is_proband = False

    # Column 7: date_of_birth
    if len(fields) > 6:
        dob_raw = fields[6].strip()
        if not _is_missing(dob_raw):
            date_of_birth = dob_raw

    # Column 8: twins — comma-separated TYPE:sample_id entries
    twin_info: list[TwinInfo] = []
    if len(fields) > 7:
        twins_raw = fields[7].strip()
        if not _is_missing(twins_raw):
            twin_info = _parse_twin_field(twins_raw, line_num)

    # Column 9: deceased
    if len(fields) > 8:
        deceased_raw = fields[8].strip()
        if not _is_missing(deceased_raw):
            is_deceased = deceased_raw == "1"

    # Column 10: proband
    if len(fields) > 9:
        proband_raw = fields[9].strip()
        if not _is_missing(proband_raw):
            is_proband = proband_raw == "1"

    # Columns 11+: extra phenotypes (only if header provided phenotype names)
    # The phenotype_names list excludes the carrier column, so we iterate
    # over all columns from index 10+ and skip the carrier column index.
    extra_phenotypes_raw: dict[str, float | None] = {}
    is_carrier = False
    is_obligate_carrier = False
    if len(fields) > 10:
        pheno_idx = 0
        for col_idx in range(10, len(fields)):
            if col_idx == carrier_col_index:
                # Parse carrier column: 0=not, 1=carrier, 2=obligate carrier
                val_raw = fields[col_idx].strip()
                if not _is_missing(val_raw):
                    try:
                        val = int(val_raw)
                    except ValueError:
                        val = 0
                    is_carrier = val == 1
                    is_obligate_carrier = val == 2
            elif phenotype_names and pheno_idx < len(phenotype_names):
                name = phenotype_names[pheno_idx]
                pheno_idx += 1
                val_raw = fields[col_idx].strip()
                if _is_pheno_missing(val_raw):
                    extra_phenotypes_raw[name] = None
                else:
                    try:
                        extra_phenotypes_raw[name] = float(val_raw)
                    except ValueError:
                        extra_phenotypes_raw[name] = None

    return Individual(
        sample_id=sample_id,
        family_id=family_id,
        sex=Sex.from_raw(sex_raw),
        phenotype=Phenotype.from_raw(phenotype_raw),
        father_id=father_id,
        mother_id=mother_id,
        date_of_birth=date_of_birth,
        twin_info=twin_info,
        is_deceased=is_deceased,
        is_proband=is_proband,
        is_carrier=is_carrier,
        is_obligate_carrier=is_obligate_carrier,
        extra_phenotypes_raw=extra_phenotypes_raw,
    )


def _parse_individual_mapped(
    fields: list[str],
    line_num: int,
    col_map: ColumnMapping,
) -> Individual:
    """Parse one row using a header-derived ColumnMapping.

    Unlike _parse_individual which assumes fixed column positions,
    this function uses the ColumnMapping to read columns by name.
    """

    def _get(idx: int | None) -> str:
        if idx is not None and idx < len(fields):
            return fields[idx].strip()
        return ""

    family_id = _get(col_map.family_id)
    sample_id = _get(col_map.sample_id)
    father_raw = _get(col_map.father)
    mother_raw = _get(col_map.mother)
    sex_raw = _get(col_map.sex)

    if col_map.phenotype is not None:
        phenotype_raw = _get(col_map.phenotype)
    elif col_map.extra_phenotype_cols:
        # No explicit phenotype column → use first extra phenotype as primary
        first_pheno_idx = col_map.extra_phenotype_cols[0][0]
        phenotype_raw = _get(first_pheno_idx)
    else:
        phenotype_raw = "0"

    father_id = None if _is_missing(father_raw) else father_raw
    mother_id = None if _is_missing(mother_raw) else mother_raw

    # Optional structural columns
    date_of_birth = None
    dob_raw = _get(col_map.dob)
    if dob_raw and not _is_missing(dob_raw):
        date_of_birth = dob_raw

    twin_info: list[TwinInfo] = []
    twins_raw = _get(col_map.twins)
    if twins_raw and not _is_missing(twins_raw):
        twin_info = _parse_twin_field(twins_raw, line_num)

    is_deceased = False
    deceased_raw = _get(col_map.deceased)
    if deceased_raw and not _is_missing(deceased_raw):
        is_deceased = deceased_raw == "1"

    is_proband = False
    proband_raw = _get(col_map.proband)
    if proband_raw and not _is_missing(proband_raw):
        is_proband = proband_raw == "1"

    is_carrier = False
    is_obligate_carrier = False
    carrier_raw = _get(col_map.carrier)
    if carrier_raw and not _is_missing(carrier_raw):
        try:
            val = int(carrier_raw)
        except ValueError:
            val = 0
        is_carrier = val == 1
        is_obligate_carrier = val == 2

    # Extra phenotypes from named columns
    extra_phenotypes_raw: dict[str, float | None] = {}
    extra_phenotypes_categorical: dict[str, str | None] = {}
    for col_idx, name in col_map.extra_phenotype_cols:
        val_raw = _get(col_idx)
        mode = col_map.phenotype_display_modes.get(name, PhenotypeDisplayMode.AUTO)
        if not val_raw or _is_pheno_missing(val_raw):
            extra_phenotypes_raw[name] = None
            if mode == PhenotypeDisplayMode.CATEGORICAL:
                extra_phenotypes_categorical[name] = None
        elif mode == PhenotypeDisplayMode.CATEGORICAL:
            # Categorical: store string value, also try numeric for raw
            extra_phenotypes_categorical[name] = val_raw.strip()
            try:
                extra_phenotypes_raw[name] = float(val_raw.replace(",", "."))
            except ValueError:
                extra_phenotypes_raw[name] = None
        else:
            try:
                extra_phenotypes_raw[name] = float(val_raw.replace(",", "."))
            except ValueError:
                extra_phenotypes_raw[name] = None

    # Annotations from a_ prefixed columns
    annotations: dict[str, str] = {}
    _ANN_MISSING = frozenset({"0", "-9", "nan", "na", ".", ""})
    for col_idx, name in col_map.annotation_cols:
        val_raw = _get(col_idx)
        if val_raw and val_raw.strip().lower() not in _ANN_MISSING:
            annotations[name] = val_raw.strip()

    if not sample_id:
        raise PedigreeParseError(f"Line {line_num}: missing sample_id")

    notation = _parse_notation({r: _get(i) for r, i in col_map.notation_cols.items()}, line_num)

    return Individual(
        sample_id=sample_id,
        family_id=family_id,
        sex=Sex.from_raw(sex_raw),
        phenotype=Phenotype.from_raw(phenotype_raw),
        father_id=father_id,
        mother_id=mother_id,
        date_of_birth=date_of_birth,
        twin_info=twin_info,
        is_deceased=is_deceased,
        is_proband=is_proband,
        is_carrier=is_carrier,
        is_obligate_carrier=is_obligate_carrier,
        extra_phenotypes_raw=extra_phenotypes_raw,
        extra_phenotypes_categorical=extra_phenotypes_categorical,
        annotations=annotations,
        notation=notation,
    )


_E = TypeVar("_E", bound=Enum)


def _codes(entries: dict[_E, tuple[str, ...]]) -> dict[str, _E]:
    """Lower-case code → enum value, from each value's accepted codes."""
    return {code: value for value, codes in entries.items() for code in codes}


_PREGNANCY_CODES = _codes(
    {
        Pregnancy.ONGOING: ("p", "pregnancy", "pregnant", "ongoing"),
        Pregnancy.MISCARRIAGE: ("sab", "miscarriage", "spontaneous_abortion"),
        Pregnancy.TERMINATION: ("top", "termination"),
        Pregnancy.ECTOPIC: ("ect", "ectopic"),
        Pregnancy.STILLBIRTH: ("sb", "stillbirth", "stillborn"),
    }
)
_ADOPTION_CODES = _codes(
    {
        Adoption.IN: ("1", "in", "adopted_in"),
        Adoption.OUT: ("2", "out", "adopted_out"),
    }
)
_CHILDLESS_CODES = _codes(
    {
        Childlessness.BY_CHOICE: ("1", "choice", "by_choice", "yes"),
        Childlessness.INFERTILITY: ("2", "infertility", "infertile"),
    }
)
_ART_CODES = _codes(
    {
        ReproductiveRole.DONOR: ("d", "donor"),
        ReproductiveRole.SURROGATE: ("s", "surrogate", "gc", "gestational_carrier"),
    }
)
_GENDER_CODES = _codes(
    {
        Gender.MAN: ("1", "m", "man", "male", "boy"),
        Gender.WOMAN: ("2", "f", "w", "woman", "female", "girl"),
        Gender.NONBINARY: (
            "3",
            "nb",
            "x",
            "nonbinary",
            "non-binary",
            "gender-diverse",
            "gender_diverse",
            "diverse",
        ),
    }
)


def _notation_value(
    raw: dict[str, str], role: str, codes: dict[str, _E], line_num: int
) -> _E | None:
    """Value of one notation column, or None when the cell is missing.

    Raises:
        PedigreeParseError: If the cell is not one of the column's codes.
    """
    text = raw.get(role, "")
    if not text or _is_missing(text):
        return None
    value = codes.get(text.strip().lower())
    if value is None:
        expected = ", ".join(sorted(codes))
        raise PedigreeParseError(
            f"Line {line_num}: unknown {role} value {text!r}; expected one of {expected}"
        )
    return value


def _coerce_notation(value: _E | str | None, role: str, codes: dict[str, _E]) -> _E | None:
    """An enum value given as such, or as one of the column's codes."""
    if value is None or isinstance(value, Enum):
        return value
    text = str(value)
    if _is_missing(text):
        return None
    found = codes.get(text.strip().lower())
    if found is None:
        raise PedigreeParseError(
            f"Unknown {role} value {value!r}; expected one of {', '.join(sorted(codes))}"
        )
    return found


def build_notation(
    *,
    pregnancy: Pregnancy | str | None = None,
    adopted: Adoption | str | None = None,
    childless: Childlessness | str | None = None,
    art: ReproductiveRole | str | None = None,
    gender: Gender | str | None = None,
    consultand: bool = False,
) -> Notation:
    """Build a Notation from enum values or the codes of the ``.ped`` columns.

    Accepts the same codes as the ``pregnancy``, ``adopted``, ``childless``,
    ``art`` and ``gender_identity`` columns (case-insensitive).

    Raises:
        PedigreeParseError: If a value is not one of the accepted codes.
    """
    return Notation(
        pregnancy=_coerce_notation(pregnancy, "pregnancy", _PREGNANCY_CODES),
        adoption=_coerce_notation(adopted, "adopted", _ADOPTION_CODES),
        childlessness=_coerce_notation(childless, "childless", _CHILDLESS_CODES),
        reproductive_role=_coerce_notation(art, "art", _ART_CODES),
        gender=_coerce_notation(gender, "gender_identity", _GENDER_CODES),
        is_consultand=consultand,
    )


def _parse_notation(raw: dict[str, str], line_num: int) -> Notation:
    """Read the notation columns of one row into a Notation.

    Codes are case-insensitive; ``0``, ``-9``, ``NA``, ``.`` and empty mean
    nothing to note.

    Raises:
        PedigreeParseError: If a value is not one of the codes of its column.
    """
    return Notation(
        pregnancy=_notation_value(raw, "pregnancy", _PREGNANCY_CODES, line_num),
        adoption=_notation_value(raw, "adopted", _ADOPTION_CODES, line_num),
        childlessness=_notation_value(raw, "childless", _CHILDLESS_CODES, line_num),
        reproductive_role=_notation_value(raw, "art", _ART_CODES, line_num),
        gender=_notation_value(raw, "gender_identity", _GENDER_CODES, line_num),
        is_consultand=raw.get("consultand", "").strip() == "1",
    )


_TWIN_TYPE_MAP = {
    "MZ": TwinType.MONOZYGOTIC,
    "DZ": TwinType.DIZYGOTIC,
    "TW": TwinType.UNKNOWN,
}


def _parse_twin_field(raw: str, line_num: int) -> list[TwinInfo]:
    """Parse a twin field value into TwinInfo objects.

    Format: comma-separated TYPE:sample_id entries.
    Example: "MZ:twin_brother" or "DZ:sib1,DZ:sib2" for triplets.

    Args:
        raw: The raw field value (already stripped, not missing).
        line_num: Line number for error messages.

    Returns:
        List of TwinInfo objects.

    Raises:
        PedigreeParseError: If the format is invalid.
    """
    twin_infos: list[TwinInfo] = []
    entries = raw.split(",")
    for entry in entries:
        entry = entry.strip()
        if not entry or _is_missing(entry):
            continue
        if ":" not in entry:
            raise PedigreeParseError(
                f"Line {line_num}: invalid twin entry '{entry}', "
                f"expected TYPE:sample_id (e.g., MZ:twin_brother)"
            )
        twin_type_str, sibling_id = entry.split(":", 1)
        twin_type_str = twin_type_str.strip().upper()
        sibling_id = sibling_id.strip()

        if twin_type_str not in _TWIN_TYPE_MAP:
            raise PedigreeParseError(
                f"Line {line_num}: unknown twin type '{twin_type_str}', expected one of MZ, DZ, TW"
            )
        if not sibling_id:
            raise PedigreeParseError(f"Line {line_num}: empty sibling_id in twin entry '{entry}'")

        twin_infos.append(
            TwinInfo(
                twin_type=_TWIN_TYPE_MAP[twin_type_str],
                sibling_id=sibling_id,
            )
        )
    return twin_infos


def _resolve_twin_relationships(family: Family) -> None:
    """Ensure twin relationships are bidirectional and consistent.

    If individual A declares MZ:B, then B should also have MZ:A in its
    twin_info (if not already present). This makes twin groups easier to
    discover during layout.

    Also validates that twin siblings exist in the family and share
    the same parents.

    Args:
        family: The family to process. Modified in place.
    """
    for ind in list(family.individuals.values()):
        for ti in ind.twin_info:
            sibling = family.individuals.get(ti.sibling_id)
            if sibling is None:
                continue

            # Check if sibling already has the reciprocal relationship
            has_reciprocal = any(t.sibling_id == ind.sample_id for t in sibling.twin_info)
            if not has_reciprocal:
                sibling.twin_info.append(TwinInfo(twin_type=ti.twin_type, sibling_id=ind.sample_id))


def _infer_unions(family: Family) -> list[Union]:
    """Infer mating unions from parent references in a family.

    Two individuals form a union if they appear together as the father
    and mother of at least one child. By convention, the male partner
    is listed first in partner_ids when sexes are known.

    Args:
        family: A Family with individuals populated but unions empty.

    Returns:
        List of inferred Union objects.
    """
    # Collect unique parent pairs and their children
    pair_children: dict[tuple[str, str], list[str]] = {}

    for ind in family.individuals.values():
        if ind.father_id and ind.mother_id:
            # Normalize order: male first if known
            father = family.individuals.get(ind.father_id)
            mother = family.individuals.get(ind.mother_id)

            if father and mother:
                if father.sex == Sex.FEMALE and mother.sex == Sex.MALE:
                    pair = (ind.mother_id, ind.father_id)
                else:
                    pair = (ind.father_id, ind.mother_id)
            else:
                pair = (ind.father_id, ind.mother_id)

            pair_children.setdefault(pair, []).append(ind.sample_id)

    return [
        Union(partner_ids=pair, children_ids=children) for pair, children in pair_children.items()
    ]


def _get_ancestors(family: Family, sample_id: str) -> set[str]:
    """Collect all ancestor sample_ids for an individual via BFS.

    Traverses parent links upward from the given individual,
    collecting all reachable ancestors. Handles missing parents
    gracefully (stops traversal at unknown parents).

    Args:
        family: The family containing the individual.
        sample_id: The individual whose ancestors to find.

    Returns:
        Set of sample_ids of all ancestors (not including the individual).
    """
    ancestors: set[str] = set()
    queue: deque[str] = deque()

    ind = family.individuals.get(sample_id)
    if ind is None:
        return ancestors

    # Seed with direct parents
    if ind.father_id and ind.father_id in family.individuals:
        queue.append(ind.father_id)
    if ind.mother_id and ind.mother_id in family.individuals:
        queue.append(ind.mother_id)

    while queue:
        current_id = queue.popleft()
        if current_id in ancestors:
            continue
        ancestors.add(current_id)

        current = family.individuals.get(current_id)
        if current is None:
            continue
        if current.father_id and current.father_id in family.individuals:
            queue.append(current.father_id)
        if current.mother_id and current.mother_id in family.individuals:
            queue.append(current.mother_id)

    return ancestors


def _detect_consanguinity(family: Family) -> None:
    """Detect and mark consanguinous unions in a family.

    A union is consanguinous if both partners share at least one
    common ancestor. This is detected by computing the ancestor set
    for each partner and checking for intersection.

    Modifies unions in place by setting is_consanguinous=True
    where applicable.

    Args:
        family: The family whose unions to check. Modified in place.
    """
    for union in family.unions:
        p1_id, p2_id = union.partner_ids
        ancestors_p1 = _get_ancestors(family, p1_id)
        ancestors_p2 = _get_ancestors(family, p2_id)

        if ancestors_p1 & ancestors_p2:
            union.is_consanguinous = True


def _validate_family(family: Family) -> list[str]:
    """Validate a family's data integrity.

    Checks:
    - Parent references point to existing individuals or are None
    - No cycles in parent-child relationships
    - Sex consistency warnings (father should be male, mother female)

    Args:
        family: The family to validate.

    Returns:
        List of distinct warning messages (non-fatal issues), in the order
        they were found. A sex inconsistency is reported once per parent,
        not once per child of that parent.

    Raises:
        PedigreeValidationError: For fatal issues like cycles.
    """
    warning_msgs: list[str] = []

    for ind in family.individuals.values():
        # Check parent references
        if ind.father_id and ind.father_id not in family.individuals:
            warning_msgs.append(
                f"Individual '{ind.sample_id}': father '{ind.father_id}' "
                f"not found in family '{family.family_id}'"
            )
        if ind.mother_id and ind.mother_id not in family.individuals:
            warning_msgs.append(
                f"Individual '{ind.sample_id}': mother '{ind.mother_id}' "
                f"not found in family '{family.family_id}'"
            )

        # Sex consistency warnings
        if ind.father_id and ind.father_id in family.individuals:
            father = family.individuals[ind.father_id]
            if father.sex == Sex.FEMALE:
                warning_msgs.append(
                    f"Individual '{ind.father_id}' listed as father but has sex=FEMALE"
                )
        if ind.mother_id and ind.mother_id in family.individuals:
            mother = family.individuals[ind.mother_id]
            if mother.sex == Sex.MALE:
                warning_msgs.append(
                    f"Individual '{ind.mother_id}' listed as mother but has sex=MALE"
                )

    # Cycle detection via DFS
    _check_cycles(family)

    return list(dict.fromkeys(warning_msgs))


def _check_cycles(family: Family) -> None:
    """Detect cycles in parent-child relationships.

    Raises:
        PedigreeValidationError: If a cycle is detected.
    """
    visited: set[str] = set()
    in_stack: set[str] = set()

    def dfs(sample_id: str) -> None:
        if sample_id in in_stack:
            raise PedigreeValidationError(
                f"Cycle detected involving individual '{sample_id}' in family '{family.family_id}'"
            )
        if sample_id in visited:
            return

        visited.add(sample_id)
        in_stack.add(sample_id)

        ind = family.individuals.get(sample_id)
        if ind:
            if ind.father_id and ind.father_id in family.individuals:
                dfs(ind.father_id)
            if ind.mother_id and ind.mother_id in family.individuals:
                dfs(ind.mother_id)

        in_stack.discard(sample_id)

    for sample_id in family.individuals:
        if sample_id not in visited:
            dfs(sample_id)


def _create_phantom_parents(family: Family) -> None:
    """Create phantom (dashed) parents for individuals with only one known parent.

    If a child has a father but no mother (or vice versa), a phantom parent
    is created with a generated ID and dashed border.

    Args:
        family: The family to modify in place.
    """
    phantoms_to_add: list[Individual] = []
    updates: list[tuple[str, str, str]] = []  # (child_id, field, phantom_id)

    for ind in list(family.individuals.values()):
        has_father = ind.father_id is not None
        has_mother = ind.mother_id is not None

        if has_father and not has_mother:
            phantom_id = f"_phantom_mother_{ind.sample_id}"
            phantoms_to_add.append(
                Individual(
                    sample_id=phantom_id,
                    family_id=family.family_id,
                    sex=Sex.FEMALE,
                    is_phantom=True,
                )
            )
            updates.append((ind.sample_id, "mother_id", phantom_id))

        elif has_mother and not has_father:
            phantom_id = f"_phantom_father_{ind.sample_id}"
            phantoms_to_add.append(
                Individual(
                    sample_id=phantom_id,
                    family_id=family.family_id,
                    sex=Sex.MALE,
                    is_phantom=True,
                )
            )
            updates.append((ind.sample_id, "father_id", phantom_id))

    for phantom in phantoms_to_add:
        family.individuals[phantom.sample_id] = phantom

    for child_id, attr, phantom_id in updates:
        setattr(family.individuals[child_id], attr, phantom_id)


def _create_duplicate_nodes(family: Family) -> None:
    """Create duplicate nodes for individuals involved in multiple unions.

    When an individual has more than one union (multiple partners), the
    pedigree is drawn with separate subtrees for each relationship.
    The individual appears once as the "original" in the first union,
    and as duplicate nodes in subsequent unions. Duplicates are connected
    to the original by dotted arcs in the rendered diagram.

    The first union (in list order) keeps the original individual.
    Subsequent unions get a duplicate node with:
    - sample_id: "_dup_{original_id}_{partner_id}"
    - is_duplicate = True
    - duplicate_of = original sample_id
    - Same sex, phenotype, deceased, proband as original

    The union's partner_ids are updated to reference the duplicate.
    Children in that union have their parent reference updated too.

    Args:
        family: The family to process. Modified in place.
    """
    # Find individuals with multiple unions
    multi_union_ids: list[str] = []
    for sid in family.individuals:
        unions = family.get_unions_for_individual(sid)
        if len(unions) > 1:
            multi_union_ids.append(sid)

    for original_id in multi_union_ids:
        original = family.individuals[original_id]
        unions = family.get_unions_for_individual(original_id)

        # Keep the first union as-is, create duplicates for the rest
        for union in unions[1:]:
            p1, p2 = union.partner_ids
            partner_id = p2 if p1 == original_id else p1

            dup_id = f"_dup_{original_id}_{partner_id}"

            # Create the duplicate individual
            dup = Individual(
                sample_id=dup_id,
                family_id=family.family_id,
                sex=original.sex,
                phenotype=original.phenotype,
                is_deceased=original.is_deceased,
                is_proband=False,  # Proband marker only on original
                is_phantom=False,
                is_duplicate=True,
                is_carrier=original.is_carrier,
                is_obligate_carrier=original.is_obligate_carrier,
                # Drawn identically to the original, phenotype bands included.
                extra_phenotypes=dict(original.extra_phenotypes),
                extra_phenotypes_raw=dict(original.extra_phenotypes_raw),
                extra_phenotypes_categorical=dict(original.extra_phenotypes_categorical),
                duplicate_of=original_id,
                notation=replace(original.notation, is_consultand=False),
            )
            family.individuals[dup_id] = dup

            # Update union partner_ids to reference the duplicate
            if p1 == original_id:
                union.partner_ids = (dup_id, p2)
            else:
                union.partner_ids = (p1, dup_id)

            # Update children's parent references
            is_father = original.sex == Sex.MALE or (
                original.sex == Sex.UNKNOWN
                and any(
                    family.individuals.get(cid) and family.individuals[cid].father_id == original_id
                    for cid in union.children_ids
                )
            )
            for child_id in union.children_ids:
                child = family.individuals.get(child_id)
                if child is None:
                    continue
                if is_father and child.father_id == original_id:
                    child.father_id = dup_id
                elif not is_father and child.mother_id == original_id:
                    child.mother_id = dup_id


def parse_pedigree(file_path: str | Path) -> Pedigree:
    """Parse a pedigree from a file.

    Auto-detects separator (tab, comma, whitespace) and header presence.

    Args:
        file_path: Path to the pedigree file.

    Returns:
        A Pedigree object containing all parsed families.

    Raises:
        PedigreeParseError: If the file cannot be parsed.
        PedigreeValidationError: If the data fails validation.
        FileNotFoundError: If the file does not exist.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Pedigree file not found: {path}")

    text = path.read_text().strip()
    if not text:
        raise PedigreeParseError(f"Pedigree file is empty: {path}")

    return _parse_text(text)


def parse_pedigree_from_string(text: str) -> Pedigree:
    """Parse a pedigree from a string.

    Args:
        text: The pedigree data as a string.

    Returns:
        A Pedigree object.

    Raises:
        PedigreeParseError: If the text cannot be parsed.
        PedigreeValidationError: If the data fails validation.
    """
    text = text.strip()
    if not text:
        raise PedigreeParseError("Pedigree text is empty")
    return _parse_text(text)


def parse_pedigree_from_dataframe(df: pd.DataFrame) -> Pedigree:
    """Parse a pedigree from a pandas DataFrame.

    The DataFrame must have at least 6 columns corresponding to:
    family_id, sample_id, paternal_id, maternal_id, sex, phenotype.
    Additional columns (7-10) are mapped to: date_of_birth, twins, deceased, proband.

    Column names are ignored — column position determines meaning.

    Args:
        df: A pandas DataFrame with at least 6 columns.

    Returns:
        A Pedigree object.

    Raises:
        PedigreeParseError: If the DataFrame has fewer than 6 columns or is empty.
        PedigreeValidationError: If the data fails validation.
    """
    try:
        import pandas  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "pandas is required for parse_pedigree_from_dataframe(). "
            "Install it with: pip install pypedigree[pandas]"
        ) from e

    if len(df.columns) < 6:
        raise PedigreeParseError(f"DataFrame must have at least 6 columns, got {len(df.columns)}")

    if len(df) == 0:
        raise PedigreeParseError("DataFrame is empty")

    # Replace NaN/NA with "0" (the standard missing value) before converting to TSV
    df_clean = df.fillna("0").astype(str)
    text = df_clean.to_csv(sep="\t", index=False, header=False)

    return _parse_text(text.strip())


def _parse_text(text: str, *, max_columns: int | None = None) -> Pedigree:
    """Internal parser that works on text content.

    Args:
        text: Non-empty, stripped pedigree text.
        max_columns: If set, truncate each row to this many columns.
            Used for PLINK .fam (6) and PLINK .ped (6) to ignore
            genotype data and prevent mis-parsing extra columns.

    Returns:
        Parsed Pedigree.
    """
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        raise PedigreeParseError("No data lines found")

    # Detect separator from first line
    separator = _detect_separator(lines[0])

    # Detect and skip header, capture extra phenotype column names
    start_line = 0
    phenotype_names: list[str] = []
    annotation_names: list[str] = []
    phenotype_display_modes: dict[str, PhenotypeDisplayMode] = {}
    carrier_col_index: int | None = None
    col_map: ColumnMapping | None = None
    use_mapped_parser = False

    if max_columns is None and _detect_header(lines[0], separator):
        header_fields = _split_line(lines[0], separator)
        col_map = _parse_header_columns(header_fields)

        if col_map.annotation_cols:
            annotation_names = [name for _, name in col_map.annotation_cols]
        if (
            col_map.extra_phenotype_cols
            or col_map.annotation_cols
            or col_map.notation_cols
            or not col_map.is_positional
        ):
            # Header has extra phenotype/annotation columns, or structural
            # columns out of the standard order → read columns by name
            use_mapped_parser = True
            phenotype_names = [name for _, name in col_map.extra_phenotype_cols]
            phenotype_display_modes = dict(col_map.phenotype_display_modes)
        elif len(header_fields) > 10:
            # Legacy path: standard 10-col header with extras at col 11+
            raw_names = [f.strip() for f in header_fields[10:]]
            for j, name in enumerate(raw_names):
                if name.lower() == "carrier":
                    carrier_col_index = 10 + j
                else:
                    stripped, mode = _strip_display_prefix(name)
                    phenotype_names.append(stripped)
                    if mode != PhenotypeDisplayMode.AUTO:
                        phenotype_display_modes[stripped] = mode

        start_line = 1
        if len(lines) < 2:
            raise PedigreeParseError("File contains only a header, no data")

    # Parse individuals and group by family
    families: dict[str, Family] = {}

    for i, line in enumerate(lines[start_line:], start=start_line + 1):
        fields = _split_line(line, separator)
        if max_columns is not None:
            fields = fields[:max_columns]

        if use_mapped_parser and col_map is not None:
            ind = _parse_individual_mapped(fields, line_num=i, col_map=col_map)
        else:
            ind = _parse_individual(
                fields,
                line_num=i,
                phenotype_names=phenotype_names,
                carrier_col_index=carrier_col_index,
            )

        if ind.family_id not in families:
            families[ind.family_id] = Family(family_id=ind.family_id)

        family = families[ind.family_id]

        if ind.sample_id in family.individuals:
            raise PedigreeParseError(
                f"Line {i}: duplicate sample_id '{ind.sample_id}' in family '{ind.family_id}'"
            )

        family.individuals[ind.sample_id] = ind

    # Post-processing for each family
    for family in families.values():
        # Create phantom parents for single-parent children
        _create_phantom_parents(family)

        # Resolve twin relationships (make bidirectional)
        _resolve_twin_relationships(family)

        # Validate
        warnings_list = _validate_family(family)
        for msg in warnings_list:
            warnings.warn(msg, PedigreeDataWarning, stacklevel=2)

        # Infer unions
        family.unions = _infer_unions(family)

        # Detect consanguinous unions (shared ancestors)
        _detect_consanguinity(family)

        # Create duplicate nodes for individuals with multiple unions
        _create_duplicate_nodes(family)

    pedigree = Pedigree(
        families=families,
        phenotype_names=phenotype_names,
        phenotype_display_modes=phenotype_display_modes,
        annotation_names=annotation_names,
    )

    # Normalize inline phenotypes if present
    if phenotype_names:
        _normalize_extra_phenotypes(pedigree, phenotype_names)

    return pedigree


def _normalize_extra_phenotypes(pedigree: Pedigree, pheno_names: list[str]) -> None:
    """Normalize extra_phenotypes_raw values to [0, 1] and store in extra_phenotypes.

    Uses display mode hints from column prefixes when available:
    - BINARY/OUTLINE: force binary interpretation
    - GRADIENT: force min-max normalization
    - CATEGORICAL: skip numeric normalization (uses string values)
    - AUTO: auto-detect from data (binary if {0,1} or {0,1,2})

    Args:
        pedigree: Pedigree with raw phenotype data. Modified in place.
        pheno_names: List of phenotype column names.
    """
    modes = pedigree.phenotype_display_modes

    # Collect all raw values per phenotype across all individuals
    all_values: dict[str, list[float]] = {n: [] for n in pheno_names}
    all_individuals: list[Individual] = []
    for family in pedigree.families.values():
        for ind in family.individuals.values():
            all_individuals.append(ind)
            for name in pheno_names:
                val = ind.extra_phenotypes_raw.get(name)
                if val is not None:
                    all_values[name].append(val)

    # Determine per-phenotype normalization
    for name in pheno_names:
        mode = modes.get(name, PhenotypeDisplayMode.AUTO)

        # Categorical: no numeric normalization needed
        if mode == PhenotypeDisplayMode.CATEGORICAL:
            for ind in all_individuals:
                ind.extra_phenotypes[name] = None  # Not used for categorical
            continue

        # Carrier: keep the 0/1/2 code itself; NaN marks a missing value.
        if mode == PhenotypeDisplayMode.CARRIER:
            for ind in all_individuals:
                raw = ind.extra_phenotypes_raw.get(name)
                ind.extra_phenotypes[name] = float("nan") if raw is None else raw
            continue

        vals = all_values[name]
        if not vals:
            continue

        # Determine effective normalization strategy
        if mode in (PhenotypeDisplayMode.BINARY, PhenotypeDisplayMode.OUTLINE):
            # Forced binary: treat as pedigree encoding if values include 2
            strategy = "binary"
        elif mode == PhenotypeDisplayMode.GRADIENT:
            strategy = "gradient"
        else:
            # AUTO: detect from data
            is_binary_01 = all(v in (0.0, 1.0) for v in vals)
            is_pedigree_encoding = all(v in (0.0, 1.0, 2.0) for v in vals) and any(
                v == 2.0 for v in vals
            )
            if is_binary_01:
                strategy = "binary_01"
            elif is_pedigree_encoding:
                strategy = "binary"
            else:
                strategy = "gradient"

        # NaN sentinel means "NA / no data" (-9), None means "unknown" (0) → grey
        _NA = float("nan")

        if strategy == "binary_01":
            for ind in all_individuals:
                raw = ind.extra_phenotypes_raw.get(name)
                if raw is None:
                    ind.extra_phenotypes[name] = _NA  # -9 → NA
                else:
                    ind.extra_phenotypes[name] = raw  # 0.0 or 1.0
        elif strategy == "binary":
            # Remap: 0 → None (unknown/grey), 1 → 0 (unaffected), 2 → 1 (affected)
            # raw=None (-9) → NaN (NA / no data)
            for ind in all_individuals:
                raw = ind.extra_phenotypes_raw.get(name)
                if raw is None:
                    ind.extra_phenotypes[name] = _NA  # -9 → NA
                elif raw >= 1.0:
                    ind.extra_phenotypes[name] = raw - 1.0  # 1→0, 2→1
                else:
                    ind.extra_phenotypes[name] = None  # 0 = unknown → grey
        else:
            # Gradient: min-max normalization to [0, 1]
            lo, hi = min(vals), max(vals)
            rng = hi - lo if hi != lo else 1.0
            for ind in all_individuals:
                raw = ind.extra_phenotypes_raw.get(name)
                if raw is None:
                    ind.extra_phenotypes[name] = _NA  # -9 → NA
                else:
                    ind.extra_phenotypes[name] = (raw - lo) / rng


def parse_phenotype_file(
    file_path: str | Path,
) -> tuple[dict[str, dict[str, float | None]], list[str]]:
    """Parse a separate phenotype file (.pheno TSV).

    Expected format: first row is header (sample_id, pheno1, pheno2, ...),
    subsequent rows are data. First column is sample_id.

    Args:
        file_path: Path to the .pheno file.

    Returns:
        Tuple of (phenotype_data, phenotype_names) where phenotype_data
        maps sample_id -> {pheno_name: raw_value}.
    """
    path = Path(file_path)
    text = path.read_text().strip()
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise PedigreeParseError(f"Phenotype file must have header + data: {path}")

    # Detect separator
    sep: str | None = None
    for s in ["\t", ","]:
        if len(lines[0].split(s)) >= 2:
            sep = s
            break

    header = lines[0].split(sep) if sep else lines[0].split()
    pheno_names = [h.strip() for h in header[1:]]

    data: dict[str, dict[str, float | None]] = {}
    for line in lines[1:]:
        parts = line.split(sep) if sep else line.split()
        sample_id = parts[0].strip()
        values: dict[str, float | None] = {}
        for j, name in enumerate(pheno_names):
            idx = j + 1
            if idx < len(parts):
                raw = parts[idx].strip()
                if _is_pheno_missing(raw):
                    values[name] = None
                else:
                    try:
                        values[name] = float(raw)
                    except ValueError:
                        values[name] = None
            else:
                values[name] = None
        data[sample_id] = values

    return data, pheno_names


def parse_annotation_file(
    file_path: str | Path,
) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Parse an external annotation file (TSV).

    Expected format: first row is header (sample_id, ann1, ann2, ...),
    subsequent rows are data. First column is sample_id.

    Args:
        file_path: Path to the annotation TSV file.

    Returns:
        Tuple of (annotation_data, annotation_names) where annotation_data
        maps sample_id -> {ann_name: value_string}.
    """
    path = Path(file_path)
    text = path.read_text().strip()
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise PedigreeParseError(f"Annotation file must have header + data: {path}")

    sep: str | None = None
    for s in ["\t", ","]:
        if len(lines[0].split(s)) >= 2:
            sep = s
            break

    header = lines[0].split(sep) if sep else lines[0].split()
    ann_names = [h.strip() for h in header[1:]]
    _ann_missing = frozenset({"0", "-9", "nan", "na", ".", ""})

    data: dict[str, dict[str, str]] = {}
    for line in lines[1:]:
        parts = line.split(sep) if sep else line.split()
        sample_id = parts[0].strip()
        values: dict[str, str] = {}
        for j, name in enumerate(ann_names):
            idx = j + 1
            if idx < len(parts):
                raw = parts[idx].strip()
                if raw.lower() not in _ann_missing:
                    values[name] = raw
        data[sample_id] = values

    return data, ann_names


def apply_annotations_to_pedigree(
    pedigree: Pedigree,
    ann_data: dict[str, dict[str, str]],
    ann_names: list[str],
) -> None:
    """Apply external annotation data to a parsed pedigree.

    Args:
        pedigree: Pedigree to update.
        ann_data: Map of sample_id -> {ann_name: value}.
        ann_names: Ordered annotation names.
    """
    for name in ann_names:
        if name not in pedigree.annotation_names:
            pedigree.annotation_names.append(name)

    for family in pedigree.families.values():
        for sid, ind in family.individuals.items():
            if sid in ann_data:
                ind.annotations.update(ann_data[sid])


def parse_fam_file(file_path: str | Path) -> Pedigree:
    """Parse a PLINK .fam file.

    The .fam format is a whitespace-delimited file with exactly 6 columns
    (FID, IID, PAT, MAT, SEX, PHENO) and no header row. This is the standard
    format used by PLINK, KING, and other genetics tools.

    Args:
        file_path: Path to the .fam file.

    Returns:
        A Pedigree object containing all parsed families.

    Raises:
        PedigreeParseError: If the file cannot be parsed.
        PedigreeValidationError: If the data fails validation.
        FileNotFoundError: If the file does not exist.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"FAM file not found: {path}")

    text = path.read_text().strip()
    if not text:
        raise PedigreeParseError(f"FAM file is empty: {path}")

    return _parse_text(text, max_columns=6)


def parse_plink_ped(file_path: str | Path) -> Pedigree:
    """Parse a PLINK .ped file, ignoring genotype columns.

    PLINK .ped files have 6 standard pedigree columns followed by genotype
    data (pairs of alleles). This function reads only the first 6 columns
    and discards genotype data. No header row is expected.

    Args:
        file_path: Path to the PLINK .ped file.

    Returns:
        A Pedigree object containing all parsed families.

    Raises:
        PedigreeParseError: If the file cannot be parsed.
        PedigreeValidationError: If the data fails validation.
        FileNotFoundError: If the file does not exist.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"PLINK .ped file not found: {path}")

    text = path.read_text().strip()
    if not text:
        raise PedigreeParseError(f"PLINK .ped file is empty: {path}")

    return _parse_text(text, max_columns=6)


def export_fam(pedigree: Pedigree, file_path: str | Path) -> Path:
    """Export a Pedigree to PLINK .fam format.

    Writes a whitespace-delimited file with 6 columns:
    FID, IID, PAT, MAT, SEX, PHENO. No header row.

    Phantom parents, duplicate nodes, and other internal constructs
    are excluded from the output.

    Args:
        pedigree: The Pedigree to export.
        file_path: Output file path.

    Returns:
        The Path of the written file.
    """
    path = Path(file_path)
    lines: list[str] = []

    for family in pedigree.families.values():
        for ind in family.individuals.values():
            if ind.is_phantom or ind.is_duplicate:
                continue

            fid = ind.family_id
            iid = ind.sample_id
            # Resolve internal IDs back to original
            father = ind.father_id
            if father:
                if father.startswith("_phantom_"):
                    father = "0"
                elif father.startswith("_dup_"):
                    fi = family.individuals.get(father)
                    father = fi.duplicate_of if fi and fi.duplicate_of else "0"
            else:
                father = "0"
            pat = father

            mother = ind.mother_id
            if mother:
                if mother.startswith("_phantom_"):
                    mother = "0"
                elif mother.startswith("_dup_"):
                    mi = family.individuals.get(mother)
                    mother = mi.duplicate_of if mi and mi.duplicate_of else "0"
            else:
                mother = "0"
            mat = mother

            sex_map = {Sex.MALE: "1", Sex.FEMALE: "2", Sex.UNKNOWN: "0"}
            sex = sex_map.get(ind.sex, "0")

            pheno_map = {Phenotype.UNAFFECTED: "1", Phenotype.AFFECTED: "2", Phenotype.UNKNOWN: "0"}
            pheno = pheno_map.get(ind.phenotype, "0")

            lines.append(f"{fid}\t{iid}\t{pat}\t{mat}\t{sex}\t{pheno}")

    path.write_text("\n".join(lines) + "\n")
    return path


def apply_phenotypes_to_pedigree(
    pedigree: Pedigree,
    pheno_data: dict[str, dict[str, float | None]],
    pheno_names: list[str],
) -> None:
    """Apply external phenotype data to individuals and normalize.

    Args:
        pedigree: The pedigree to update. Modified in place.
        pheno_data: Map of sample_id -> {pheno_name: raw_value}.
        pheno_names: Ordered list of phenotype names.
    """
    pedigree.phenotype_names = pheno_names

    for family in pedigree.families.values():
        for ind in family.individuals.values():
            sid = ind.duplicate_of if ind.is_duplicate else ind.sample_id
            if sid in pheno_data:
                ind.extra_phenotypes_raw = dict(pheno_data[sid])

    _normalize_extra_phenotypes(pedigree, pheno_names)
