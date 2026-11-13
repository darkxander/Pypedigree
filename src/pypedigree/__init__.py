"""pypedigree - Generate pedigree (family tree) diagrams from pedigree files."""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _metadata_version

from pypedigree.api import PedigreeBuilder
from pypedigree.config import (
    DEFAULT_PHENOTYPE_PALETTE,
    PedigreeStyle,
    apply_pgs_spacing,
    is_pgs_column,
    resolve_style,
)
from pypedigree.errors import (
    PedigreeDataWarning,
    PedigreeError,
    PedigreeLayoutError,
    PedigreeParseError,
    PedigreeStyleWarning,
    PedigreeValidationError,
    PedigreeWarning,
)
from pypedigree.layout import compute_layout, compute_layout_circular
from pypedigree.models import (
    Adoption,
    Childlessness,
    DuplicateArc,
    EdgeLayout,
    EdgeType,
    Family,
    FamilyLayout,
    Gender,
    Individual,
    LineStyle,
    NodeLayout,
    Notation,
    Pedigree,
    Phenotype,
    PhenotypeDisplayMode,
    Pregnancy,
    ReproductiveRole,
    Sex,
    TwinConnector,
    TwinType,
    Union,
)
from pypedigree.parser import (
    apply_annotations_to_pedigree,
    apply_phenotypes_to_pedigree,
    build_notation,
    export_fam,
    parse_annotation_file,
    parse_fam_file,
    parse_pedigree,
    parse_pedigree_from_dataframe,
    parse_pedigree_from_string,
    parse_phenotype_file,
    parse_plink_ped,
)
from pypedigree.renderer import render_family, render_pedigree
from pypedigree.svg_renderer import render_interactive_svg

try:
    __version__ = _metadata_version("pypedigree")
except PackageNotFoundError:  # pragma: no cover - running from a source tree
    __version__ = "0.0.0.dev0"

__all__ = [
    # Models
    "Individual",
    "Union",
    "Family",
    "Pedigree",
    "Sex",
    "Phenotype",
    "PhenotypeDisplayMode",
    "TwinType",
    "Notation",
    "Pregnancy",
    "Adoption",
    "Childlessness",
    "ReproductiveRole",
    "Gender",
    "NodeLayout",
    "EdgeLayout",
    "EdgeType",
    "LineStyle",
    "TwinConnector",
    "DuplicateArc",
    "FamilyLayout",
    # Parser
    "parse_pedigree",
    "parse_pedigree_from_dataframe",
    "parse_pedigree_from_string",
    "parse_phenotype_file",
    "apply_phenotypes_to_pedigree",
    "parse_annotation_file",
    "apply_annotations_to_pedigree",
    "parse_fam_file",
    "parse_plink_ped",
    "export_fam",
    "build_notation",
    # Layout
    "compute_layout",
    "compute_layout_circular",
    # Renderer
    "render_pedigree",
    "render_family",
    "render_interactive_svg",
    # Config
    "PedigreeStyle",
    "DEFAULT_PHENOTYPE_PALETTE",
    "resolve_style",
    "apply_pgs_spacing",
    "is_pgs_column",
    # API
    "PedigreeBuilder",
    # Errors
    "PedigreeError",
    "PedigreeParseError",
    "PedigreeValidationError",
    "PedigreeLayoutError",
    # Warnings
    "PedigreeWarning",
    "PedigreeStyleWarning",
    "PedigreeDataWarning",
    # Version
    "__version__",
]
