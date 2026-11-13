"""Custom exceptions for pypedigree."""


class PedigreeError(Exception):
    """Base exception for all pypedigree errors."""


class PedigreeParseError(PedigreeError):
    """Raised when a pedigree file cannot be parsed."""


class PedigreeValidationError(PedigreeError):
    """Raised when pedigree data fails validation checks."""


class PedigreeLayoutError(PedigreeError):
    """Raised when the layout engine encounters an unresolvable situation."""


class PedigreeWarning(UserWarning):
    """Category for non-fatal pypedigree diagnostics.

    Warnings use this category so callers can filter them independently of
    other libraries' warnings, e.g.::

        warnings.simplefilter("error", PedigreeWarning)
    """


class PedigreeStyleWarning(PedigreeWarning):
    """Raised when a style definition contains something unusable."""


class PedigreeDataWarning(PedigreeWarning):
    """Raised when pedigree data is questionable but still renderable."""
