class PGFragError(Exception):
    """Base error for invalid PG-Frag inputs or unavailable capabilities."""


class ValidationError(PGFragError):
    pass


class UnsupportedCapability(PGFragError):
    pass
