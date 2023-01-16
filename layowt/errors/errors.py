"""Errors and Warnings."""

class LayowtError(Exception):
    """Root exception class"""

class EmptyLayoutError(AttributeError):
    """Raised when a layout geometry is empty."""

class LayoutSchemaError(LayowtError):
    """Raised when a database schema does not have all the required system tables for Layout upload."""