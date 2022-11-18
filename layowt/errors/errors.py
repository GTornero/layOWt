"""Errors and Warnings."""

class LayowtError(Exception):
    """Root exception class"""

class EmptyLayoutError(AttributeError):
    """Raised when a layout geometry is empty."""