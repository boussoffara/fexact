# type: ignore[attr-defined]
"""Awesome `fisher-acm` is a Python cli/package created with https://github.com/TezRomacH/python-package-template"""

try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:  # pragma: no cover
    from importlib_metadata import version, PackageNotFoundError


try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
