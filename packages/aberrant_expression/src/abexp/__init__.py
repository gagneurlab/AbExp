"""The Python packages of the AbExp pipeline.

This module imports none of the subpackages, so that importing one subpackage does not import the dependencies of
the others.
"""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version('aberrant-expression')
except PackageNotFoundError:
    # not installed, e.g. imported from a checkout through PYTHONPATH
    __version__ = 'unknown'
