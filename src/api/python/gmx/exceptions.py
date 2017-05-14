"""Exceptions and Warnings raised by gmx module operations.

Errors, warnings, and other exceptions used in the Gromacs
gmx Python package are defined in the gmx.exceptions submodule.

The Gromacs gmx Python package defines a root exception,
gmx.exceptions.Error, from which all Exceptions thrown from
within the module should derive. If a published component of
the gmx package throws an execption that cannot be caught
as a gmx.exceptions.Error, please report the bug.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


__all__ = ['Error', 'UsageError']

class Error(Exception):
    """Base exception for gmx.exceptions classes."""

class UsageError(Error):
    """Unsupported syntax or call signatures.

    Generic usage error for Gromacs gmx module.
    """

class FileError(Error):
    """Problem with a file or filename."""
