"""Utility functions supporting the Gromacs Python interface.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

__all__ = []

import gmx.fileio as fileio
from gmx.exceptions import UsageError
from gmx.exceptions import FileError

def _filetype(filename):
    """Use Gromacs file I/O module to identify known file types.

    Attempts to open the provided file. If successful, returns the gmx.io class
    associated with the identified file type. Relative file paths
    should be interpreted with regard to the execution context, but
    are currently interpreted relative to the directory set in the
    ``_current_dir`` attribute of the imported gmx.fileio module
    whereever the Python script is running, which can be reassigned
    by the user.

    This function cannot be part of the public interface until the
    meaning of file paths is clarified regarding the presence of
    execution contexts.

    Args:
        filename (str): the name of a file for Gromacs to inspect

    Returns:
        gmx.io.File subclass

    Raises:
        FileError is filename is not a readable file.
    """
    # Get an absolute filename.
    try:
        filepath = os.path.join(fileio._current_dir, filename)
        if os.path.isfile(filepath):
            return fileio.TprFile
        else:
            raise UsageError("Argument is not a readable file.")
    except:
        # This would be a good place to chain exceptions. Let's see where
        # Python 2/3 compatibility falls.
        raise UsageError("Cannot process filename argument.")

def _test():
        import doctest
        doctest.testmod()

if __name__ == "__main__":
    _test()
