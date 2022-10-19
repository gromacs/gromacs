#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2022- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

"""Provide some additional utilities."""
import functools
import os
from types import MappingProxyType

from gmxapi.operation import function_wrapper


@function_wrapper(output={"path": str})
def join_path(first: str, second: str, output=None):
    """Get a Future path for use in data flow.

    This is useful when a base path or filename is not known until runtime.

    Attributes:
        output.path (str): *first* and *second*, joined by the native filesystem path separator.
    """
    output.path = os.path.join(first, second)


@functools.lru_cache()
def config():
    """Get the GROMACS configuration detected during installation.

    Returns read-only dictionary proxy to file written during installation.
    The :py:class:`~typing.Mapping` contains information about the supporting
    GROMACS installation that was used to configure the Python package
    installation. The exact keys in the Mapping is not formally specified,
    and mostly for internal use.

    .. versionadded:: 0.4

    """
    import json
    from importlib.resources import open_text

    with open_text("gmxapi", "gmxconfig.json") as textfile:
        _config = json.load(textfile)
    return MappingProxyType(_config)
