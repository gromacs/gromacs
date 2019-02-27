# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.
"""
Provide version and release information.

Attributes:
    major (int): gmxapi major version number.
    minor (int): gmxapi minor version number.
    patch (int): gmxapi patch level number.
    release (bool): True if imported gmx module is an officially tagged release, else False.

"""

# TODO: configure with CMake
# __version__ = "@PROJECT_VERSION@"
# major = @PROJECT_VERSION_MAJOR@
# minor = @PROJECT_VERSION_MINOR@
# patch = @PROJECT_VERSION_PATCH@

from gmxapi.exceptions import Error as GmxapiError

__version__ = "0.1.0-dev0"
major = 0
minor = 1
patch = 0

# Note: this is not automatically updated. See RELEASE.txt and https://github.com/kassonlab/gmxapi/issues/152
release = False

_named_features = []


class FeatureError(GmxapiError):
    """Module exception to indicate missing features."""


def api_is_at_least(major_version, minor_version=0, patch_version=0):
    """Allow client to check whether installed module supports the requested API level.

    Arguments:
        major_version (int): gmxapi major version number.
        minor_version (int): optional gmxapi minor version number (default: 0).
        patch_version (int): optional gmxapi patch level number (default: 0).

    Returns:
        True if installed gmx package is greater than or equal to the input level

    Note that if gmxapi.version.release is False, the package is not guaranteed to correctly or
    fully support the reported API level.
    """
    if not isinstance(major_version, int) or not isinstance(minor_version, int) or not isinstance(patch_version, int):
        raise TypeError('Version levels must be provided as integers.')
    if major >= major_version:
        return True
    elif major == major_version and minor >= minor_version:
        return True
    elif major == major_version and minor == minor_version and patch >= patch_version:
        return True
    else:
        return False


def has_feature(name='', enable_exception=False):
    """Query whether a named feature is available in the installed package.

    Between updates to the API specification, new features or experimental aspects
    may be introduced into the package and need to be detectable. This function
    is intended to facilitate code testing and resolving differences between
    development branches. Users should refer to the documentation for the package
    modules and API level.

    The function can be used to get a boolean result or can be used to raise an
    exception in code block by setting `enable_exception=True`

    Returns:
        True if named feature is recognized by the installed package, else False.

    Raises:
        gmxapi.version.FeatureError if `enable_exception == True` and feature is not found.

    """
    if name in _named_features:
        return True
    else:
        if enable_exception:
            raise FeatureError('Feature {} not available.'.format(str(name)))
        return False
