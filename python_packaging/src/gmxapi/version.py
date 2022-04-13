#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2019- The GROMACS Authors
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

"""
gmxapi version and release information.

The ``gmxapi.__version__`` attribute contains a :pep:`version string <440>`.
The more general way to access the package version is with the
:py:mod:`pkg_resources <https://setuptools.readthedocs.io/en/latest/pkg_resources.html>` module::

    pkg_resources.get_distribution('gmxapi').version

`gmxapi.version` module functions `api_is_at_least()` and `has_feature()`
support additional convenience and introspection.

.. versionchanged:: 0.2

    This module no longer provides public data attributes.
    Instead, use the module functions or
    :py:mod:`packaging.version <https://packaging.pypa.io/en/latest/version/>`.

.. seealso::

    Consider https://packaging.pypa.io/en/latest/version/ for programmatic
    handling of the version string. For example::

        from packaging.version import parse
        gmxapi_version = pkg_resources.get_distribution('gmxapi').version
        if parse(gmxapi_version).is_prerelease:
            print('The early bird gets the worm.')

.. todo:: Use pkg_resources.get_distribution('gmxapi').version and
          "development installations" instead of relying on or publicizing
          a __version__ attribute.
"""
import warnings

from .exceptions import FeatureNotAvailableError

# TODO(#3851): Version management policy and procedures.
_major = 0
_minor = 4
_micro = 0
_suffix = 'a1'

# Reference https://www.python.org/dev/peps/pep-0440/
# and https://packaging.pypa.io/en/latest/version/
__version__ = '{major}.{minor}.{micro}{suffix}'.format(major=_major,
                                                       minor=_minor,
                                                       micro=_micro,
                                                       suffix=_suffix)

# Named features describe functionality or behavior introduced since the last
# major release, and should be described in gmxapi documentation or issue
# tracking system. Note that, as features become part of the specification,
# conditionals depending on them should be phased out of the package source. At
# major releases, the named feature list should be reset to empty. Optionally,
# we could raise a DeprecationWarning for calls to has_feature() for features
# that have become part of the specification, at least for a few minor release or
# a few years, to avoid introducing errors to client code.
#
# Bugs and bug fixes may be indicated with names consisting of tracked issue URLs.

# Named features for gmxapi 0.x (pre-1.0 versions).
_named_features_0 = [[]] * (_minor + 1)

# Features added since the initial gmxapi prototype, targeted for version 0.1.
# Functional requirements were described in issues #2045 and #2893. See also
# https://gitlab.com/gromacs/gromacs/-/blob/release-2020/python_packaging/roadmap.rst
_named_features_0[0] = ['fr1', 'fr3', 'fr7', 'fr15']
# Features named since the finalization of the 0.1 specification with GROMACS 2020.
_named_features_0[1] = []

_named_features_0[2] = [
    'container_futures',
    'mdrun_checkpoint_output',
    'mdrun_runtime_args',
]

_named_features_0[3] = [
    'cli_env_kwarg',
]

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
    if not isinstance(major_version, int) or not isinstance(minor_version, int) or not isinstance(
            patch_version,
            int):
        raise TypeError('Version levels must be provided as integers.')
    if _major > major_version:
        return True
    elif _major == major_version and _minor >= minor_version:
        return True
    elif _major == major_version and _minor == minor_version and _micro >= patch_version:
        return True
    else:
        return False


def has_feature(name='', enable_exception=False) -> bool:
    """Query whether a named feature is available in the installed package.

    Between updates to the API specification, new features or experimental aspects
    may be introduced into the package and need to be detectable. This function
    is intended to facilitate code testing and resolving differences between
    development branches. Users should refer to the documentation for the package
    modules and API level.

    The primary use case is, in conjunction with `api_is_at_least()`, to allow
    client code to robustly identify expected behavior and API support through
    conditional execution and branching. Note that behavior is strongly
    specified by the API major version number. Features that have become part of
    the specification and bug-fixes referring to previous major versions should
    not be checked with *has_feature()*. Using *has_feature()* with old feature
    names will produce a DeprecationWarning for at least one major version, and
    client code should be updated to avoid logic errors in future versions.

    For convenience, setting ``enable_exception = True`` causes the function to
    instead raise a gmxapi.exceptions.FeatureNotAvailableError for unrecognized feature names.
    This allows extension code to cleanly produce a gmxapi exception instead of
    first performing a boolean check. Also, some code may be unexecutable for
    more than one reason, and sometimes it is cleaner to catch all
    `gmxapi.exceptions.Error` exceptions for a code block, rather than to
    construct complex conditionals.

    Returns:
        True if named feature is recognized by the installed package, else False.

    Raises:
        gmxapi.exceptions.FeatureNotAvailableError: If ``enable_exception == True`` and feature is not found.

    """
    # First, issue a warning if the feature name is subject to removal because
    # of the history of the API specification.
    for version in range(_minor):
        # For sufficiently advanced API versions, we want to warn that old
        # feature checks lose meaning and should no longer be checked.
        # We provide a suggestion with the API version that absorbed their
        # specification.
        if name in _named_features_0[version]:
            warnings.warn(
                f'Old feature name. Use `api_is_at_least({_major}, {version + 1})` instead of `has_feature({name})`.',
                category=DeprecationWarning,
                stacklevel=2
            )

    # Check whether the feature is listed in the API specification amendments.
    if any(name in features for features in _named_features_0):
        return True
    else:
        if enable_exception:
            raise FeatureNotAvailableError('Feature {} not available.'.format(str(name)))
        return False
