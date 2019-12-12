#
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
import warnings

__version__ = "0.1.0"

# TODO: (pending infrastructure and further discussion) Configure with CMake.
# __version__ = "@PROJECT_VERSION@"
# major = @PROJECT_VERSION_MAJOR@
# minor = @PROJECT_VERSION_MINOR@
# patch = @PROJECT_VERSION_PATCH@

from gmxapi.exceptions import FeatureNotAvailableError

major = 0
minor = 1
patch = 0

# Note: this is not automatically updated. See RELEASE.txt and https://github.com/kassonlab/gmxapi/issues/152
release = True

# Features added since the initial gmxapi prototype, targeted for version 0.1.
_named_features_0_0 = ['fr1', 'fr3', 'fr7', 'fr15']
# Features named since the finalization of the 0.1 specification with GROMACS 2020.
_named_features_0_1 = []
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
#
# Features consisting of 'fr' and a numeric suffix are the functional requirements
# described in roadmap.rst, as described at https://redmine.gromacs.org/issues/2893
#
# fr1: wrap importable Python code.
# fr2: output proxy establishes execution dependency (superseded by fr3)
# fr3: output proxy can be used as input
# fr4: dimensionality and typing of named data causes generation of correct work topologies
# fr5: explicit many-to-one or many-to-many data flow
# fr7: Python bindings for launching simulations
# fr8: gmx.mdrun understands ensemble work
# fr9: MD plugins
# fr10: fused operations for use in looping constructs
# fr11: Python access to TPR file contents
# fr12: Simulation checkpoint handling
# fr13: ``run`` module function simplifies user experience
# fr14: Easy access to GROMACS run time parameters
# fr15: Simulation input modification
# fr16: Create simulation input from simulation output
# fr17: Prepare simulation input from multiple sources
# fr18: GROMACS CLI tools receive improved Python-level support over generic commandline_operations
# fr19: GROMACS CLI tools receive improved C++-level support over generic commandline_operations
# fr20: Python bindings use C++ API for expressing user interface
# fr21 User insulated from filesystem paths
# fr22 MPI-based ensemble management from Python
# fr23 Ensemble simulations can themselves use MPI


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
    if major > major_version:
        return True
    elif major == major_version and minor >= minor_version:
        return True
    elif major == major_version and minor == minor_version and patch >= patch_version:
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
    if api_is_at_least(0, 2):
        # For sufficiently advanced API versions, we want to warn that old
        # feature checks lose meaning and should no longer be checked.
        # We provide a suggestion with the API version that absorbed their
        # specification.
        if name in _named_features_0_0:
            warnings.warn(
                'Old feature name. Use `api_is_at_least(0, 1)` instead of `has_feature({})`.'.format(name),
                category=DeprecationWarning
            )

    # Check whether the feature is listed in the API specification amendments.
    if name in _named_features_0_0 + _named_features_0_1:
        return True
    else:
        if enable_exception:
            raise FeatureNotAvailableError('Feature {} not available.'.format(str(name)))
        return False
