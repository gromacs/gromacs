#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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

"""Reusable definitions for test modules.

Provides utilities and pytest fixtures for gmxapi and GROMACS tests.

To load these facilities in a pytest environment, set a `pytest_plugins`
variable in a conftest.py
(Reference https://docs.pytest.org/en/latest/writing_plugins.html#requiring-loading-plugins-in-a-test-module-or-conftest-file)

    pytest_plugins = "gmxapi.testsupport"

.. seealso:: https://docs.pytest.org/en/latest/plugins.html#findpluginname

.. todo:: Consider moving this to a separate optional package.
"""

import json
import logging
import os
import shutil
import tempfile
import warnings
from contextlib import contextmanager
from enum import Enum
from typing import Union

import pytest

mpi_status = 'Test requires mpi4py managing 2 MPI ranks.'
skip_mpi = False
try:
    from mpi4py import MPI

    if not MPI.Is_initialized():
        skip_mpi = True
        mpi_status += ' MPI is not initialized'
    elif MPI.COMM_WORLD.Get_size() < 2:
        skip_mpi = True
        mpi_status += ' MPI context is too small.'
except ImportError:
    skip_mpi = True
    mpi_status += ' mpi4py is not available.'


def pytest_configure(config):
    config.addinivalue_line("markers", "withmpi_only: test requires mpi4py managing 2 MPI ranks.")


def pytest_runtest_setup(item):
    # Handle the withmpi_only marker.
    for _ in item.iter_markers(name='withmpi_only'):
        if skip_mpi:
            pytest.skip(mpi_status)
        # The API uses iteration because markers may be duplicated, but we only
        # care about whether 'withmpi_only' occurs at all.
        break


def pytest_addoption(parser):
    """Add command-line user options for the pytest invocation."""
    parser.addoption(
        '--rm',
        action='store',
        default='always',
        choices=['always', 'never', 'success'],
        help='Remove temporary directories "always", "never", or on "success".'
    )
    parser.addoption(
        '--threads',
        type=int,
        help='Maximum number of threads per process per gmxapi session.'
    )


class RmOption(Enum):
    """Enumerate allowable values of the --rm option."""
    always = 'always'
    never = 'never'
    success = 'success'


@pytest.fixture(scope='session')
def remove_tempdir(request) -> RmOption:
    """pytest fixture to get access to the --rm CLI option."""
    arg = request.config.getoption('--rm')
    return RmOption(arg)


@pytest.fixture(scope='session')
def gmxconfig():
    from .commandline import _config
    config = _config()
    yield config


@pytest.fixture(scope='session')
def mdrun_kwargs(request, gmxconfig):
    """pytest fixture to provide a mdrun_kwargs dictionary for the mdrun ResourceManager.
    """
    from gmxapi.simulation.mdrun import ResourceManager as _ResourceManager
    if gmxconfig is None:
        raise RuntimeError('--threads argument requires a usable gmxconfig.json')
    arg = request.config.getoption('--threads')
    if arg is None:
        return {}
    mpi_type = gmxconfig['gmx_mpi_type']
    if mpi_type is not None and mpi_type == "tmpi":
        kwargs = {'threads': int(arg)}
    else:
        kwargs = {}
    # TODO: (#3718) Normalize the handling of run-time arguments.
    _ResourceManager.mdrun_kwargs = dict(**kwargs)
    return kwargs


@contextmanager
def scoped_chdir(dir):
    oldpath = os.getcwd()
    os.chdir(dir)
    try:
        yield dir
        # If the `with` block using scoped_chdir produces an exception, it will
        # be raised at this point in this function. We want the exception to
        # propagate out of the `with` block, but first we want to restore the
        # original working directory, so we skip `except` but provide a `finally`.
    finally:
        os.chdir(oldpath)


@contextmanager
def _cleandir(remove_tempdir: Union[str, RmOption]):
    """Context manager for a clean temporary working directory.

    Arguments:
        remove_tempdir (RmOption): whether to remove temporary directory "always",
                                   "never", or on "success"

    Raises:
        ValueError: if remove_tempdir value is not valid.

    The context manager will issue a warning for each temporary directory that
    is not removed.
    """
    if not isinstance(remove_tempdir, RmOption):
        remove_tempdir = RmOption(remove_tempdir)

    newpath = tempfile.mkdtemp()

    def remove():
        shutil.rmtree(newpath)

    def warn():
        warnings.warn('Temporary directory not removed: {}'.format(newpath))

    # Initialize callback function reference
    if remove_tempdir == RmOption.always:
        callback = remove
    else:
        callback = warn

    try:
        with scoped_chdir(newpath):
            yield newpath
        # If we get to this line, the `with` block using _cleandir did not throw.
        # Clean up the temporary directory unless the user specified `--rm never`.
        # I.e. If the user specified `--rm success`, then we need to toggle from `warn` to `remove`.
        if remove_tempdir != RmOption.never:
            callback = remove
    finally:
        callback()


@pytest.fixture
def cleandir(remove_tempdir: RmOption):
    """Provide a clean temporary working directory for a test.

    Example usage:

        import os
        import pytest

        @pytest.mark.usefixtures("cleandir")
        def test_cwd_starts_empty():
            assert os.listdir(os.getcwd()) == []
            with open("myfile", "w") as f:
                f.write("hello")

        def test_cwd_also_starts_empty(cleandir):
            assert os.listdir(os.getcwd()) == []
            assert os.path.abspath(os.getcwd()) == os.path.abspath(cleandir)
            with open("myfile", "w") as f:
                f.write("hello")

        @pytest.mark.usefixtures("cleandir")
        class TestDirectoryInit(object):
            def test_cwd_starts_empty(self):
                assert os.listdir(os.getcwd()) == []
                with open("myfile", "w") as f:
                    f.write("hello")

            def test_cwd_also_starts_empty(self):
                assert os.listdir(os.getcwd()) == []
                with open("myfile", "w") as f:
                    f.write("hello")

    Ref: https://docs.pytest.org/en/latest/fixture.html#using-fixtures-from-classes-modules-or-projects
    """
    with _cleandir(remove_tempdir) as newdir:
        yield newdir


@pytest.fixture(scope='session')
def gmxcli():
    from .commandline import cli_executable
    command = cli_executable()

    if command is None:
        message = "Tests need 'gmx' command line tool, but could not find it on the path."
        raise RuntimeError(message)
    try:
        assert os.access(command, os.X_OK)
    except Exception as E:
        raise RuntimeError('"{}" is not an executable gmx wrapper program'.format(command)) from E
    yield command
