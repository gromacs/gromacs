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

"""Reusable definitions for test modules.

Provides utilities and pytest fixtures for gmxapi and GROMACS tests.

To load these facilities in a pytest environment, set a `pytest_plugins`
variable in a conftest.py
(Reference https://docs.pytest.org/en/latest/writing_plugins.html#requiring-loading-plugins-in-a-test-module-or-conftest-file)

    pytest_plugins = "gmxapi.testsupport"

.. seealso:: https://docs.pytest.org/en/latest/plugins.html#findpluginname

.. todo:: Consider moving this to a separate optional package.
"""

import os
import warnings
from contextlib import contextmanager

import pytest

mpi_status = "Test requires mpi4py managing 2 MPI ranks."
skip_mpi = False
rank_number = 0
comm_size = 1
rank_tag = ""
comm = None
try:
    from mpi4py import MPI

    if not MPI.Is_initialized():
        skip_mpi = True
        mpi_status += " MPI is not initialized"
    else:
        comm = MPI.COMM_WORLD
        if comm.Get_size() < 2:
            skip_mpi = True
            mpi_status += " MPI context is too small."
        else:
            rank_number = comm.Get_rank()
            comm_size = comm.Get_size()
except ImportError:
    MPI = None
    skip_mpi = True
    mpi_status += " mpi4py is not available."


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "withmpi_only: test requires mpi4py managing 2 MPI ranks."
    )


def pytest_runtest_setup(item):
    # Handle the withmpi_only marker.
    for _ in item.iter_markers(name="withmpi_only"):
        if skip_mpi:
            pytest.skip(mpi_status)
        # The API uses iteration because markers may be duplicated, but we only
        # care about whether 'withmpi_only' occurs at all.
        break


def pytest_addoption(parser):
    """Add command-line user options for the pytest invocation."""
    # TODO(#4345): Remove `--rm` option after 2023 release.
    parser.addoption(
        "--rm",
        action="store",
        default="always",
        choices=["always", "never", "success"],
        help="No longer used. See https://docs.pytest.org/en/latest/how-to/tmp_path.html",
    )
    parser.addoption(
        "--threads",
        type=int,
        help="Maximum number of threads per process per gmxapi session.",
    )
    parser.addoption(
        "--pydevd",
        type=str,
        action="store",
        nargs="?",
        const="pydevd_pycharm",
        help="Attempt to connect to PyDev.Debugger using the indicated module.",
    )
    parser.addoption(
        "--pydevd-host", type=str, default="localhost", help="Set the pydevd host."
    )
    parser.addoption(
        "--pydevd-port", type=int, default=12345, help="Set the pydevd port."
    )
    parser.addoption(
        "--pydevd-rank",
        type=int,
        default=0,
        help="In MPI, specify the rank to attach to the debug server.",
    )


@pytest.fixture(scope="session", autouse=True)
def pydev_debug(request):
    """If requested, try to connect to a PyDev.Debugger backend at
    host.docker.internal:12345.

    Note: the IDE run configuration must be started before launching pytest.
    """
    pydevd_module = request.config.getoption("--pydevd")
    if pydevd_module:
        host = request.config.getoption("--pydevd-host")
        port = request.config.getoption("--pydevd-port")
        participating_rank = request.config.getoption("--pydevd-rank")
        if rank_number == participating_rank:
            try:
                import importlib

                pydevd = importlib.import_module(pydevd_module)
                settrace = getattr(pydevd, "settrace", None)
                if settrace is None:
                    raise RuntimeError(
                        f"PyDevD interface not supported. {pydevd} does not have `settrace`."
                    )
                return settrace(
                    host, port=port, stdoutToServer=True, stderrToServer=True
                )
            except ImportError:
                warnings.warn(f"{pydevd_module} not found. Ignoring `--pydevd` option.")


@pytest.fixture(scope="session")
def gmxconfig():
    from .commandline import _config

    config = _config()
    yield config


@pytest.fixture(scope="function")
def mdrun_kwargs(request, gmxconfig):
    """pytest fixture to provide a mdrun_kwargs dictionary for the mdrun ResourceManager.

    .. versionchanged:: 0.4.1
        Use function scope to avoid unexpectedly mutating the dictionary between tests.

    """
    if gmxconfig is None:
        raise RuntimeError("--threads argument requires a usable gmxconfig.json")
    arg = request.config.getoption("--threads")
    if arg is None:
        return {}
    mpi_type = gmxconfig["gmx_mpi_type"]
    if mpi_type is not None and mpi_type == "tmpi":
        kwargs = {"threads": int(arg)}
    else:
        kwargs = {}
    # TODO: (#3718) Normalize the handling of run-time arguments.
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


@pytest.fixture
def cleandir(tmp_path):
    """Provide a clean temporary working directory for a test.

    Temporary directory is created and managed as described at
    https://docs.pytest.org/en/latest/how-to/tmp_path.html

    TODO:
        Avoid or suppress ``PytestWarning: (rm_rf) error removing ...``
        from duplicated temporary directory management on multiple MPI ranks.
    """
    with scoped_chdir(tmp_path):
        yield tmp_path


@pytest.fixture(scope="session")
def gmxcli():
    from .commandline import cli_executable

    command = cli_executable()

    if command is None:
        message = (
            "Tests need 'gmx' command line tool, but could not find it on the path."
        )
        raise RuntimeError(message)
    try:
        assert os.access(command, os.X_OK)
    except Exception as E:
        raise RuntimeError(
            '"{}" is not an executable gmx wrapper program'.format(command)
        ) from E
    yield command
