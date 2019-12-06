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

"""Reusable definitions for test modules.

Provides utilities and pytest fixtures for gmxapi and GROMACS tests.

To load these facilities in a pytest environment, set a `pytest_plugins`
variable in a conftest.py
(Reference https://docs.pytest.org/en/latest/writing_plugins.html#requiring-loading-plugins-in-a-test-module-or-conftest-file)

    pytest_plugins = "gmxapi.testsupport"

.. seealso:: https://docs.pytest.org/en/latest/plugins.html#findpluginname

.. todo:: Consider moving this to a separate optional package.
"""

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
