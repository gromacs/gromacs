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

"""Test gromacs.mdrun operation.

Factory produces deferred execution operation.

TODO: Factory expects input for conformation, topology, simulation parameters, simulation state.

TODO: Factory accepts additional keyword input to indicate binding
 to the "potential" interface.
"""

import logging
import os
import pytest

import gmxapi as gmx
from gmxapi.testsupport import withmpi_only

gmx.logger.setLevel(logging.WARNING)

# Configure the `logging` module before and non-built-in packages start to use it.
logging.getLogger().setLevel(logging.WARNING)
# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)
# create formatter and add it to the handler
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s: %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
logging.getLogger().addHandler(ch)


@pytest.mark.usefixtures('cleandir')
def test_run_from_tpr(spc_water_box):
    assert os.path.exists(spc_water_box)

    md = gmx.mdrun(spc_water_box)
    md.run()
