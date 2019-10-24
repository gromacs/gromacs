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

# Configure the `logging` module before proceeding any further.
gmx.logger.setLevel(logging.WARNING)

try:
    from mpi4py import MPI
    rank_number = MPI.COMM_WORLD.Get_rank()
except ImportError:
    rank_number = 0
    rank_tag = ''
    MPI = None
else:
    rank_tag = 'rank{}:'.format(rank_number)

# Use this formatter to improve the caplog log records.
formatter = logging.Formatter(rank_tag + '%(name)s:%(levelname)s: %(message)s')

# For additional console logging, create and attach a stream handler.
# For example:
#    ch = logging.StreamHandler()
#    ch.setFormatter(formatter)
#    logging.getLogger().addHandler(ch)


@pytest.mark.usefixtures('cleandir')
def test_run_from_tpr(spc_water_box):
    assert os.path.exists(spc_water_box)

    md = gmx.mdrun(spc_water_box)
    md.run()
    # TODO: better handling of output on unused MPI ranks.


@withmpi_only
@pytest.mark.usefixtures('cleandir')
def test_run_trivial_ensemble(spc_water_box, caplog):
    from mpi4py import MPI
    current_rank = MPI.COMM_WORLD.Get_rank()
    with caplog.at_level(logging.DEBUG):
        caplog.handler.setFormatter(formatter)
        with caplog.at_level(logging.WARNING, 'gmxapi'), \
                caplog.at_level(logging.DEBUG, 'gmxapi.mdrun'), \
                caplog.at_level(logging.DEBUG, 'gmxapi.modify_input'), \
                caplog.at_level(logging.DEBUG, 'gmxapi.read_tpr'), \
                caplog.at_level(logging.DEBUG, 'gmxapi.simulation'):

            tpr_filename = spc_water_box
            ensemble_width = 2
            simulation_input = gmx.read_tpr([tpr_filename] * ensemble_width)
            assert simulation_input.output.ensemble_width == ensemble_width
            assert len(simulation_input.output._simulation_input.result()) == ensemble_width
            md = gmx.mdrun(simulation_input)
            assert md.output.ensemble_width == ensemble_width
            md.run()

            output_directory = md.output._work_dir.result()
            logging.info('output_directory result: {}'.format(str(output_directory)))
            assert len(output_directory) == 2

            # Note that the 'cleandir' test fixture will clean up the output directory on
            # other ranks, so only check the current rank. Generally, our behavior
            # is undefined if the client removes the working directory while the job
            # is in progress. We can consider adding some sort of synchronization at
            # the end of the job if running in temporary directories becomes an
            # important use case outside of testing.
            assert output_directory[0] != output_directory[1]
            assert os.path.exists(output_directory[current_rank])
            assert os.path.exists(md.output.trajectory.result()[current_rank])


@pytest.mark.usefixtures('cleandir')
def test_run_from_read_tpr_op(spc_water_box, caplog):
    with caplog.at_level(logging.DEBUG):
        caplog.handler.setFormatter(formatter)
        with caplog.at_level(logging.DEBUG, 'gmxapi'):
            simulation_input = gmx.read_tpr(spc_water_box)
            md = gmx.mdrun(input=simulation_input)

            md.run()
            if rank_number == 0:
                assert os.path.exists(md.output.trajectory.result())


@pytest.mark.usefixtures('cleandir')
def test_run_from_modify_input_op(spc_water_box, caplog):
    with caplog.at_level(logging.DEBUG):

        simulation_input = gmx.read_tpr(spc_water_box)
        modified_input = gmx.modify_input(input=simulation_input, parameters={'nsteps': 4})
        md = gmx.mdrun(input=modified_input)

        md.run()
