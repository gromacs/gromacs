#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019,2020, by the GROMACS development team, led by
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

"""Test gmxapi functionality described in roadmap.rst."""

import logging

import pytest

import gmxapi as gmx

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

@pytest.mark.usefixtures('cleandir')
def test_fr15(spc_water_box, caplog):
    """FR15: Simulation input modification.

    * *gmx.modify_input produces new (tpr) simulation input in data flow operation*
      (requires interaction with library development)
    """
    try:
        from mpi4py import MPI
        current_rank = MPI.COMM_WORLD.Get_rank()
        ensemble_size = MPI.COMM_WORLD.Get_size()
    except ImportError:
        current_rank = 0
        ensemble_size = 1

    with caplog.at_level(logging.DEBUG):
        caplog.handler.setFormatter(formatter)
        with caplog.at_level(logging.WARNING, 'gmxapi'), \
                 caplog.at_level(logging.DEBUG, 'gmxapi.mdrun'), \
                 caplog.at_level(logging.DEBUG, 'gmxapi.modify_input'), \
                 caplog.at_level(logging.DEBUG, 'gmxapi.read_tpr'), \
                 caplog.at_level(logging.DEBUG, 'gmxapi.simulation'):

            initial_input = gmx.read_tpr([spc_water_box for _ in range(ensemble_size)])
            parameters = list([{'ld-seed': i} for i in range(ensemble_size)])
            param_sweep = gmx.modify_input(input=initial_input,
                                           parameters=parameters
                                           )
            md = gmx.mdrun(param_sweep)
            # TODO: (#3179) Handle ensembles of size 1 more elegantly.
            if md.output.ensemble_width > 1:
                result_list = md.output.parameters['ld-seed'].result()
            else:
                result_list = [md.output.parameters['ld-seed'].result()]
            for expected, actual in zip(parameters, result_list):
                assert expected['ld-seed'] == actual
