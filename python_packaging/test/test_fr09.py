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

"""Test gmxapi functionality described in roadmap.rst."""

import pytest

import gmxapi as gmx
from gmxapi.version import has_feature

@pytest.mark.skipif(not has_feature('fr9'),
                   reason="Feature level not met.")
def test_fr9():
    """FR9: *gmx.mdrun supports interface for binding MD plugins*

    (requires interaction with library development)
    """
    import sample_restraint

    starting_structure = 'input_conf.gro'
    topology_file = 'input.top'
    run_parameters = 'params.mdp'

    initial_tpr = gmx.commandline_operation(
        'gmx',
        'grompp',
        input={
            '-f': run_parameters,
            '-c': starting_structure,
            '-p': topology_file
        },
        output={'-o': gmx.OutputFile('.tpr')})

    simulation_input = gmx.read_tpr(initial_tpr.output.file['-o'])

    # Prepare a simple harmonic restraint between atoms 1 and 4
    restraint_params = {'sites': [1, 4],
                        'R0': 2.0,
                        'k': 10000.0}

    restraint = sample_restraint.harmonic_restraint(input=restraint_params)

    md = gmx.mdrun(input=simulation_input, potential=sample_restraint)

    #md.run()
