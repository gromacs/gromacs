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

@pytest.mark.skipif(not has_feature('fr15'),
                   reason="Feature level not met.")
def test_fr15():
    """FR15: Simulation input modification.

    * *gmx.modify_input produces new (tpr) simulation input in data flow operation*
      (requires interaction with library development)
    * gmx.make_input dispatches appropriate preprocessing for file or in-memory simulation input.
    """
    initial_input = gmx.read_tpr([tpr_filename for _ in range(10)])
    tau_t = list([i/10. for i in range(10)])
    param_sweep = gmx.modify_input(input=initial_input,
                                   parameters={
                                       'tau_t': tau_t
                                   }
                                   )
    md = gmx.mdrun(param_sweep)
    for tau_expected, tau_actual in zip(tau_t, md.output.params['tau_t'].extract()):
        assert tau_expected == tau_actual
