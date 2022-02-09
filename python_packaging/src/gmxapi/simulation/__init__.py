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

"""GROMACS simulation subpackage for gmxapi.

Provides operations for configuring and running molecular simulations.

The initial version of this module is a port of the gmxapi 0.0.7 facilities from
https://github.com/kassonlab/gmxapi and is not completely integrated with the
gmxapi 0.1 specification. Operation execution is dispatched to the old execution
manager for effective ensemble handling and C++ MD module binding. This should
be an implementation detail that is not apparent to the typical user, but it is
worth noting that chains of gmxapi.simulation module operations will be
automatically bundled for execution as gmxapi 0.0.7 style API sessions. Run time
options and file handling will necessarily change as gmxapi data flow handling
evolves.

In other words, if you rely on behavior not specified explicitly in the user
documentation, please keep an eye on the module documentation when updating
gmxapi and please participate in the ongoing discussions for design and
implementation.
"""

__all__ = ['abc',
           'mdrun',
           'modify_input',
           'read_tpr']

from gmxapi.simulation import abc
from .mdrun import mdrun
from .read_tpr import read_tpr
from .modify_input import modify_input
