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

"""gmxapi Python package for GROMACS.

This package provides Python access to GROMACS molecular simulation tools.
Operations can be connected flexibly to allow high performance simulation and
analysis with complex control and data flows. Users can define new operations
in C++ or Python with the same tool kit used to implement this package.

"""

__all__ = ['commandline_operation',
           'concatenate_lists',
           'function_wrapper',
           'join_arrays',
           'logger',
           'logical_not',
           'make_constant',
           'mdrun',
           'modify_input',
           'ndarray',
           'read_tpr',
           'subgraph',
           'while_loop',
           'NDArray',
           '__version__']

from ._logging import logger
from .version import __version__

# Import utilities
from .operation import computed_result, function_wrapper
# Import public types
from .datamodel import NDArray
# Import the public operations
from .datamodel import ndarray
from .operation import concatenate_lists, join_arrays, logical_not, make_constant
from .commandline import commandline_operation
from .simulation import mdrun, modify_input, read_tpr
# TODO: decide where this lives
from .operation import subgraph
# TODO: decide where this lives
from .operation import while_loop
