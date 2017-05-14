#!/usr/bin/env/python
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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

"""Provide Python access to Gromacs

The gmx Python module provides an interface suitable for scripting Gromacs
workflows, interactive use, or connectivity to external Python-based APIs.

The API allows interaction with Gromacs that is decoupled from command-line
interfaces, terminal I/O, and filesystem access. Computation and data management
are managed by the API until/unless the user explicitly requests specific data,
such as for writing to a local file or manipulating with a tool that does not
implement the Gromacs API.

When data must be retrieved from Gromacs, efforts are made to do so as efficiently
as possible, so the user should consult documentation for the specific Gromacs
objects they are interested in regarding efficient access, if performance is
critical. For instance, exporting Gromacs data directly to a numpy array can be
much faster and require less memory than exporting to a Python list or tuple,
and using available iterators directly can save a lot of memory versus creating
an array and then iterating over it in two separate steps.

For more efficient iterative access to Gromacs data, such as analyzing a simulation
in progress, or applying Python analysis scripts in a trajectory analysis workflow,
consider using an appropriate call-back or, better yet, creating a C++ plugin
that can be inserted directly in the tool chain.

For more advanced use, the module provides means to access or manipulate Gromacs
more granularly than the command-line tool. This allows rapid prototyping of
new methods, debugging of unexpected simulation behavior, and adaptive workflows.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TODO: what should happen when ``from gmx import *``?
__all__ = ['System', 'exceptions']

import gmx
from gmx import util
from gmx import fileio
import gmx.core

class System(object):
    """Gromacs simulation system objects.

    A System object connects all of the objects necessary to describe a molecular
    system to be simulated.
    """

    @staticmethod
    def _from_file(inputrecord):
        """Process a file to create a System object.

        Calls an appropriate helper function to parse a file in the current context and create a System object. If no Context is currently bound, a local default context is created and bound.

        TODO: clarify the file location relative to the execution context.
        Until then, this helper method should not be part of the public
        interface.

        Args:
            inputrecord (str): input file name

        Returns:
            gmx.System
        """
        system = gmx.System()
        if util._filetype(inputrecord) is fileio.TprFile:
            system.runner = gmx.core._runner_from_tpr(inputrecord)
        return system

    def run(self):
        """Helper method to launch execution.

        If the System is attached to a Context, the Context is initialized, if
        necessary, and its instance run()
        method is called. If there is not yet a Context, one is created and then
        used.

        Returns:
            Gromacs status object.
        """
        return 0
