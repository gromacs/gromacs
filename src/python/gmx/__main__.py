#!/usr/bin/env python
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

"""
Provides the gmx.__main__() method for the gmx package main module.
Defines behavior when module is invoked as a script or with
``python -m gmx``
"""

import sys
import numpy

from .io import TrajectoryFile

# Create the Python proxy to the caching gmx::TrajectoryAnalysisModule object.
filename = sys.argv[1]
mytraj = TrajectoryFile(filename, 'r')

# Implicitly create the Runner object and get an iterator based on selection.
#frames = mytraj.select(...)

# Iterating on the module advances the Runner.
# Since the Python interpreter is explicitly asking for data,
# the runner must now be initialized and begin execution.
# mytraj.runner.initialize(context, options)
# mytraj.runner.next()
frames = mytraj.select('real selections not yet implemented')
try:
    frame = next(frames)
    print(frame)
    print(frame.position)
    print(frame.position.extract())
    print("{} atoms in frame".format(frame.position.extract().N))
    print(numpy.array(frame.position.extract(), copy=False))
except StopIteration:
    print("no frames")

# Subsequent iterations only need to step the runner and return a frame.
for frame in frames:
    print(frame.position.extract())

# The generator yielding frames has finished, so the runner has been released.
# The module caching the frames still exists and could still be accessed or
# given to a new runner with a new selection.
#for frame in mytraj.select(...):
#   # do some other stuff
