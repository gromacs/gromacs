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

"""Provide the high-level interface to the file i/o behaviors the gmx package.

"""

# TODO: fix namespace polution
import gmx.core

class TrajectoryFile:
    """Provides an interface to Gromacs supported trajectory file formats.

    TrajectoryFile objects are Trajectory data sources or sinks (depending on
    access mode) that can be used in a Gromacs tool chain. If data is requested
    from an object before it has been explicitly attached to a tool chain or
    runner, these are implicitly created. If no execution context is specified,
    a local single-threaded context is created.

    Thus, stand-alone TrajectoryFile objects can serve as a rudimentary interface
    to file types supported by Gromacs.

    If the file mode is 'r', the object created supports the Python iterator
    protocol for reading one frame at a time.

    Other file modes are not yet supported.

    Example usage:

       import gmx

       # Create the Python proxy to the caching gmx::TrajectoryAnalysisModule object.
       mytraj = gmx.io.TrajectoryFile(filename, 'r')

       # Implicitly create the Runner object and get an iterator based on selection.
       frames = mytraj.select(...)

       # Iterating on the module advances the Runner.
       # Since the Python interpreter is explicitly asking for data,
       # the runner must now be initialized and begin execution.
       # mytraj.runner.initialize(context, options)
       # mytraj.runner.next()
       next(frames)

       # Subsequent iterations only need to step the runner and return a frame.
       for frame in frames:
           # do some stuff

       # The generator yielding frames has finished, so the runner has been released.
       # The module caching the frames still exists and could still be accessed or
       # given to a new runner with a new selection.
       for frame in mytraj.select(...):
           # do some other stuff

    """

    # TODO: decide how to represent file mode and how/if to reflect in the C++ API
    READ = 1

    def __init__(self, filename, mode=None):
        """Prepare filename for the requested mode of access.

        Args:
            filename (str): filesystem path to the file to be opened.
            mode (str): file access mode. 'r' for read.

        """
        if mode == 'r':
            self.mode = TrajectoryFile.READ
            self.filename = filename
            self._cpp_module = gmx.core.CachingTafModule()
        #elif mode == 'w':
        #    self.mode = TrajectoryFile.WRITE
        else:
            raise ValueError("Trajectory file access mode not supported.")

    def select(self, selection=None):
        """Generator to read atom positions frame-by-frame.

        An analysis module runner is implicitly created when the iterator
        is returned and destroyed when the iterator raises StopIteration.

        Args:
            selection: atom selection to retrieve from trajectory file

        Returns:
            iterator to frames

        """

        # Implementation details:
        #
        # 1. Python asks for first frame.
        # 2. Runner starts.
        # 3. Runner calls analyze_frame via next().
        # 4. Runner returns control to Python interpreter and yields frame.
        # 5. Python accesses the frame by interacting with the module.
        # 6. Python asks for next frame.
        # 7. Runner calls finish for the current frame and analyze for the next.
        # 8. Runner returns control to Python.
        # 9. When Runner runs out of frames, the Python generator is done.
        # 10. When Python leaves the context object or otherwise destroys the runner, it cleans up.

        #assert(isinstance(selection, gmx.core.Selection))

        # Create runner and bind module
        runner = gmx.core.TafRunner(self._cpp_module)

        # Create options object with which to initialize runner
        options = gmx.core.Options(filename=self.filename)

        # Initialize runner and module
        runner.initialize(options)
        while runner.next():
            frame = self._cpp_module.frame()
            yield frame
