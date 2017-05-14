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

# There is not yet a build target for running these Python unit tests.
# For the moment, running them directly works fine with
# ``PYTHONPATH=/path/to/module python -m unittest test_something``.
# Note that the trailing ``py`` is removed because the name is interpreted
# as a module. We can rely on sphinx.ext.doctest for a lot, but I
# think it is more important to have Python tests work from distutils
# than from CMake, but I'm deferring the topic for now.
#
# We can do some testing on the C++ code for the Python bindings, but
# we might want to just test the bindings from Python.

#gmx_add_unit_test(PythonApiUnitTests trajectoryanalysis-test
#
#                  $<TARGET_OBJECTS:analysisdata-test-shared>)

#add_executable(test_pywrappers ${UNITTEST_TARGET_OPTIONS} test_selection.cpp)
#target_link_libraries(test_pywrappers core ${GMX_EXE_LINKER_FLAGS} ${GMX_STDLIB_LIBRARIES})
