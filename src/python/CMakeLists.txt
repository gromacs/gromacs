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

# Set the rpath of the
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

pybind11_add_module(core gmx_core.cpp data.cpp options.cpp runner.cpp trajectory.cpp)
target_link_libraries(core PRIVATE libgromacs)

# use, i.e. don't skip the full RPATH for the build tree
set_target_properties(core PROPERTIES SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set_target_properties(core PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)
set_target_properties(core PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
set_target_properties(core PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# Install target should copy Python source files into lib/python{version}/gmx
# and properly install core.so into same.

# Create the 'install' target, which is different than the setuptools
# installation location.
set(PYGMX_INSTALL_DIR ${LIB_INSTALL_DIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} CACHE PATH "Installation directory for python module")
install(DIRECTORY "gmx" DESTINATION ${PYGMX_INSTALL_DIR})
install(CODE "MESSAGE(\"Python module installed to ${CMAKE_INSTALL_PREFIX}/${PYGMX_INSTALL_DIR}.
You will need to update your PYTHONPATH environment variable for python to find it.\")")
# The compiled core module is a submodule of the 'gmx' package.
install(TARGETS core DESTINATION ${PYGMX_INSTALL_DIR}/gmx)

# TODO:
# Allow build and installation via setuptools from files configured by CMake.
#
# For user installation of a freshly built gmx module, we need to copy the right
# source files, configure setup.py, and provide a directory from which the user
# can run `pip install`.
# To do: we should probably check somehow that the version of pip matches the
# configured Python executable.

# Perform CMake variable substitution
#configure_file(
#    "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in"
#    "${CMAKE_BINARY_DIR}/setup.py"
#)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
