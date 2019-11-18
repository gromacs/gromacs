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

# Note that this file is heavily influenced by
# https://github.com/pybind/pybind11/blob/master/tests/CMakeLists.txt
# (Future editors would be advised to return to this source for inspiration.)
# The above copyright notice can be considered to apply to differences between
# the current file and the file at
# https://github.com/pybind/pybind11/blob/a3f4a0e8ab38f6b67174a162a5f0324e900a9c6b/tests/CMakeLists.txt
# Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>

set(GMXPY_TEST_FILES
    test_commandline.py
    test_exceptions.py
    test_mdrun.py
    test_operation.py
    )

# Make sure pytest is found or produce a fatal error
if(NOT GMXAPI_PYTEST_FOUND)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import pytest; print(pytest.__version__)"
                    RESULT_VARIABLE pytest_not_found OUTPUT_VARIABLE pytest_version ERROR_QUIET)
    if(pytest_not_found)
        message(FATAL_ERROR "Running the tests requires pytest. Please install it manually"
                " (try: ${PYTHON_EXECUTABLE} -m pip install pytest)")
    elseif(pytest_version VERSION_LESS 3.0)
        message(FATAL_ERROR "Running the tests requires pytest >= 3.0. Found: ${pytest_version}"
                "Please update it (try: ${PYTHON_EXECUTABLE} -m pip install -U pytest)")
    endif()
    set(GMXAPI_PYTEST_FOUND TRUE CACHE INTERNAL "Suppress checking for Python pytest module.")
endif()

# TODO: (#2961) Let the Python package discover required resources instead of passing GMXBIN
get_target_property(GMXBIN gmx RUNTIME_OUTPUT_DIRECTORY)
add_custom_target(gmxapi_pytest
                  COMMAND ${CMAKE_COMMAND} -E env GMXBIN=${GMXBIN}
                    ${PYTHON_EXECUTABLE} -m pytest --log-cli-level ERROR ${CMAKE_CURRENT_SOURCE_DIR}
                  DEPENDS _gmxapi
                  WORKING_DIRECTORY ${GMXAPI_PYTHON_STAGING_DIR})
# The current test fixtures require the `gmx` tool-wrapper executable.
add_dependencies(gmxapi_pytest gmx)

add_dependencies(check gmxapi_pytest)
