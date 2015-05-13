#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015, by the GROMACS development team, led by
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

# Adapted from code posted on cmake-users by Mark Moll
find_package(PythonInterp)
function(find_python_module module)
    string(TOUPPER ${module} module_upper)
    if(NOT PYTHONMODULE_${module_upper})
        if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
	    set(${module}_FIND_REQUIRED TRUE)
	endif()
        if (NOT PYTHON_EXECUTABLE)
            message(STATUS "Cannot find python module ${module} because no python executable is known")
        else()
	    # A module's location is usually a directory, but for binary modules
	    # it's a .so file.
	    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
	        "import re, ${module}; print re.compile('/__init__.py.*').sub('',${module}.__file__)"
	        RESULT_VARIABLE _${module}_status 
	        OUTPUT_VARIABLE _${module}_location
	        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()
	if(NOT _${module}_status)
	    set(PYTHONMODULE_${module_upper} ${_${module}_location} CACHE STRING 
		"Location of Python module ${module}")
	endif()
    endif()
    find_package_handle_standard_args(PYTHONMODULE_${module} DEFAULT_MSG PYTHONMODULE_${module_upper})
endfunction()
