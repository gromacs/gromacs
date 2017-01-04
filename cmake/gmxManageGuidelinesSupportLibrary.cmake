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

# GROMACS used to have an optional dependency on the GNU Scientific
# Library, which shares the initialism "GSL," so this option is named
# fully.
option(GMX_EXTERNAL_GUIDELINES_SUPPORT_LIBRARY "Use external Microsoft Guidelines Support Library instead of compiling the version bundled with GROMACS." OFF)
mark_as_advanced(GMX_EXTERNAL_GUIDELINES_SUPPORT_LIBRARY)
if(GMX_EXTERNAL_GUIDELINES_SUPPORT_LIBRARY)
    check_cxx_source_compiles("#include <gsl/string_span>\nint main(){ gsl::cstring_span<> s;}" HAVE_GUIDELINES_SUPPORT_LIBRARY)
    if(NOT HAVE_GUIDELINES_SUPPORT_LIBRARY)
        message(FATAL_ERROR "An external Microsoft Guidelines Support Library could not be found in the CMake search path, please adjust")
    endif()
else()
    # Use the version bundled with GROMACS
    include_directories(BEFORE SYSTEM "src/external/gsl")
endif()

foreach(_build_type DEBUG RELWITHASSERT)
    set(CMAKE_C_FLAGS_${_build_type} "${CMAKE_C_FLAGS_${_build_type}} -DGSL_THROW_ON_CONTRACT_VIOLATION")
    set(CMAKE_CXX_FLAGS_${_build_type} "${CMAKE_CXX_FLAGS_${_build_type}} -DGSL_THROW_ON_CONTRACT_VIOLATION")
endforeach()
