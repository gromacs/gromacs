#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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
#
# - Finds OpenMP support
# Copied from cmake 2.8.7 prior versions were assuming C++
# compiler is enabled. Can be removed when 2.8.7 becomes 
# required or C++ becomes required.
#
# This module can be used to detect OpenMP support in a compiler.
# If the compiler supports OpenMP, the flags required to compile with
# openmp support are set.
#
# The following variables are set:
#   OpenMP_C_FLAGS - flags to add to the C compiler for OpenMP support
#   OpenMP_CXX_FLAGS - flags to add to the CXX compiler for OpenMP support
#   OPENMP_FOUND - true if openmp is detected
#
# Supported compilers can be found at http://openmp.org/wp/openmp-compilers/

#=============================================================================
# Copyright 2009 Kitware, Inc.
# Copyright 2008-2009 André Rigland Brodtkorb <Andre.Brodtkorb@ifi.uio.no>
# Copyright 2012 Rolf Eike Beer <eike@sf-mail.de>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

set(_OPENMP_REQUIRED_VARS)

function(_OPENMP_FLAG_CANDIDATES LANG)
  set(OpenMP_FLAG_CANDIDATES
    #GNU
    "-fopenmp"
    #Microsoft Visual Studio
    "/openmp"
    #Intel windows
    "-Qopenmp"
    #PathScale, Intel
    "-openmp"
    #Empty, if compiler automatically accepts openmp
    " "
    #Sun
    "-xopenmp"
    #HP
    "+Oopenmp"
    #IBM XL C/c++
    "-qsmp"
    #Portland Group, MIPSpro
    "-mp"
  )

  set(OMP_FLAG_GNU "-fopenmp")
  set(OMP_FLAG_HP "+Oopenmp")
  if(WIN32)
    set(OMP_FLAG_Intel "-Qopenmp")
  else()
    set(OMP_FLAG_Intel "-openmp")
  endif()
  set(OMP_FLAG_MIPSpro "-mp")
  set(OMP_FLAG_MSVC "/openmp")
  set(OMP_FLAG_PathScale "-openmp")
  set(OMP_FLAG_PGI "-mp")
  set(OMP_FLAG_SunPro "-xopenmp")
  set(OMP_FLAG_XL "-qsmp")

  # Move the flag that matches the compiler to the head of the list,
  # this is faster and doesn't clutter the output that much. If that
  # flag doesn't work we will still try all.
  if(OMP_FLAG_${CMAKE_${LANG}_COMPILER_ID})
    list(REMOVE_ITEM OpenMP_FLAG_CANDIDATES "${OMP_FLAG_${CMAKE_${LANG}_COMPILER_ID}}")
    list(INSERT OpenMP_FLAG_CANDIDATES 0 "${OMP_FLAG_${CMAKE_${LANG}_COMPILER_ID}}")
  endif()

  set(OpenMP_${LANG}_FLAG_CANDIDATES "${OpenMP_FLAG_CANDIDATES}" PARENT_SCOPE)
endfunction(_OPENMP_FLAG_CANDIDATES)

# sample openmp source code to test
set(OpenMP_C_TEST_SOURCE 
"
#include <omp.h>
int main() { 
#ifdef _OPENMP
  return 0; 
#else
  breaks_on_purpose
#endif
}
")

# check c compiler
if(CMAKE_C_COMPILER_LOADED)
  # if these are set then do not try to find them again,
  # by avoiding any try_compiles for the flags
  if(OpenMP_C_FLAGS)
    unset(OpenMP_C_FLAG_CANDIDATES)
  else()
    _OPENMP_FLAG_CANDIDATES("C")
    include(CheckCSourceCompiles)
  endif()

  foreach(FLAG ${OpenMP_C_FLAG_CANDIDATES})
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(OpenMP_FLAG_DETECTED CACHE)
    message(STATUS "Try OpenMP C flag = [${FLAG}]")
    check_c_source_compiles("${OpenMP_C_TEST_SOURCE}" OpenMP_FLAG_DETECTED)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(OpenMP_FLAG_DETECTED)
      set(OpenMP_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif(OpenMP_FLAG_DETECTED)
  endforeach(FLAG ${OpenMP_C_FLAG_CANDIDATES})

  set(OpenMP_C_FLAGS "${OpenMP_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for OpenMP parallization")

  list(APPEND _OPENMP_REQUIRED_VARS OpenMP_C_FLAGS)
  unset(OpenMP_C_FLAG_CANDIDATES)
endif()

# check cxx compiler
if(CMAKE_CXX_COMPILER_LOADED)
  # if these are set then do not try to find them again,
  # by avoiding any try_compiles for the flags
  if(OpenMP_CXX_FLAGS)
    unset(OpenMP_CXX_FLAG_CANDIDATES)
  else()
    _OPENMP_FLAG_CANDIDATES("CXX")
    include(CheckCXXSourceCompiles)

    # use the same source for CXX as C for now
    set(OpenMP_CXX_TEST_SOURCE ${OpenMP_C_TEST_SOURCE})
  endif()

  foreach(FLAG ${OpenMP_CXX_FLAG_CANDIDATES})
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(OpenMP_FLAG_DETECTED CACHE)
    message(STATUS "Try OpenMP CXX flag = [${FLAG}]")
    check_cxx_source_compiles("${OpenMP_CXX_TEST_SOURCE}" OpenMP_FLAG_DETECTED)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(OpenMP_FLAG_DETECTED)
      set(OpenMP_CXX_FLAGS_INTERNAL "${FLAG}")
      break()
    endif(OpenMP_FLAG_DETECTED)
  endforeach(FLAG ${OpenMP_CXX_FLAG_CANDIDATES})

  set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS_INTERNAL}"
    CACHE STRING "C++ compiler flags for OpenMP parallization")

  list(APPEND _OPENMP_REQUIRED_VARS OpenMP_CXX_FLAGS)
  unset(OpenMP_CXX_FLAG_CANDIDATES)
  unset(OpenMP_CXX_TEST_SOURCE)
endif()

if(_OPENMP_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(OpenMP
                                    REQUIRED_VARS ${_OPENMP_REQUIRED_VARS})

  mark_as_advanced(${_OPENMP_REQUIRED_VARS})

  unset(_OPENMP_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindOpenMP requires C or CXX language to be enabled")
endif()
