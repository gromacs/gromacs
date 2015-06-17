#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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

# - Define macro to check if isfinite or _isfinite exists
#
#  gmx_test_isfinite(VARIABLE)
#
#  VARIABLE will be set to true if isfinite exists
#
#  gmx_test__isfinite(VARIABLE)
#
#  VARIABLE will be set to true if _isfinite exists
#
#  gmx_test__finite(VARIABLE) - disabled since it doesn't seem to work the way the MSVC docs suggest
#
#  VARIABLE will be set to true if _finite exists
#

MACRO(gmx_test_isfinite VARIABLE)

  if(NOT DEFINED isfinite_compile_ok)
    MESSAGE(STATUS "Checking for isfinite")

    set(CMAKE_REQUIRED_INCLUDES "math.h")
    if (HAVE_LIBM)
        set(CMAKE_REQUIRED_LIBRARIES "m")
    endif()
    check_c_source_compiles(
      "#include <math.h>
int main(void) {
  float f;
  isfinite(f);
  return 0;
}" isfinite_compile_ok)

    if(isfinite_compile_ok)
        MESSAGE(STATUS "Checking for isfinite - yes")
    else()
        MESSAGE(STATUS "Checking for isfinite - no")
    endif()
    set(isfinite_compile_ok "${isfinite_compile_ok}" CACHE INTERNAL "Result of isfinite check")
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()

  if(isfinite_compile_ok)
    set(${VARIABLE} ${isfinite_compile_ok}
                "Result of test for isfinite")
  endif()

ENDMACRO(gmx_test_isfinite VARIABLE)

MACRO(gmx_test__isfinite VARIABLE)

  if(NOT DEFINED _isfinite_compile_ok)
    MESSAGE(STATUS "Checking for _isfinite")

    set(CMAKE_REQUIRED_INCLUDES "math.h")
    if (HAVE_LIBM)
        set(CMAKE_REQUIRED_LIBRARIES "m")
    endif()
    check_c_source_compiles(
      "#include <math.h>
int main(void) {
  float f;
  _isfinite(f);
}" _isfinite_compile_ok)

    if(_isfinite_compile_ok)
        MESSAGE(STATUS "Checking for _isfinite - yes")
    else()
        MESSAGE(STATUS "Checking for _isfinite - no")
    endif()
    set(_isfinite_compile_ok "${_isfinite_compile_ok}" CACHE INTERNAL "Result of _isfinite check")
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()

  if(_isfinite_compile_ok)
    set(${VARIABLE} ${_isfinite_compile_ok}
                "Result of test for _isfinite")
  endif()

ENDMACRO(gmx_test__isfinite VARIABLE)

# Necessary for MSVC
MACRO(gmx_test__finite VARIABLE)

  if(NOT DEFINED _finite_compile_ok)
    MESSAGE(STATUS "Checking for _finite")

    set(CMAKE_REQUIRED_INCLUDES "float.h")
    check_c_source_compiles(
      "#include <float.h>
int main(void) {
  float f;
  _finite(f);
}" _finite_compile_ok)

    if(_finite_compile_ok)
        MESSAGE(STATUS "Checking for _finite - yes")
    else()
        MESSAGE(STATUS "Checking for _finite - no")
    endif()
    set(_finite_compile_ok "${_finite_compile_ok}" CACHE INTERNAL "Result of _finite check")
    set(CMAKE_REQUIRED_INCLUDES)
  endif()

  if(_finite_compile_ok)
    set(${VARIABLE} ${_finite_compile_ok}
                "Result of test for _finite")
  endif()

ENDMACRO(gmx_test__finite VARIABLE)
