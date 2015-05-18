#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

# - Find FFTW 3
# Find the native FFTW headers and libraries.
#
#  ${FFTW}_INCLUDE_DIRS - where to find FFTW headers
#  ${FFTW}_LIBRARIES    - List of libraries when using FFTW.
#  ${FFTW}_PKG          - The name of the pkg-config package needed
#  ${FFTW}_HAVE_SIMD    - True if FFTW was built with SIMD support
#  ${FFTW}_HAVE_AVX     - True if FFTW was built with AVX support
#  ${FFTW}_FOUND        - True if FFTW was found
#  where ${FFTW} is FFTW or FFTWF

list(LENGTH FFTW_FIND_COMPONENTS FFTW_NUM_COMPONENTS_WANTED)
if(${FFTW_NUM_COMPONENTS_WANTED} LESS 1)
  message(FATAL_ERROR "No FFTW component to search given")
elseif(${FFTW_NUM_COMPONENTS_WANTED} GREATER 1)
  message(FATAL_ERROR "We only support finding one FFTW component at the time, go and implement it ;-)")
elseif(${FFTW_FIND_COMPONENTS} MATCHES "^fftw(f)?$")
  if (NOT FFTW_FIND_VERSION OR FFTW_FIND_VERSION EQUAL 3) #find FFTW3 by default
    string(TOUPPER "${FFTW_FIND_COMPONENTS}" FFTW)
    string(REGEX REPLACE "fftw" "fftw3" ${FFTW}_PKG "${FFTW_FIND_COMPONENTS}")
    set(${FFTW}_FUNCTION_PREFIX "${FFTW_FIND_COMPONENTS}")
  else()
    message(FATAL_ERROR "We only support finding FFTW version 3, go and implement it ;-)")
  endif()
else()
  message(FATAL_ERROR "We do not support finding ${FFTW_FIND_COMPONENTS}, go and implement it ;-)")
endif()

find_package(PkgConfig QUIET)
if(NOT __pkg_config_checked_PC_${FFTW} OR NOT ${FFTW}_LIBRARY)
  pkg_check_modules(PC_${FFTW} "${${FFTW}_PKG}")
  if(NOT PC_${FFTW}_FOUND)
    message(STATUS "pkg-config could not detect ${${FFTW}_PKG}, trying generic detection")
  endif()
endif()

find_path(${FFTW}_INCLUDE_DIR "fftw3.h" HINTS ${PC_${FFTW}_INCLUDE_DIRS})
find_library(${FFTW}_LIBRARY NAMES "${${FFTW}_PKG}" HINTS ${PC_${FFTW}_LIBRARY_DIRS})

set(${FFTW}_LIBRARIES "${${FFTW}_LIBRARY}")
set(${FFTW}_INCLUDE_DIRS "${${FFTW}_INCLUDE_DIR}")

#better error message than find_package_handle_standard_args
if (${FFTW}_LIBRARY AND ${FFTW}_INCLUDE_DIR)
  set(${FFTW}_FOUND TRUE)
elseif (NOT ${FFTW}_LIBRARY)
  message("Could not find ${${FFTW}_PKG} library named lib${${FFTW}_PKG}, please specify its location in CMAKE_PREFIX_PATH or ${FFTW}_LIBRARY by hand (e.g. -D${FFTW}_LIBRARY='/path/to/lib${${FFTW}_PKG}.so')")
elseif (NOT ${FFTW}_INCLUDE_DIR)
  message("Could not the ${${FFTW}_PKG} header fftw3.h, please specify its path in ${FFTW}_INCLUDE_DIR by hand (e.g. -D${FFTW}_INCLUDE_DIR='/path/to/include')")
endif()

if (${FFTW}_FOUND)
  #The user could specify trash in ${FFTW}_LIBRARY, so test if we can link it
  include(CheckLibraryExists)
  include(gmxOptionUtilities)
  if (HAVE_LIBM)
    #adding MATH_LIBRARIES here to allow static libs, this does not harm us as we are anyway using it
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif ()
  gmx_check_if_changed(FFTW_LIBRARY_CHANGED ${FFTW}_LIBRARIES)
  if (FFTW_LIBRARY_CHANGED)
    unset(FOUND_${FFTW}_PLAN CACHE)
  endif()
  check_library_exists("${${FFTW}_LIBRARIES}" "${${FFTW}_FUNCTION_PREFIX}_plan_r2r_1d" "" FOUND_${FFTW}_PLAN)
  if(NOT FOUND_${FFTW}_PLAN)
    message(FATAL_ERROR "Could not find ${${FFTW}_FUNCTION_PREFIX}_plan_r2r_1d in ${${FFTW}_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what went wrong. If you are using a static lib (.a) make sure you have specified all dependencies of ${${FFTW}_PKG} in ${FFTW}_LIBRARY by hand (e.g. -D${FFTW}_LIBRARY='/path/to/lib${${FFTW}_PKG}.so;/path/to/libm.so') !")
  endif()

  # Check for FFTW3 compiled with --enable-sse
  foreach(SSE_FUNCTION ${${FFTW}_FUNCTION_PREFIX}_have_simd_sse)
    if (FFTW_LIBRARY_CHANGED)
      unset(${FFTW}_HAVE_${SSE_FUNCTION} CACHE)
    endif()
    check_library_exists("${${FFTW}_LIBRARIES}" "${SSE_FUNCTION}" "" ${FFTW}_HAVE_${SSE_FUNCTION})
    if(${FFTW}_HAVE_${SSE_FUNCTION})
      set(${FFTW}_HAVE_SSE TRUE)
      break()
    endif()
  endforeach()

  # Check for FFTW3 compiled with --enable-sse2
  foreach(SSE2_FUNCTION ${${FFTW}_FUNCTION_PREFIX}_have_simd_sse2)
    if (FFTW_LIBRARY_CHANGED)
      unset(${FFTW}_HAVE_${SSE2_FUNCTION} CACHE)
    endif()
    check_library_exists("${${FFTW}_LIBRARIES}" "${SSE2_FUNCTION}" "" ${FFTW}_HAVE_${SSE2_FUNCTION})
    if(${FFTW}_HAVE_${SSE2_FUNCTION})
      set(${FFTW}_HAVE_SSE2 TRUE)
      break()
    endif()
  endforeach()

  # Check for FFTW3 with 128-bit AVX compiled with --enable-avx
  foreach(AVX_128_FUNCTION ${${FFTW}_FUNCTION_PREFIX}_have_simd_avx_128)
    if (FFTW_LIBRARY_CHANGED)
      unset(${FFTW}_HAVE_${AVX_128_FUNCTION} CACHE)
    endif()
    check_library_exists("${${FFTW}_LIBRARIES}" "${AVX_128_FUNCTION}" "" ${FFTW}_HAVE_${AVX_128_FUNCTION})
    if(${FFTW}_HAVE_${AVX_128_FUNCTION})
      set(${FFTW}_HAVE_AVX_128 TRUE)
      break()
    endif()
  endforeach()

  # Check for FFTW3 with 128-bit AVX2 compiled with --enable-avx2
  foreach(AVX2_128_FUNCTION ${${FFTW}_FUNCTION_PREFIX}_have_simd_avx2_128)
    if (FFTW_LIBRARY_CHANGED)
      unset(${FFTW}_HAVE_${AVX2_128_FUNCTION} CACHE)
    endif()
    check_library_exists("${${FFTW}_LIBRARIES}" "${AVX2_128_FUNCTION}" "" ${FFTW}_HAVE_${AVX2_128_FUNCTION})
    if(${FFTW}_HAVE_${AVX2_128_FUNCTION})
      set(${FFTW}_HAVE_AVX2_128 TRUE)
      break()
    endif()
  endforeach()

  #in 3.3 sse function name has changed
  foreach(SIMD_FCT ${${FFTW}_FUNCTION_PREFIX}_have_simd_sse2;${${FFTW}_FUNCTION_PREFIX}_have_simd_avx;${${FFTW}_FUNCTION_PREFIX}_have_simd_altivec;${${FFTW}_FUNCTION_PREFIX}_have_simd_neon;${${FFTW}_FUNCTION_PREFIX}_have_sse2;${${FFTW}_FUNCTION_PREFIX}_have_sse;${${FFTW}_FUNCTION_PREFIX}_have_altivec)
    if (FFTW_LIBRARY_CHANGED)
      unset(${FFTW}_HAVE_${SIMD_FCT} CACHE)
    endif()
    check_library_exists("${${FFTW}_LIBRARIES}" "${SIMD_FCT}" "" ${FFTW}_HAVE_${SIMD_FCT})
    if(${FFTW}_HAVE_${SIMD_FCT})
      set(${FFTW}_HAVE_SIMD TRUE)
      break()
    endif()
  endforeach()
  #Verify FFTW is compiled with fPIC (necessary for shared libraries)
  if (CMAKE_OBJDUMP AND CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS AND NOT CYGWIN)
      execute_process(COMMAND ${CMAKE_OBJDUMP} --reloc ${${FFTW}_LIBRARY} OUTPUT_VARIABLE ${FFTW}_OBJDUMP)
      if (${${FFTW}_OBJDUMP} MATCHES "R_X86_64" #Should always be true for static libraries. Checks that objdump works properly and that the library isn't dynamic
              AND NOT ${${FFTW}_OBJDUMP} MATCHES "R_X86_64_PLT32")
          message(FATAL_ERROR "The FFTW library ${${FFTW}_LIBRARY} cannot be used with shared libraries. Provide a different FFTW library by setting ${FFTW}_LIBRARY. If you don't have a different one, recompile FFTW with \"--enable-shared\" or \"--with-pic\". Or disable shared libraries for GROMACS by setting BUILD_SHARED_LIBS to \"no\". Note: Disabling shared libraries requires up to 10x as much disk space.")
      endif()
  endif()
  set(CMAKE_REQUIRED_LIBRARIES)
endif ()

mark_as_advanced(${FFTW}_INCLUDE_DIR ${FFTW}_LIBRARY)
