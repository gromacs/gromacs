# - Find FFTW 3
# Find the native FFTW headers and libraries.
#
#  ${FFTW}_INCLUDE_DIRS - where to find FFTW headers
#  ${FFTW}_LIBRARIES    - List of libraries when using FFTW.
#  ${FFTW}_PKG          - The name of the pkg-config package needed
#  ${FFTW}_HAVE_SIMD    - True if FFTW was build with SIMD support
#  ${FFTW}_FOUND        - True if FFTW was found
#  where ${FFTW} is FFTW or FFTWF
#
#  This file is part of Gromacs        Copyright (c) 2012

#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.

#  To help us fund GROMACS development, we humbly ask that you cite
#  the research papers on the package. Check out http://www.gromacs.org

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

find_package(PkgConfig)
if(NOT __pkg_config_checked_PC_${FFTW})
  pkg_check_modules(PC_${FFTW} "${${FFTW}_PKG}")
endif(NOT __pkg_config_checked_PC_${FFTW})

find_path(${FFTW}_INCLUDE_DIR "fftw3.h" HINTS ${PC_${FFTW}_INCLUDE_DIRS})
find_library(${FFTW}_LIBRARY NAMES "${${FFTW}_PKG}" HINTS ${PC_${FFTW}_LIBRARY_DIRS})

set(${FFTW}_LIBRARIES "${${FFTW}_LIBRARY}")
set(${FFTW}_INCLUDE_DIRS "${${FFTW}_INCLUDE_DIR}")

#better error message than find_package_handle_standard_args
if (${FFTW}_LIBRARY AND ${FFTW}_INCLUDE_DIR)
  set(${FFTW}_FOUND TRUE)
elseif (NOT ${FFTW}_LIBRARY)
  message("Could not find ${${FFTW}_PKG} library named lib${${FFTW}_PKG}, please specify its location in ${FFTW}_LIBRARY by hand (e.g. -D${FFTW}_LIBRARY='/path/to/lib${${FFTW}_PKG}.so')")
elseif (NOT ${FFTW}_INCLUDE_DIR)
  message("Could not the ${${FFTW}_PKG} header fftw3.h, please specify its path in ${FFTW}_INCLUDE_DIR by hand (e.g. -D${FFTW}_INCLUDE_DIR='/path/to/include')")
endif()

if (${FFTW}_FOUND)
  #The user could specify trash in ${FFTW}_LIBRARY, so test if we can link it
  include(CheckLibraryExists)
  if (HAVE_LIBM)
    #adding MATH_LIBRARIES here to allow static libs, this does not harm us as we are anyway using it
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif (HAVE_LIBM)
  check_library_exists("${${FFTW}_LIBRARIES}" "${${FFTW}_FUNCTION_PREFIX}_plan_r2r_1d" "" FOUND_${FFTW}_PLAN)
  if(NOT FOUND_${FFTW}_PLAN)
    message(FATAL_ERROR "Could not find ${${FFTW}_FUNCTION_PREFIX}_plan_r2r_1d in ${${FFTW}_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what went wrong. If you are using a static lib (.a) make sure you have specified all dependencies of ${${FFTW}_PKG} in ${FFTW}_LIBRARY by hand (e.g. -D${FFTW}_LIBRARY='/path/to/lib${${FFTW}_PKG}.so;/path/to/libm.so') !")
  endif(NOT FOUND_${FFTW}_PLAN)
  #in 3.3 sse function name has changed
  foreach(SIMD_FCT ${${FFTW}_FUNCTION_PREFIX}_have_simd_sse2;${${FFTW}_FUNCTION_PREFIX}_have_simd_avx;${${FFTW}_FUNCTION_PREFIX}_have_simd_altivec;${${FFTW}_FUNCTION_PREFIX}_have_simd_neon;${${FFTW}_FUNCTION_PREFIX}_have_sse2;${${FFTW}_FUNCTION_PREFIX}_have_sse;${${FFTW}_FUNCTION_PREFIX}_have_altivec)
    check_library_exists("${${FFTW}_LIBRARIES}" "${SIMD_FCT}" "" ${FFTW}_HAVE_${SIMD_FCT})
    if(${FFTW}_HAVE_${SIMD_FCT})
      set(${FFTW}_HAVE_SIMD TRUE CACHE  BOOL "If ${${FFTW}_PKG} was built with SIMD support")
      break()
    endif(${FFTW}_HAVE_${SIMD_FCT})
  endforeach()
  set(CMAKE_REQUIRED_LIBRARIES)
endif (${FFTW}_FOUND)
set(${FFTW}_HAVE_SIMD FALSE CACHE BOOL "If ${${FFTW}_PKG} was built with SIMD support")

mark_as_advanced(${FFTW}_INCLUDE_DIR ${FFTW}_LIBRARY ${FFTW}_HAVE_SIMD)
