# - Find FFTW 2/3 single double
# Find the native FFTW headers and libraries.
#
#  FFTW_INCLUDE_DIRS - where to find FFTW headers
#  FFTW_LIBRARIES    - List of libraries when using FFTW.
#  FFTW_PKG          - The name of the pkg-config package needed
#  FFTW_HAVE_SSE     - True if FFTW was build with SSE support
#  FFTW_FOUND        - True if FFTW was found
#
#  This file is part of Gromacs        Copyright (c) 2012

#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.

#  To help us fund GROMACS development, we humbly ask that you cite
#  the research papers on the package. Check out http://www.gromacs.org

find_package(PkgConfig)
list(LENGTH FFTW_FIND_COMPONENTS FFTW_NUM_COMPONENTS_WANTED)
if(${FFTW_NUM_COMPONENTS_WANTED} LESS 1)
  message(FATAL_ERROR "No FFTW component to search given")
elseif(${FFTW_NUM_COMPONENTS_WANTED} GREATER 1)
  message(FATAL_ERROR "We only support finding one FFTW component at the time, go and implement it ;-)")
elseif(${FFTW_FIND_COMPONENTS} MATCHES "^fftw(f)?$")
  if (${FFTW_FIND_VERSION} EQUAL 3 OR NOT FFTW_FIND_VERSION) #find FFTW3 by default
    string(REGEX REPLACE "fftw" "fftw3" FFTW_PKG "${FFTW_FIND_COMPONENTS}")
    set(FFTW_HEADER "fftw3.h")
    set(FFTW_FUNCTION "${FFTW_FIND_COMPONENTS}_plan_r2r_1d")
  #elseif(${FFTW_FIND_VERSION} EQUAL 2)
  #  set(FFTW_PKG "${FFTW_FIND_COMPONENTS}")
  #  set(FFTW_HEADER "${FFTW_FIND_COMPONENTS}.h")
  else()
    message(FATAL_ERROR "We only support finding FFTW version 3, go and implement it ;-)")
  endif()
else()
  message(FATAL_ERROR "We do not support finding ${FFTW_FIND_COMPONENTS}, go and implement it ;-)")
endif()

if(NOT __pkg_config_checked_PC_FFTW_${FFTW_PKG})
  pkg_check_modules(PC_FFTW_${FFTW_PKG} "${FFTW_PKG}")
endif(NOT __pkg_config_checked_PC_FFTW_${FFTW_PKG})

if(FFTW_LIBRARY)
  set(FFTW_LIBRARY_${FFTW_PKG} "${FFTW_LIBRARY}" CACHE INTERNAL "Path to ${FFTW_PKG} library" FORCE)
endif(FFTW_LIBRARY)

if(FFTW_INCLUDE_DIR)
  set(FFTW_INCLUDE_DIR_${FFTW_PKG} "${FFTW_INCLUDE_DIR}" CACHE INTERNAL "Path to ${FFTW_HEADER}" FORCE)
endif(FFTW_INCLUDE_DIR)

#we use _${FFTW_PKG} variables to have different cache entries for fftw3 and ffw3f
find_path(FFTW_INCLUDE_DIR_${FFTW_PKG} "${FFTW_HEADER}" HINTS ${PC_FFTW_${FFTW_PKG}_INCLUDE_DIRS})
find_library(FFTW_LIBRARY_${FFTW_PKG} NAMES "${FFTW_PKG}" HINTS ${PC_FFTW_${FFTW_PKG}_LIBRARY_DIRS})

#make _${FFTW_PKG} variables INTERNAL to avoid confusion in cmake-gui
set(FFTW_LIBRARY_${FFTW_PKG} ${FFTW_LIBRARY_${FFTW_PKG}} CACHE INTERNAL "Path to ${FFTW_PKG} library" FORCE)
set(FFTW_INCLUDE_DIR_${FFTW_PKG} ${FFTW_INCLUDE_DIR_${FFTW_PKG}} CACHE INTERNAL "Path to ${FFTW_HEADER}" FORCE)

set(FFTW_LIBRARY "${FFTW_LIBRARY_${FFTW_PKG}}" CACHE FILEPATH "Path to ${FFTW_PKG} library" FORCE)
set(FFTW_INCLUDE_DIR "${FFTW_INCLUDE_DIR_${FFTW_PKG}}" CACHE DIRECTORY "Path to ${FFTW_HEADER}" FORCE)

# set default find_package outcome variables
set(FFTW_LIBRARIES "${FFTW_LIBRARY_${FFTW_PKG}}")
set(FFTW_INCLUDE_DIRS "${FFTW_INCLUDE_DIR_${FFTW_PKG}}")

set(FFTW_FOUND FALSE)
if (FFTW_LIBRARY_${FFTW_PKG} AND FFTW_INCLUDE_DIR_${FFTW_PKG})
  set(FFTW_FOUND TRUE)
elseif (NOT FFTW_LIBRARY_${FFTW_PKG})
  message("Could not find ${FFTW_PKG} library named lib${FFTW_PKG}, please specify its location in FFTW_LIBRARY by hand (e.g. -DFFTW_LIBRARY='/path/to/lib${FFTW_PKG}.so')")
elseif (NOT FFTW_INCLUDE_DIR_${FFTW_PKG})
  message("Could not the ${FFTW_PKG} header ${FFTW_HEADER}, please specify its path in FFTW_INCLUDE_DIR by hand (e.g. -DFFTW_INCLUDE_DIR='/path/to/include')")
endif()

set(FFTW_HAVE_SSE FALSE CACHE BOOL "If ${FFTW_PKG} was built with SSE support")
if (FFTW_FOUND AND HAVE_LIBM AND NOT "${FFTW_LIBRARY_PREVIOUS}" STREQUAL "${FFTW_LIBRARY}")
  #The user could specify trash in FFTW_LIBRARY, so test if we can link it
  include(CheckLibraryExists)
  #adding MATH_LIBRARIES here to allow static libs, this does not harm us as we are anyway using it
  set(CMAKE_REQUIRED_LIBRARIES m)
  unset(FOUND_FFTW_PLAN_${FFTW_PKG} CACHE)
  check_library_exists("${FFTW_LIBRARIES}" "${FFTW_FUNCTION}" "" FOUND_FFTW_PLAN_${FFTW_PKG})
  if(NOT FOUND_FFTW_PLAN_${FFTW_PKG})
    message(FATAL_ERROR "Could not find ${FFTW_FUNCTION} in ${FFTW_LIBRARY_${FFTW_PKG}}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what went wrong. If you are using a static lib (.a) make sure you have specified all dependencies of ${FFTW_PKG} in FFTW_LIBRARY by hand (e.g. -DFFTW_LIBRARY='/path/to/lib${FFTW_PKG}.so;/path/to/libm.so') !")
  endif(NOT FOUND_FFTW_PLAN_${FFTW_PKG})
  #in 3.3 sse function name has changed
  foreach(SSE_FCT fftwf_have_simd_sse2;fftw_have_simd_sse2;fftwf_have_sse;fftw_have_sse2)
    unset(FFTW_HAVE_${SSE_FCT}_${FFTW_PKG} CACHE)
    check_library_exists("${FFTW_LIBRARIES}" "${SSE_FCT}" "" FFTW_HAVE_${SSE_FCT}_${FFTW_PKG})
    if(FFTW_HAVE_${SSE_FCT}_${FFTW_PKG})
      set(FFTW_HAVE_SSE_${FFTW_PKG} "${FFTW_HAVE_${SSE_FCT}_${FFTW_PKG}}")
      break()
    endif(FFTW_HAVE_${SSE_FCT}_${FFTW_PKG})
  endforeach()
  set(FFTW_HAVE_SSE "${FFTW_HAVE_SSE_${FFTW_PKG}}" CACHE BOOL "If ${FFTW_PKG} was built with SSE support" FORCE)
  set(CMAKE_REQUIRED_LIBRARIES)
  set(FFTW_LIBRARY_PREVIOUS "${FFTW_LIBRARY}" CACHE INTERNAL "Value of FFTW_LIBRARY when executing the linking check that last time" FORCE)
endif (FFTW_FOUND AND HAVE_LIBM AND NOT "${FFTW_LIBRARY_PREVIOUS}" STREQUAL "${FFTW_LIBRARY}")

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY FFTW_HAVE_SSE)
