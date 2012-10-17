# - Find FFTW 3 double
# Find the native FFTW headers and libraries.
#
#  FFTW_INCLUDE_DIRS - where to find FFTW headers
#  FFTW_LIBRARIES    - List of libraries when using FFTW.
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
if(NOT __pkg_config_checked_PC_FFTW)
  pkg_check_modules(PC_FFTW "fftw3")
endif(NOT __pkg_config_checked_PC_FFTW)

find_path(FFTW_INCLUDE_DIR "fftw3.h" HINTS ${PC_FFTW_INCLUDE_DIRS})
find_library(FFTW_LIBRARY NAMES "fftw3" HINTS ${PC_FFTW_LIBRARY_DIRS})

set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
set(FFTW_INCLUDE_DIRS "${FFTW_INCLUDE_DIR}")

#better error message than find_package_handle_standard_args
if (FFTW_LIBRARY AND FFTW_INCLUDE_DIR)
  set(FFTW_FOUND TRUE)
elseif (NOT FFTW_LIBRARY)
  message("Could not find fftw3 library named libfftw3, please specify its location in FFTW_LIBRARY by hand (e.g. -DFFTW_LIBRARY='/path/to/libfftw3.so')")
elseif (NOT FFTW_INCLUDE_DIR)
  message("Could not the fftw3 header fftw3.h, please specify its path in FFTW_INCLUDE_DIR by hand (e.g. -DFFTW_INCLUDE_DIR='/path/to/include')")
endif()

set(FFTW_HAVE_SSE FALSE CACHE BOOL "If fftw3 was built with SSE support")
if (FFTW_FOUND AND HAVE_LIBM)
  #The user could specify trash in FFTW_LIBRARY, so test if we can link it
  include(CheckLibraryExists)
  #adding MATH_LIBRARIES here to allow static libs, this does not harm us as we are anyway using it
  set(CMAKE_REQUIRED_LIBRARIES m)
  check_library_exists("${FFTW_LIBRARIES}" "fftw_plan_r2r_1d" "" FOUND_FFTW_PLAN)
  if(NOT FOUND_FFTW_PLAN)
    message(FATAL_ERROR "Could not find fftw_plan_r2r_1d in ${FFTW_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what went wrong. If you are using a static lib (.a) make sure you have specified all dependencies of fftw3 in FFTW_LIBRARY by hand (e.g. -DFFTW_LIBRARY='/path/to/libfftw3.so;/path/to/libm.so') !")
  endif(NOT FOUND_FFTW_PLAN)
  #in 3.3 sse function name has changed
  foreach(SSE_FCT fftw_have_simd_sse2;fftw_have_sse2)
    check_library_exists("${FFTW_LIBRARIES}" "${SSE_FCT}" "" FFTW_HAVE_${SSE_FCT})
    if(FFTW_HAVE_${SSE_FCT})
      set(FFTW_HAVE_SSE "${FFTW_HAVE_${SSE_FCT}}")
      break()
    endif(FFTW_HAVE_${SSE_FCT})
  endforeach()
  set(CMAKE_REQUIRED_LIBRARIES)
endif (FFTW_FOUND AND HAVE_LIBM)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY FFTW_HAVE_SSE)
