# - Find FFTW 3 single
# Find the native FFTWF headers and libraries.
#
#  FFTWF_INCLUDE_DIRS - where to find FFTWF headers
#  FFTWF_LIBRARIES    - List of libraries when using FFTWF.
#  FFTWF_HAVE_SSE     - True if FFTWF was build with SSE support
#  FFTWF_FOUND        - True if FFTWF was found
#
#  This file is part of Gromacs        Copyright (c) 2012

#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.

#  To help us fund GROMACS development, we humbly ask that you cite
#  the research papers on the package. Check out http://www.gromacs.org

find_package(PkgConfig)
if(NOT __pkg_config_checked_PC_FFTWF)
  pkg_check_modules(PC_FFTWF "fftw3f")
endif(NOT __pkg_config_checked_PC_FFTWF)

find_path(FFTWF_INCLUDE_DIR "fftw3.h" HINTS ${PC_FFTWF_INCLUDE_DIRS})
find_library(FFTWF_LIBRARY NAMES "fftw3f" HINTS ${PC_FFTWF_LIBRARY_DIRS})

set(FFTWF_LIBRARIES "${FFTWF_LIBRARY}")
set(FFTWF_INCLUDE_DIRS "${FFTWF_INCLUDE_DIR}")

#better error message than find_package_handle_standard_args
if (FFTWF_LIBRARY AND FFTWF_INCLUDE_DIR)
  set(FFTWF_FOUND TRUE)
elseif (NOT FFTWF_LIBRARY)
  message("Could not find fftw3f library named libfftw3f, please specify its location in FFTWF_LIBRARY by hand (e.g. -DFFTWF_LIBRARY='/path/to/libfftw3f.so')")
elseif (NOT FFTWF_INCLUDE_DIR)
  message("Could not the fftw3f header fftw3.h, please specify its path in FFTWF_INCLUDE_DIR by hand (e.g. -DFFTWF_INCLUDE_DIR='/path/to/include')")
endif()

if (FFTWF_FOUND AND HAVE_LIBM)
  #The user could specify trash in FFTWF_LIBRARY, so test if we can link it
  include(CheckLibraryExists)
  #adding MATH_LIBRARIES here to allow static libs, this does not harm us as we are anyway using it
  set(CMAKE_REQUIRED_LIBRARIES m)
  check_library_exists("${FFTWF_LIBRARIES}" "fftwf_plan_r2r_1d" "" FOUND_FFTWF_PLAN)
  if(NOT FOUND_FFTWF_PLAN)
    message(FATAL_ERROR "Could not find fftwf_plan_r2r_1d in ${FFTWF_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what went wrong. If you are using a static lib (.a) make sure you have specified all dependencies of fftw3f in FFTWF_LIBRARY by hand (e.g. -DFFTWF_LIBRARY='/path/to/libfftw3f.so;/path/to/libm.so') !")
  endif(NOT FOUND_FFTWF_PLAN)
  #in 3.3 sse function name has changed
  foreach(SSE_FCT fftwf_have_simd_sse2;fftwf_have_sse2)
    check_library_exists("${FFTWF_LIBRARIES}" "${SSE_FCT}" "" FFTWF_HAVE_${SSE_FCT})
    if(FFTWF_HAVE_${SSE_FCT})
      set(FFTWF_HAVE_SSE TRUE BOOL "If fftw3f was built with SSE support") 
      break()
    endif(FFTWF_HAVE_${SSE_FCT})
  endforeach()
  set(CMAKE_REQUIRED_LIBRARIES)
endif (FFTWF_FOUND AND HAVE_LIBM)
set(FFTWF_HAVE_SSE FALSE CACHE BOOL "If fftw3f was built with SSE support")

mark_as_advanced(FFTWF_INCLUDE_DIR FFTWF_LIBRARY FFTWF_HAVE_SSE)
