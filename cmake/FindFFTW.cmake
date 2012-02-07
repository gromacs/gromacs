# - Find fftw 2/3 single double
# Find the native FFTW headers and libraries.
#
#  FFTW_INCLUDE_DIRS - where to find fftw headers
#  FFTW_LIBRARIES    - List of libraries when using fftw3.
#  FFTW_PKG          - The name of the pkg-config package needed
#  FFTW_FOUND        - True if fftw
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
  message(FATAL_ERROR "No fftw component to search given")
elseif(${FFTW_NUM_COMPONENTS_WANTED} GREATER 1)
  message(FATAL_ERROR "We only support finding one fftw component at the time, go and implement it ;-)")
elseif(${FFTW_FIND_COMPONENTS} MATCHES "^fftw(f)?$")
  if (${FFTW_FIND_VERSION} EQUAL 3 OR NOT FFTW_FIND_VERSION) #find fftw3 by default
    string(REGEX REPLACE "fftw" "fftw3" FFTW_PKG "${FFTW_FIND_COMPONENTS}")
    set(FFTW_HEADER "fftw3.h")
    set(FFTW_FUNCTION "${FFTW_FIND_COMPONENTS}_plan_r2r_1d")
  #elseif(${FFTW_FIND_VERSION} EQUAL 2)
  #  set(FFTW_PKG "${FFTW_FIND_COMPONENTS}")
  #  set(FFTW_HEADER "${FFTW_FIND_COMPONENTS}.h")
  else()
    message(FATAL_ERROR "We only support finding fftw version 3, go and implement it ;-)")
  endif()
else()
  message(FATAL_ERROR "We do not support finding ${FFTW_FIND_COMPONENTS}, go and implement it ;-)")
endif()

pkg_check_modules(PC_FFTW "${FFTW_PKG}")
find_path(FFTW_INCLUDE_DIR_${FFTW_PKG} "${FFTW_HEADER}" HINTS ${PC_FFTW_INCLUDE_DIRS})

find_library(FFTW_LIBRARY_${FFTW_PKG} NAMES "${FFTW_PKG}" HINTS ${PC_FFTW_LIBRARY_DIRS} )

#do that so that the user can use FFTW_LIBRARY, FFTW_INCLUDE_DIR to provide custom stuff 
set(FFTW_LIBRARY ${FFTW_LIBRARY_${FFTW_PKG}})
set(FFTW_INCLUDE_DIR ${FFTW_INCLUDE_DIR_${FFTW_PKG}})

set(FFTW_LIBRARIES ${FFTW_LIBRARY} )
set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_INCLUDE_DIR )

if (FFTW_FOUND AND HAVE_LIBM)
  include(CheckLibraryExists)
  #adding MATH_LIBRARIES here to allow static libs, this does not harm us as we are anyway using it
  check_library_exists("${FFTW_LIBRARIES};m" "${FFTW_FUNCTION}" "" FOUND_FFTW_PLAN)
  if(NOT FOUND_FFTW_PLAN)
    message(FATAL_ERROR "Could not find ${FFTW_FUNCTION} in ${FFTW_LIBRARY}, take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what was going wrong. If you are using a static lib (.a) make sure you have specified all dependencies of ${FFTW_PKG} in FFTW_LIBRARY by hand (i.e. -DFFTW_LIBRARY='/path/to/libfftw3.so;/path/to/libm.so') !")
  endif(NOT FOUND_FFTW_PLAN)
endif (FFTW_FOUND AND HAVE_LIBM)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY )
