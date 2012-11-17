#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.

#  To help us fund GROMACS development, we humbly ask that you cite
#  the research papers on the package. Check out http://www.gromacs.org
#
#
# - Find FFTW 3
# Find the native FFTW headers and libraries.
#
# Usage:
#
#  find_package(FFTW COMPONENTS [fftw|fftwf])
#
#  fftw or fftwf selects double or single precision FFTW3, respectively
#
# Input:
#
#  FFTW_FIND_COMPONENTS - contains fftw or fftwf component to find (set by the find_package() command)
#  FFTW_PATH            - optional user-specified path prefix in which the FFTW include and lib directories reside
#
# Outputs:
#
#  ${FFTW}_INCLUDE_DIRS - where to find FFTW headers (cache, advanced)
#  ${FFTW}_LIBRARIES    - List of libraries when using FFTW (cache, advanced)
#  ${FFTW}_PKG          - The name of a suitable pkg-config dependency of GROMACS
#  ${FFTW}_HAVE_SIMD    - True if FFTW was build with SIMD support
#  ${FFTW}_FOUND        - True if FFTW was found
#
#
# Note that user-visible (i.e. cache) variables need to be qualified with respect to precision
# so that the user can change the precision subsequently and that will trigger the
# appropriate detection. Internal variables do not need to be qualified, because we never treat
# more than one precision at a time.
#
#  This file is part of GROMACS        Copyright (c) 2012

#
# Is the component valid for GROMACS?
#
list(LENGTH FFTW_FIND_COMPONENTS FFTW_NUM_COMPONENTS_WANTED)
if(${FFTW_NUM_COMPONENTS_WANTED} LESS 1)
  message(FATAL_ERROR "find_package(FFTW) requires you specify a component: fftw or fftwf.")
elseif(${FFTW_NUM_COMPONENTS_WANTED} GREATER 1)
  message(FATAL_ERROR "find_package(FFTW) can only find one FFTW component at a time.")
elseif(${FFTW_FIND_COMPONENTS} MATCHES "^fftw(f)?$")
  if (NOT FFTW_FIND_VERSION OR FFTW_FIND_VERSION EQUAL 3) # Find FFTW3 by default
    string(TOUPPER "${FFTW_FIND_COMPONENTS}" FFTW)
    string(REGEX REPLACE "fftw" "fftw3" FFTW_COMPONENT "${FFTW_FIND_COMPONENTS}")
    set(${FFTW}_FUNCTION_PREFIX "${FFTW_FIND_COMPONENTS}")
  else()
    message(FATAL_ERROR "find_package(FFTW) can only find FFTW version 3.")
  endif()
else()
  message(FATAL_ERROR "find_package(FFTW) does not support finding component ${FFTW_FIND_COMPONENTS}. Choose fftw or fftwf.")
endif()

#
# If we've got pkgconfig and we haven't already got a location for FFTW, try pkgconfig for clues
#
find_package(PkgConfig)
if(NOT __pkg_config_checked_PKGCONFIG_${FFTW} OR NOT ${FFTW}_LIBRARY)
  pkg_check_modules(PKGCONFIG_${FFTW} "${FFTW_COMPONENT}")
endif(NOT __pkg_config_checked_PKGCONFIG_${FFTW} OR NOT ${FFTW}_LIBRARY)

#
# Find path and library
#
find_path(${FFTW}_INCLUDE_DIR "fftw3.h" HINTS ${FFTW_PATH} ${PKGCONFIG_${FFTW}_INCLUDE_DIRS} DOC "Path to fftw3.h include file")
find_library(${FFTW}_LIBRARY NAMES "${FFTW_COMPONENT}" HINTS ${${FFTW}_PATH} ${PKGCONFIG_${FFTW}_LIBRARY_DIRS} PATH_SUFFICES lib lib64 DOC "Library for linking FFTW3")

#
# Set the variables we'll export to the parent scope (and we'd add any dependencies
# here if they existed).
#
set(${FFTW}_LIBRARIES "${${FFTW}_LIBRARY}")
set(${FFTW}_INCLUDE_DIRS "${${FFTW}_INCLUDE_DIR}")

#
# We could use find_package_handle_standard_args() here, but our error
# messages are better for our users.
#
if (${FFTW}_LIBRARY AND ${FFTW}_INCLUDE_DIR)
  set(${FFTW}_FOUND TRUE)
else()
  if (NOT ${FFTW}_LIBRARY)
    message("Could not find ${FFTW_COMPONENT} library named lib${FFTW_COMPONENT}. If it is installed, please specify its path prefix on the command line (e.g. cmake .. -DFFTW_PATH=/usr/local, or in CMAKE_PREFIX_PATH.")
  endif()
  if (NOT ${FFTW}_INCLUDE_DIR)
    message("Could not the ${FFTW_COMPONENT} header fftw3.h. If it is installed, please specify its path prefix on the command line (e.g. cmake .. -DFFTW_PATH=/usr/local, or in CMAKE_PREFIX_PATH.")
  endif()
endif()

#
# Are we satisfied with the FFTW we have found?
#
if (${FFTW}_FOUND)
  # The user could specify trash in ${FFTW}_LIBRARY or ${FFTW}_INCLUDE_PATH, and
  # find_xxx() won't over-ride that, so test if we can compile and link properly.
  if (HAVE_LIBM)
    # Adding MATH_LIBRARIES here to allow static libs, which we'll link later for our own use, anyway
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif (HAVE_LIBM)

  # Does the header file exist?
  set(CMAKE_REQUIRED_INCLUDES ${${FFTW}_INCLUDE_DIR})
  check_include_files(fftw3.h FOUND_FFTW3_H)
  if(NOT FOUND_FFTW3_H)
    message(FATAL_ERROR "FFTW header file ${${FFTW}_INCLUDE_DIR}/fftw3.h does not exist, so GROMACS will not able to link to it. Please specify a valid location where FFTW can be found.")
  endif()

  # Does the library define a relevant function?
  set(CMAKE_REQUIRED_LIBRARIES "${${FFTW}_LIBRARIES}")
  include(CheckSymbolExists)
  check_symbol_exists("${${FFTW}_FUNCTION_PREFIX}_plan_r2r_1d" fftw3.h FOUND_FFTW_PLAN)
  unset(CMAKE_REQUIRED_INCLUDES)
  unset(CMAKE_REQUIRED_LIBRARIES)

  if(NOT FOUND_FFTW_PLAN)
    message(FATAL_ERROR "FFTW function ${${FFTW}_FUNCTION_PREFIX}_plan_r2r_1d in ${${FFTW}_LIBRARY} did not exist in a test program that used the header file ${${FFTW}_INCLUDE_DIRS}/fftw3.h, so GROMACS will not be able to use this library. Take a look at the error message in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log to find out what went wrong. If you are using a static lib (.a), make sure you have specified all dependencies of ${FFTW_COMPONENT} in FFTW_LIBRARY by hand (e.g. cmake .. -DFFTW_LIBRARY='/path/to/lib${FFTW_COMPONENT}.so;/path/to/libm.so') !")
  endif(NOT FOUND_FFTW_PLAN)

  # In FFTW 3.3, the SSE function name has changed
  foreach(SIMD_FUNCTION ${${FFTW}_FUNCTION_PREFIX}_have_simd_sse2;${${FFTW}_FUNCTION_PREFIX}_have_simd_avx;${${FFTW}_FUNCTION_PREFIX}_have_simd_altivec;${${FFTW}_FUNCTION_PREFIX}_have_simd_neon;${${FFTW}_FUNCTION_PREFIX}_have_sse2;${${FFTW}_FUNCTION_PREFIX}_have_sse;${${FFTW}_FUNCTION_PREFIX}_have_altivec)
    include(CheckLibraryExists)
    check_library_exists("${${FFTW}_LIBRARIES}" "${SIMD_FUNCTION}" "" ${FFTW}_HAVE_${SIMD_FUNCTION})
    if(${FFTW}_HAVE_${SIMD_FUNCTION})
      set(${FFTW}_HAVE_SIMD TRUE BOOL "If ${FFTW_COMPONENT} was built with SIMD support")
      break()
    endif(${FFTW}_HAVE_${SIMD_FUNCTION})
  endforeach()

  # Verify that FFTW is compiled with -fPIC (sometimes necessary for GROMACS to link shared libraries)
  if (CMAKE_OBJDUMP AND CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND BUILD_SHARED_LIBS)
      execute_process(COMMAND ${CMAKE_OBJDUMP} --reloc ${${FFTW}_LIBRARY} OUTPUT_VARIABLE FFTW_OBJDUMP)
      if (${FFTW_OBJDUMP} MATCHES "R_X86_64" # Should always be true for static libraries. Checks that objdump works properly and that the library isn't dynamic
              AND NOT ${FFTW_OBJDUMP} MATCHES "R_X86_64_PLT32")
          message(FATAL_ERROR "The FFTW library ${${FFTW}_LIBRARY} cannot be used to link shared libraries in GROMACS. Either provide a different FFTW library by setting ${FFTW}_LIBRARY, recompile FFTW with \"--enable-shared\" or \"--with-pic\", or disable the building of shared libraries in GROMACS by setting BUILD_SHARED_LIBS to \"OFF\". Note: Disabling shared libraries can require up to ten times as much disk space when installing GROMACS.")
      endif()
  endif()
  set(CMAKE_REQUIRED_LIBRARIES)

endif(${FFTW}_FOUND)

set(${FFTW}_HAVE_SIMD FALSE BOOL "Whether FFTW was built with SIMD support")
set(${FFTW}_PKG ${FFTW_COMPONENT})
mark_as_advanced(${FFTW}_INCLUDE_DIR ${FFTW}_LIBRARY)
